#!/usr/bin/env python
from __future__ import print_function
import os,sys
import glob
import subprocess
from stio import fourframe
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import EarthLocation,AltAz,GCRS,FK5
from astropy.time import Time
import astropy.units as u
from scipy import optimize
import warnings
from astropy.utils.exceptions import AstropyWarning
import ppgplot

def residual(a,x,y,z):
    return z-(a[0]+a[1]*x+a[2]*y)

def offsets(ra0,de0,ra,de):
    w=wcs.WCS(naxis=2)
    w.wcs.crpix=np.array([0.0,0.0])
    w.wcs.crval=u.Quantity([ra0,de0])
    w.wcs.cd=np.array([[1.0/3600.0,0.0],[0.0,1.0/3600.0]])
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    w.wcs.set_pv([(2,1,45.0)])

    world=np.stack((ra,de),axis=-1)
    pix=w.wcs_world2pix(world,1)
    return pix[:,0],pix[:,1]
    

def read_pixel_catalog(fname):
    d=np.loadtxt(fname)
    x,y,mag=d[:,0],d[:,1],d[:,2]
    
    return x,y,mag

def read_tycho2_catalog(maxmag=9.0):
    hdu=fits.open("/home/bassa/code/python/satellite/astrometry-tycho2/build/tyc2.fits")

    ra=hdu[1].data.field('RA')*u.deg
    de=hdu[1].data.field('DEC')*u.deg
    mag=hdu[1].data.field('MAG_VT')

    c=mag<maxmag
    print("Selected %d stars from Tycho-2 catalog with m<%.1f"%(np.sum(c),maxmag))
    return ra[c],de[c],mag[c]

def match_catalogs(ra,de,x,y,xs,ys,rmin):
    res=[]
    for i in range(len(x)):
        dx,dy=xs-x[i],ys-y[i]
        r=np.sqrt(dx*dx+dy*dy)
        if np.min(r)<rmin:
            j=np.argmin(r)
            res.append([ra[i],de[i],xs[j],ys[j]])
            
    return np.array(res)

def fit_wcs(res,x0,y0):
    ra,de,x,y=res[:,0],res[:,1],res[:,2],res[:,3]
    ra0,de0=np.mean(ra),np.mean(de)
    dx,dy=x-x0,y-y0

    for i in range(5):
        w=wcs.WCS(naxis=2)
        w.wcs.crpix=np.array([0.0,0.0])
        w.wcs.crval=np.array([ra0,de0])
        w.wcs.cd=np.array([[1.0,0.0],[0.0,1.0]])
        w.wcs.ctype=["RA---TAN","DEC--TAN"]
        w.wcs.set_pv([(2,1,45.0)])

        world=np.stack((ra,de),axis=-1)
        pix=w.wcs_world2pix(world,1)
        rx,ry=pix[:,0],pix[:,1]

        ax,cov_q,infodict,mesg,ierr=optimize.leastsq(residual,[0.0,0.0,0.0],args=(dx,dy,rx),full_output=1)
        ay,cov_q,infodict,mesg,ierr=optimize.leastsq(residual,[0.0,0.0,0.0],args=(dx,dy,ry),full_output=1)

        ra0,de0=w.wcs_pix2world([[ax[0],ay[0]]],1)[0]

    w=wcs.WCS(naxis=2)
    w.wcs.crpix=np.array([x0,y0])
    w.wcs.crval=np.array([ra0,de0])
    w.wcs.cd=np.array([[ax[1],ax[2]],[ay[1],ay[2]]])
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    w.wcs.set_pv([(2,1,45.0)])

    whdr={"CRPIX1":x0,"CRPIX2":y0,"CRVAL1":ra0,"CRVAL2":de0,"CD1_1":ax[1],"CD1_2":ax[2],"CD2_1":ay[1],"CD2_2":ay[2],"CTYPE1":"RA---TAN","CTYPE2":"DEC--TAN","CUNIT1":"DEG","CUNIT2":"DEG"}

    return whdr

def generate_star_catalog(file):
    # Get sextractor location
    env=os.getenv("ST_DATADIR")

    # Format command
    command="sextractor %s -c %s/sextractor/default.sex"%(file,env)

    # Run sextractor
    output=subprocess.check_output(command,shell=True,stderr=subprocess.STDOUT)

    # Rename file
    if os.path.exists("test.cat"):
        os.rename("test.cat",file+".cat")

    return

def estimate_pixel_locations(ff,loc,az,alt,sx,sy):
    # Get RA/Dec/parallactic angle
    t=Time(ff.mjd,format='mjd',scale='utc')

    # Compute position
    p=AltAz(az=az,alt=alt,obstime=t,location=loc)

    # Transform to GCRS and FK5
    pgcrs=p.transform_to(GCRS)
    pfk5=p.transform_to(FK5(equinox=t))

    # Compute parallactic angle
    h=t.sidereal_time('mean','greenwich')+loc.lon-pfk5.ra
    q=np.arctan2(np.sin(h),np.tan(loc.lat)*np.cos(pfk5.dec)-np.sin(pfk5.dec)*np.cos(h))

    # Compute offsets
    d=np.arccos(np.sin(pgcrs.dec)*np.sin(de)+np.cos(pgcrs.dec)*np.cos(de)*np.cos(pgcrs.ra-ra))
    c=(d<30.0*u.deg)
    rx,ry=offsets(pgcrs.ra,pgcrs.dec,ra[c],de[c])

    # Estimate positions
    xp=-(np.cos(q)*rx-np.sin(q)*ry)/sx+ff.nx/2.0
    yp=(np.sin(q)*rx+np.cos(q)*ry)/sy+ff.ny/2.0

    return xp,yp,mag[c]


    
if __name__ == "__main__":
    # Set warnings
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    warnings.simplefilter('ignore', AstropyWarning)

    rmin,rmax=1.0,10.0
    mmin,mmax=5.0,10.0
    
    # ppgplot arrays
    tr=np.array([-0.5,1.0,0.0,-0.5,0.0,1.0])
    tr=np.array([0.,1.0,0.0,0.,0.0,1.0])
    heat_l=np.array([0.0,0.2,0.4,0.6,1.0])
    heat_r=np.array([0.0,0.5,1.0,1.0,1.0])
    heat_g=np.array([0.0,0.0,0.5,1.0,1.0])
    heat_b=np.array([0.0,0.0,0.0,0.3,1.0])
    
    # Get files
    files=sorted(glob.glob("2*.fits"))

    # Read Tycho2 catalog
    ra,de,mag=read_tycho2_catalog(10.0)

    # Set location
    loc=EarthLocation(lat=52.8344*u.deg,lon=6.3785*u.deg,height=10*u.m)

    # Default position
    az,alt=31.21*u.deg,29.60*u.deg
    
    # Loop over files
    for file in files[0:2]:
        # Generate star catalog
        if not os.path.exists(file+".cat"):
            generate_star_catalog(file)

        # Read pixel catalog
        xs,ys,ms=read_pixel_catalog(file+".cat")
        
        # Read fourframe
        ff=fourframe(file)

        xp,yp,mag=estimate_pixel_locations(ff,loc,az,alt,36,33);

        rp=rmax-(rmax-rmin)*(mag-mmin)/(mmax-mmin)
        
        # Plot
        ppgplot.pgopen("/xs")
        ppgplot.pgsfs(2)
        ppgplot.pgpap(0.0,1.0)
        ppgplot.pgsvp(0.1,0.95,0.1,0.8)
            
        ppgplot.pgsch(0.8)
        ppgplot.pgmtxt("T",6.0,0.0,0.0,"UT Date: %.23s  COSPAR ID: %04d"%(ff.nfd,ff.site_id))
        if (3600.0*ff.crres[0]<1e-3) | (3600.0*ff.crres[1]<1e-3) | (ff.crres[0]/ff.sx>2.0) | (ff.crres[1]/ff.sy>2.0):
            ppgplot.pgsci(2)
        else:
            ppgplot.pgsci(1)
        ppgplot.pgmtxt("T",4.8,0.0,0.0,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')"%(ff.crval[0],3600.0*ff.crres[0],ff.crval[1],3600.0*ff.crres[1]))
        ppgplot.pgsci(1)
        ppgplot.pgmtxt("T",3.6,0.0,0.0,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d"%(ff.wx,ff.wy,3600.0*ff.sx,3600.0*ff.sy))
        ppgplot.pgmtxt("T",2.4,0.0,0.0,"Stat: %5.1f+-%.1f (%.1f-%.1f)"%(np.mean(ff.zavg),np.std(ff.zavg),ff.zavgmin,ff.zavgmax))
        
        ppgplot.pgsch(1.0)
        ppgplot.pgwnad(0.0,ff.nx,0.0,ff.ny)
        ppgplot.pglab("x (pix)","y (pix)"," ")
        ppgplot.pgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5)
            
        ppgplot.pgimag(ff.zavg,ff.nx,ff.ny,0,ff.nx-1,0,ff.ny-1,ff.zavgmax,ff.zavgmin,tr)
        ppgplot.pgbox("BCTSNI",0.,0,"BCTSNI",0.,0)
        ppgplot.pgstbg(1)

        ppgplot.pgsci(3)
        ppgplot.pgpt(xs,ys,4)
        ppgplot.pgsci(4)
        for i in range(len(rp)):
            ppgplot.pgcirc(xp[i],yp[i],rp[i])
        #        ppgplot.pgpt(xp,yp,24)


        ppgplot.pgend()
        
            
    sys.exit()
            
    # Read fits
    hdu=fits.open(fname)

    # Fix to force WCS with two axes only
    hdu[0].header["NAXIS"]=2
    w=wcs.WCS(hdu[0].header)

    # Read data
    data=hdu[0].data
    zavg,zstd,zmax,znum=data
    ny,nx=zavg.shape

    # Get extrema
    vmin=np.mean(zavg)-2.0*np.std(zavg)
    vmax=np.mean(zavg)+3.0*np.std(zavg)

    # Read pixel catalog
    xs,ys,ms=read_pixel_catalog(fname+".cat")


    
    # Convert to pixel coordinates
    xc,yc=w.all_world2pix(ra,de,1)
    c=(xc>0) & (xc<nx) & (yc>0) & (yc<ny)
    xt,yt=xc[c],yc[c]
    rat,det=ra[c],de[c]
            
    # Match catalogs
    res=match_catalogs(rat,det,xt,yt,xs,ys,10.0)

    whdr=fit_wcs(res,nx//2,ny//2)
    
    #hdu=fits.PrimaryHDU(header=w.to_header(),data=data)
    #hdr=fits.Header()
    hdr=hdu[0].header
    for k,v in whdr.items():
        hdr[k]=v
    hdu=fits.PrimaryHDU(header=hdr,data=data)
    hdu.writeto("test.fits",overwrite=True)
    
