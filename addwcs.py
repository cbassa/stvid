#!/usr/bin/env python
from __future__ import print_function
import os
import subprocess
import numpy as np
from astropy.io import fits
from astropy import wcs
import astropy.units as u
import warnings
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import SkyCoord,FK5,ICRS
from astropy.time import Time
import ppgplot
from scipy import optimize
import glob

class tycho2_catalog:
    """Tycho2 catalog"""

    def __init__(self,maxmag=9.0):
        hdu=fits.open("/home/bassa/code/python/satellite/astrometry-tycho2/build/tyc2.fits")

        ra=hdu[1].data.field('RA')*u.deg
        dec=hdu[1].data.field('DEC')*u.deg
        mag=hdu[1].data.field('MAG_VT')

        c=mag<maxmag
        self.ra=ra[c]
        self.dec=dec[c]
        self.mag=mag[c]
        
class pixel_catalog:
    """Pixel catalog"""

    def __init__(self,fname):
        d=np.loadtxt(fname)
        self.x=d[:,0]
        self.y=d[:,1]
        self.mag=d[:,2]
        self.ra=np.empty_like(self.x)
        self.dec=np.empty_like(self.x)
        self.imag=np.empty_like(self.x)
        self.flag=np.zeros_like(self.x)

def estimate_wcs_from_reference(ref,fname):
    # Read header of reference
    hdu=fits.open(ref)
    hdu[0].header["NAXIS"]=2
    w=wcs.WCS(hdu[0].header)

    # Get time and position from reference
    tref=Time(hdu[0].header["MJD-OBS"],format="mjd",scale="utc")
    pref=SkyCoord(ra=w.wcs.crval[0],dec=w.wcs.crval[1],unit="deg",frame="icrs").transform_to(FK5(equinox=tref))

    # Read time from target
    hdu=fits.open(fname)
    t=Time(hdu[0].header["MJD-OBS"],format="mjd",scale="utc")

    # Correct wcs
    dra=t.sidereal_time("mean","greenwich")-tref.sidereal_time("mean","greenwich")
    p=FK5(ra=pref.ra+dra,dec=pref.dec,equinox=t).transform_to(ICRS)
    w.wcs.crval=np.array([p.ra.degree,p.dec.degree])

    return w

def match_catalogs(ast_catalog,pix_catalog,w,rmin):
    # Select stars towards pointing center
    ra,dec=w.wcs.crval*u.deg
    d=np.arccos(np.sin(dec)*np.sin(ast_catalog.dec)+np.cos(dec)*np.cos(ast_catalog.dec)*np.cos(ra-ast_catalog.ra))
    c=(d<30.0*u.deg)
    ra,dec,mag=ast_catalog.ra[c],ast_catalog.dec[c],ast_catalog.mag[c]

    # Convert RA/Dec to pixels
    pix=w.wcs_world2pix(np.stack((ra,dec),axis=-1),0)
    xs,ys=pix[:,0],pix[:,1]

    for i in range(len(pix_catalog.x)):
        dx=xs-pix_catalog.x[i]
        dy=ys-pix_catalog.y[i]
        r=np.sqrt(dx*dx+dy*dy)
        if np.min(r)<rmin:
            j=np.argmin(r)
            pix_catalog.ra[i]=ra[j].value
            pix_catalog.dec[i]=dec[j].value
            pix_catalog.imag[i]=mag[j]
            pix_catalog.flag[i]=1

    return xs,ys

def residual(a,x,y,z):
    return z-(a[0]+a[1]*x+a[2]*y)

def fit_wcs(w,pix_catalog):
    x0,y0=w.wcs.crpix
    ra0,dec0=w.wcs.crval
    dx,dy=pix_catalog.x-x0,pix_catalog.y-y0

    # Iterate to remove outliers
    nstars=np.sum(pix_catalog.flag==1)
    for j in range(10):
        w=wcs.WCS(naxis=2)
        w.wcs.crpix=np.array([0.0,0.0])
        w.wcs.cd=np.array([[1.0,0.0],[0.0,1.0]])
        w.wcs.ctype=["RA---TAN","DEC--TAN"]
        w.wcs.set_pv([(2,1,45.0)])
        c=pix_catalog.flag==1

        # Iterate to move crval to crpix location
        for i in range(5):
            w.wcs.crval=np.array([ra0,dec0])
            
            world=np.stack((pix_catalog.ra,pix_catalog.dec),axis=-1)
            pix=w.wcs_world2pix(world,1)
            rx,ry=pix[:,0],pix[:,1]
            
            ax,cov_q,infodict,mesg,ierr=optimize.leastsq(residual,[0.0,0.0,0.0],args=(dx[c],dy[c],rx[c]),full_output=1)
            ay,cov_q,infodict,mesg,ierr=optimize.leastsq(residual,[0.0,0.0,0.0],args=(dx[c],dy[c],ry[c]),full_output=1)
        
            ra0,dec0=w.wcs_pix2world([[ax[0],ay[0]]],1)[0]

        # Compute residuals
        drx=ax[0]+ax[1]*dx+ax[2]*dy-rx
        dry=ay[0]+ay[1]*dx+ay[2]*dy-ry
        dr=np.sqrt(drx*drx+dry*dry)
        rms=np.sqrt(np.sum(dr[c]**2)/len(dr[c]))

        dr[~c]=1.0
        c=(dr<2.0*rms)
        pix_catalog.flag[~c]=0

        # Break if converged
        if np.sum(c)==nstars:
            break
        nstars=np.sum(c)

    # Compute residuals
    rmsx=np.sqrt(np.sum(drx[c]**2)/len(drx[c]))
    rmsy=np.sqrt(np.sum(dry[c]**2)/len(dry[c]))

    # Store header
    w=wcs.WCS(naxis=2)
    w.wcs.crpix=np.array([x0,y0])
    w.wcs.crval=np.array([ra0,dec0])
    w.wcs.cd=np.array([[ax[1],ax[2]],[ay[1],ay[2]]])
    w.wcs.ctype=["RA---TAN","DEC--TAN"]
    w.wcs.set_pv([(2,1,45.0)])

    whdr={"CRPIX1":x0,"CRPIX2":y0,"CRVAL1":ra0,"CRVAL2":dec0,"CD1_1":ax[1],"CD1_2":ax[2],"CD2_1":ay[1],"CD2_2":ay[2],"CTYPE1":"RA---TAN","CTYPE2":"DEC--TAN","CUNIT1":"DEG","CUNIT2":"DEG","CRRES1":rmsx,"CRRES2":rmsy}    

    return w,whdr,rmsx,rmsy,rms

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

if __name__ == "__main__":
    # Set warnings
    warnings.filterwarnings("ignore", category=UserWarning, append=True)
    warnings.simplefilter("ignore", AstropyWarning)

    # Files
    ref="2018-03-14T19:00:03.086.fits"

    files=sorted(glob.glob("2*.fits"))
#    fname="2018-03-14T19:59:54.552.fits"
#    fname="2018-03-14T19:00:03.086.fits"
#    fname="../20180225a/tmp/2018-02-25T20:29:55.944.fits"

    for fname in files[0:10]:
        # Generate star catalog
        if not os.path.exists(fname+".cat"):
            generate_star_catalog(fname)

        # Estimated WCS
        w=estimate_wcs_from_reference(ref,fname)

        # Read catalogs
        ast_catalog=tycho2_catalog(10.0)
        pix_catalog=pixel_catalog(fname+".cat")

        # Match catalogs
        x,y=match_catalogs(ast_catalog,pix_catalog,w,10.0)

        # Fit transformation
        w,whdr,rmsx,rmsy,rms=fit_wcs(w,pix_catalog)

        nused=np.sum(pix_catalog.flag==1)
        nstars=len(pix_catalog.flag)

        print("%s %8.4f %8.4f %3d/%3d %6.1f %6.1f %6.1f"%(fname,w.wcs.crval[0],w.wcs.crval[1],nused,nstars,3600.0*rmsx,3600.0*rmsy,3600.0*np.sqrt(rmsx**2+rmsy**2)))

    
#        hdu=fits.open(fname)
#        img=hdu[0].data[0]
#        ny,nx=img.shape

    
#        ppgplot.pgopen("/xs")
#        ppgplot.pgwnad(0.0,nx,0.0,ny)

#        vmin=np.mean(img)-2.0*np.std(img)
#        vmax=np.mean(img)+3.0*np.std(img)
    
#        ppgplot.pgimag(img,nx,ny,0,nx-1,0,ny-1,vmax,vmin,np.array([-0.5,1.0,0.0,-0.5,0.0,1.0]))

#        ppgplot.pgbox("BCTSNI",0.,0,"BCTSNI",0.,0)

#        ppgplot.pgsci(3)
#        ppgplot.pgpt(pix_catalog.x,pix_catalog.y,4)
#        ppgplot.pgsci(4)
#        ppgplot.pgpt(x,y,24)
#        ppgplot.pgsci(2)
#        ppgplot.pgpt(pix_catalog.x[c],pix_catalog.y[c],6)
#        ppgplot.pgend()
    
