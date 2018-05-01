#!/usr/bin/env python
from __future__ import print_function
import os
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord,FK5,ICRS 
from astropy.time import Time
from scipy import optimize

# Class for the Tycho 2 catalog
class tycho2_catalog:
    """Tycho2 catalog"""

    def __init__(self,maxmag=9.0):
        hdu=fits.open(os.path.join(os.getenv("ST_DATADIR"),"data/tyc2.fits"))

        ra=hdu[1].data.field('RA')*u.deg
        dec=hdu[1].data.field('DEC')*u.deg
        mag=hdu[1].data.field('MAG_VT')

        c=mag<maxmag
        self.ra=ra[c]
        self.dec=dec[c]
        self.mag=mag[c]

# Estimate the WCS from a reference file
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

# Match the astrometry and pixel catalog
def match_catalogs(ast_catalog,pix_catalog,w,rmin):
    # Select stars towards pointing center
    ra,dec=w.wcs.crval*u.deg
    d=np.arccos(np.sin(dec)*np.sin(ast_catalog.dec)+np.cos(dec)*np.cos(ast_catalog.dec)*np.cos(ra-ast_catalog.ra))
    c=(d<30.0*u.deg)
    ra,dec,mag=ast_catalog.ra[c],ast_catalog.dec[c],ast_catalog.mag[c]

    # Convert RA/Dec to pixels
    pix=w.wcs_world2pix(np.stack((ra,dec),axis=-1),0)
    xs,ys=pix[:,0],pix[:,1]

    # Loop over stars
    nmatch=0
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
            nmatch+=1

    return nmatch

# Residual function
def residual(a,x,y,z):
    return z-(a[0]+a[1]*x+a[2]*y)

# Fit transformation
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

    return w,rmsx,rmsy,rms

def add_wcs(fname,w,rmsx,rmsy):
    # Read fits
    hdu=fits.open(fname)

    whdr={"CRPIX1":w.wcs.crpix[0],"CRPIX2":w.wcs.crpix[1],"CRVAL1":w.wcs.crval[0],"CRVAL2":w.wcs.crval[1],"CD1_1":w.wcs.cd[0,0],"CD1_2":w.wcs.cd[0,1],"CD2_1":w.wcs.cd[1,0],"CD2_2":w.wcs.cd[1,1],"CTYPE1":"RA---TAN","CTYPE2":"DEC--TAN","CUNIT1":"DEG","CUNIT2":"DEG","CRRES1":rmsx,"CRRES2":rmsy}    
    
    # Add keywords
    hdr=hdu[0].header
    for k,v in whdr.items():
        hdr[k]=v
    hdu=fits.PrimaryHDU(header=hdr,data=hdu[0].data)
    hdu.writeto(fname,overwrite=True,output_verify="ignore")

    return
    
def calibrate_from_reference(fname,ref,pix_catalog):
    # Estimated WCS
    w=estimate_wcs_from_reference(ref,fname)

    # Default rms values
    rmsx=0.0
    rmsy=0.0
    
    # Read catalogs
    if (pix_catalog.nstars>4):
        ast_catalog=tycho2_catalog(10.0)
        
        # Match catalogs
        nmatch=match_catalogs(ast_catalog,pix_catalog,w,10.0)
    
        # Fit transformation
        if nmatch>4:
            w,rmsx,rmsy,rms=fit_wcs(w,pix_catalog)
        
    # Add wcs
    add_wcs(fname,w,rmsx,rmsy)

    return w,rmsx,rmsy
