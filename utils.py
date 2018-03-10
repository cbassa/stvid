from astropy.coordinates import SkyCoord,EarthLocation,AltAz,FK5,get_sun
from astropy.time import Time
import astropy.units as u
import numpy as np
from scipy import interpolate

# Sunrise/sunset algorithm from Astronomical Algorithms by Jean Meeus
def get_sunset_and_sunrise(tnow,loc,refalt):
    # Get time
    nmjd=64
    mjd0=np.floor(tnow.mjd)
    mnow=tnow.mjd-mjd0
    mjd=np.linspace(mjd0-1.0,mjd0+2.0,nmjd)
    t=Time(mjd,format='mjd',scale='utc')
    
    # Get sun position
    psun=get_sun(t)
    
    # Correct for precession by converting to FK5
    pos=SkyCoord(ra=psun.ra.degree,dec=psun.dec.degree,frame='icrs',unit='deg').transform_to(FK5(equinox=t))
    
    # Compute altitude extrema
    de=np.mean(pos.dec)
    minalt=np.arcsin(np.sin(loc.latitude)*np.sin(de)-np.cos(loc.latitude)*np.cos(de))
    maxalt=np.arcsin(np.sin(loc.latitude)*np.sin(de)+np.cos(loc.latitude)*np.cos(de))
   
    # Never sets, never rises?
    if minalt>refalt:
        return "sun never sets",t[0],t[0]
    elif maxalt<refalt:
        return "sun never rises",t[0],t[0]

    # Set up interpolating function
    fra=interpolate.interp1d(t.mjd,pos.ra.degree)
    fde=interpolate.interp1d(t.mjd,pos.dec.degree)

    # Get GMST
    gmst0=Time(mjd0,format='mjd',scale='utc').sidereal_time('mean','greenwich')

    # Get transit time
    mtransit=np.mod((fra(mjd0)*u.degree-loc.longitude-gmst0)/(360.0*u.degree),1.0)
    while True:
        gmst=gmst0+360.985647*u.deg*mtransit
        ra=fra(mjd0+mtransit)*u.deg
        ha=gmst+loc.longitude-ra
        mtransit-=ha/(360.0*u.deg)
        if np.abs(ha.degree)<1e-9:
            break

    # Set transit
    ttransit=Time(mjd0+mtransit.value,format='mjd',scale='utc')

    # Hour angle offset
    ha0=np.arccos((np.sin(refalt)-np.sin(loc.latitude)*np.sin(np.mean(pos.dec)))/(np.cos(loc.latitude)*np.cos(np.mean(pos.dec))))

    # Get rise time
    mrise=mtransit-ha0/(360.0*u.deg)
    if mrise<mnow:
        mrise+=1.0
    while True:
        gmst=gmst0+360.985647*u.deg*mrise
        ra=fra(mjd0+mrise)*u.deg
        de=fde(mjd0+mrise)*u.deg
        ha=gmst+loc.longitude-ra
        alt=np.arcsin(np.sin(loc.latitude)*np.sin(de)+np.cos(loc.latitude)*np.cos(de)*np.cos(ha))
        dm=(alt-refalt)/(360.0*u.deg*np.cos(de)*np.cos(loc.latitude)*np.sin(ha))
        mrise+=dm
        if np.abs(dm)<1e-9:
            break

    # Set rise time
    trise=Time(mjd0+mrise.value,format='mjd',scale='utc')

    # Get set time
    mset=mtransit+ha0/(360.0*u.deg)
    if mset<mnow:
        mset+=1.0
    while True:
        gmst=gmst0+360.985647*u.deg*mset
        ra=fra(mjd0+mset)*u.deg
        de=fde(mjd0+mset)*u.deg
        ha=gmst+loc.longitude-ra
        alt=np.arcsin(np.sin(loc.latitude)*np.sin(de)+np.cos(loc.latitude)*np.cos(de)*np.cos(ha))
        dm=(alt-refalt)/(360.0*u.deg*np.cos(de)*np.cos(loc.latitude)*np.sin(ha))
        mset+=dm
        if np.abs(dm)<1e-9:
            break

    # Set set time
    tset=Time(mjd0+mset.value,format='mjd',scale='utc')

    return "sun rises and sets",tset,trise
