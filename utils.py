#!/usr/bin/env python
# Common Utilities

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, FK5, get_sun
from astropy.time import Time
import astropy.units as u
import numpy as np
from scipy import interpolate


# Sunrise/sunset algorithm from Astronomical Algorithms by Jean Meeus
def get_sunset_and_sunrise(tnow, loc, refalt):
    # Get time
    nmjd = 64
    mjd0 = np.floor(tnow.mjd)
    mnow = tnow.mjd-mjd0
    mjd = np.linspace(mjd0-1.0, mjd0+3.0, nmjd)
    t = Time(mjd, format='mjd', scale='utc')

    # Get sun position
    psun = get_sun(t)

    # Correct for precession by converting to FK5
    pos = SkyCoord(ra=psun.ra.degree,
                   dec=psun.dec.degree,
                   frame='icrs',
                   unit='deg').transform_to(FK5(equinox=t))

    # Compute altitude extrema
    de = np.mean(pos.dec)
    minalt = np.arcsin(np.sin(loc.lat)*np.sin(de)-np.cos(loc.lat)*np.cos(de))
    maxalt = np.arcsin(np.sin(loc.lat)*np.sin(de)+np.cos(loc.lat)*np.cos(de))

    # Never sets, never rises?
    if minalt > refalt:
        return "sun never sets", t[0], t[0]
    elif maxalt < refalt:
        return "sun never rises", t[0], t[0]

    # Prevent discontinuities in right ascension
    dra = pos[-1].ra-pos[0].ra
    ra = pos.ra.degree
    if dra < 0.0:
        c = pos.ra.degree < 180.0
        ra[c] = pos[c].ra.degree+360.0

    # Set up interpolating function
    fra = interpolate.interp1d(t.mjd, ra)
    fde = interpolate.interp1d(t.mjd, pos.dec.degree)

    # Get GMST
    gmst0 = Time(mjd0,
                 format='mjd',
                 scale='utc').sidereal_time('mean', 'greenwich')

    # Get transit time
    mtransit = np.mod((fra(mjd0)*u.degree-loc.lon-gmst0)/(360.0*u.degree), 1.0)
    while True:
        gmst = gmst0+360.985647*u.deg*mtransit
        ra = fra(mjd0+mtransit)*u.deg
        ha = gmst+loc.lon-ra
        mtransit -= ha/(360.0*u.deg)
        if np.abs(ha.degree) < 1e-9:
            break

    # Hour angle offset
    ha0 = np.arccos((np.sin(refalt)
                     - np.sin(loc.lat)
                     * np.sin(np.mean(pos.dec)))
                    / (np.cos(loc.lat)
                    * np.cos(np.mean(pos.dec))))

    # Get set time
    mset = mtransit+ha0/(360.0*u.deg)
    while True:
        gmst = gmst0+360.985647*u.deg*mset
        ra = fra(mjd0+mset)*u.deg
        de = fde(mjd0+mset)*u.deg
        ha = gmst+loc.lon-ra
        alt = np.arcsin(np.sin(loc.lat)
                        * np.sin(de)
                        + np.cos(loc.lat)
                        * np.cos(de)
                        * np.cos(ha))
        dm = (alt-refalt)/(360.0
                           * u.deg
                           * np.cos(de)
                           * np.cos(loc.lat)
                           * np.sin(ha))
        mset += dm

        # Break loop or find sunset on next day
        if np.abs(dm) < 1e-9:
            if mset >= mnow:
                break
            else:
                mset += 1.0

    # Set set time
    tset = Time(mjd0+mset.value, format='mjd', scale='utc')

    # Get rise time
    mrise = mtransit-ha0/(360.0*u.deg)
    while True:
        gmst = gmst0+360.985647*u.deg*mrise
        ra = fra(mjd0+mrise)*u.deg
        de = fde(mjd0+mrise)*u.deg
        ha = gmst+loc.lon-ra
        alt = np.arcsin(np.sin(loc.lat)
                        * np.sin(de)
                        + np.cos(loc.lat)
                        * np.cos(de)
                        * np.cos(ha))
        dm = (alt-refalt)/(360.0
                           * u.deg
                           * np.cos(de)
                           * np.cos(loc.lat)
                           * np.sin(ha))
        mrise += dm

        # Break loop or find sunrise on next day
        if np.abs(dm) < 1e-9:
            if mrise >= mnow:
                break
            else:
                mrise += 1.0

    # Set rise time
    trise = Time(mjd0+mrise.value, format='mjd', scale='utc')

    return "sun rises and sets", tset, trise

