#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord,FK5,AltAz,EarthLocation
from astropy.time import Time

table=ascii.read("imgstat.csv",format="csv")

t=Time(table['mjd'],format="mjd",scale="utc")
pos=SkyCoord(ra=table['ra'],dec=table['de'],frame="icrs",unit="deg")

loc=EarthLocation(lat=52.8344*u.deg,lon=6.3785*u.deg,height=10*u.m)

pa=pos.transform_to(AltAz(obstime=t,location=loc))

mjd0=np.floor(np.min(table['mjd']))

plt.figure(figsize=(20,10))
plt.subplot(411)

plt.plot(table['mjd']-mjd0,table['mean'],label='Brightness')
plt.plot(table['mjd']-mjd0,table['std'],label='Variation')
plt.ylabel("ADU")
plt.legend()
plt.subplot(412)
plt.plot(table['mjd']-mjd0,pa.az.degree)
plt.ylabel("Azimuth (deg)")
plt.subplot(413)
plt.plot(table['mjd']-mjd0,pa.alt.degree)
plt.ylabel("Altitude (deg)")
plt.subplot(414)
plt.plot(table['mjd']-mjd0,table['rmsx'],label='RA')
plt.plot(table['mjd']-mjd0,table['rmsy'],label='Dec')
plt.ylim(0,60)
plt.ylabel("Residual (arcseconds)")
plt.xlabel("MJD - %.0f"%mjd0)
plt.legend()
plt.savefig("imgstat.png")
