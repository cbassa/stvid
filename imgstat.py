#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
import configparser
import argparse
import os

# Read commandline options
conf_parser = argparse.ArgumentParser(description='Plot image statistics')
conf_parser.add_argument("-c", "--conf_file",
                         help="Specify configuration file. If no file" +
                         " is specified 'configuration.ini' is used.",
                         metavar="FILE")
conf_parser.add_argument("-i", "--input",
                         help="Specify file to be processed. If no file" +
                         " is specified ./imgstat.csv will be used.",
                         metavar='FILE', default="./imgstat.csv")
conf_parser.add_argument("-d", "--directory",
                         help="Specify directory of observations. If no" +
                         " directory is specified parent will be used.",
                         metavar='DIR', dest='file_dir', default=".")
conf_parser.add_argument("-o", "--output",
                         help="Specify output file. Default is 'imgstat.png'",
                         metavar='FILE', default="./imgstat.png")

args = conf_parser.parse_args()

# Process commandline options and parse configuration
cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
if args.conf_file:
    cfg.read([args.conf_file])
else:
    cfg.read('configuration.ini')

# Move to processing directory
os.chdir(args.file_dir)
    
table = ascii.read(args.input, format="csv")

t = Time(table['mjd'], format="mjd", scale="utc")
pos = SkyCoord(ra=table['ra'], dec=table['de'], frame="icrs", unit="deg")

# Set location
loc = EarthLocation(lat=cfg.getfloat('Common', 'observer_lat')*u.deg,
                    lon=cfg.getfloat('Common', 'observer_lon')*u.deg,
                    height=cfg.getfloat('Common', 'observer_el')*u.m)

pa = pos.transform_to(AltAz(obstime=t, location=loc))

mjd0 = np.floor(np.min(table['mjd']))

plt.figure(figsize=(20, 10))
plt.subplot(411)

plt.plot(table['mjd']-mjd0, table['mean'], label='Brightness')
plt.plot(table['mjd']-mjd0, table['std'], label='Variation')
plt.ylabel("ADU")
plt.legend()
plt.subplot(412)
plt.plot(table['mjd']-mjd0, pa.az.degree)
plt.ylabel("Azimuth (deg)")
plt.subplot(413)
plt.plot(table['mjd']-mjd0, pa.alt.degree)
plt.ylabel("Altitude (deg)")
plt.subplot(414)
plt.plot(table['mjd']-mjd0, table['rmsx'], label='RA')
plt.plot(table['mjd']-mjd0, table['rmsy'], label='Dec')
plt.ylim(0, 60)
plt.ylabel("Residual (arcseconds)")
plt.xlabel("MJD - %.0f" % mjd0)
plt.legend()
plt.savefig(args.output)
