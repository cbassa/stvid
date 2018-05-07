#!/usr/bin/env python
from __future__ import print_function
import glob
import numpy as np
from stvid.stio import fourframe
from stvid.stars import generate_star_catalog
from stvid.astrometry import calibrate_from_reference
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import EarthLocation
import warnings
import configparser
import argparse
import os

if __name__ == "__main__":

    # Read commandline options
    conf_parser = argparse.ArgumentParser(description='Process captured' +
                                                      ' video frames.')
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify configuration file. If no file" +
                             " is specified 'configuration.ini' is used.",
                             metavar="FILE")
    conf_parser.add_argument("-d", "--directory",
                             help="Specify directory of observations. If no" +
                             " directory is specified parent will be used.",
                             metavar='DIR', dest='file_dir', default=".")

    args = conf_parser.parse_args()

    # Process commandline options and parse configuration
    cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    if args.conf_file:
        cfg.read([args.conf_file])
    else:
        cfg.read('configuration.ini')

    # Set warnings
    warnings.filterwarnings("ignore", category=UserWarning, append=True)
    warnings.simplefilter("ignore", AstropyWarning)

    # Set location
    loc = EarthLocation(lat=cfg.getfloat('Common', 'observer_lat')*u.deg,
                        lon=cfg.getfloat('Common', 'observer_lon')*u.deg,
                        height=cfg.getfloat('Common', 'observer_el')*u.m)

    # Get files
    files = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits")))

    # Statistics file
    fstat = open(os.path.join(args.file_dir, "imgstat.csv"), "w")
    fstat.write("fname,mjd,ra,de,rmsx,rmsy,mean,std,nstars,nused\n")

    # Loop over files
    for fname in files:
        # Generate star catalog
        pix_catalog = generate_star_catalog(fname)

        # Calibrate astrometry
        calibrate_from_reference(fname,
                                 os.path.join(args.file_dir, "test.fits"),
                                 pix_catalog)

        # Stars available and used
        nused = np.sum(pix_catalog.flag == 1)
        nstars = pix_catalog.nstars

        # Get properties
        ff = fourframe(fname)

        print(("%s,%.8lf,%.6f,%.6f,%.3f,%.3f," +
              "%.3f,%.3f,%d,%d") % (ff.fname, ff.mjd, ff.crval[0],
                                    ff.crval[1], 3600*ff.crres[0],
                                    3600*ff.crres[1], np.mean(ff.zavg),
                                    np.std(ff.zavg), nstars, nused))
        fstat.write(("%s,%.8lf,%.6f,%.6f,%.3f,%.3f,%.3f," +
                     "%.3f,%d,%d\n") % (ff.fname, ff.mjd, ff.crval[0],
                                        ff.crval[1], 3600*ff.crres[0],
                                        3600*ff.crres[1], np.mean(ff.zavg),
                                        np.std(ff.zavg), nstars, nused))

    fstat.close()
