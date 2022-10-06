#!/usr/bin/env python3
import os
import sys
import glob
import time

import argparse
import configparser

import numpy as np

import warnings

from termcolor import colored

from stvid import calibration

#from stvid.fourframe import FourFrame
#from stvid.fourframe import Observation
#from stvid.fourframe import AstrometricCatalog

#from stvid.stars import pixel_catalog
#from stvid.stars import store_calibration
#from stvid.stars import generate_star_catalog
#from stvid.astrometry import calibrate_from_reference
#from stvid.astrometry import generate_reference_with_anet

from astropy.utils.exceptions import AstropyWarning

if __name__ == "__main__":
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description="Process captured" +
                                          " video frames.")
    conf_parser.add_argument("-c",
                             "--conf_file",
                             help="Specify configuration file. If no file" +
                             " is specified 'configuration.ini' is used.",
                             metavar="FILE")
    conf_parser.add_argument("-d",
                             "--directory",
                             help="Specify directory of observations. If no" +
                             " directory is specified parent will be used.",
                             metavar="DIR",
                             dest="file_dir",
                             default=".")
    args = conf_parser.parse_args()
    
    # Read configuration file
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ":"))
    conf_file = args.conf_file if args.conf_file else "configuration.ini"
    result = cfg.read([conf_file])
    if not result:
        print("Could not read config file: %s\nExiting..." % conf_file)
        sys.exit()

    # Set warnings
    warnings.filterwarnings("ignore", category=UserWarning, append=True)
    warnings.simplefilter("ignore", AstropyWarning)

    
    # Observer settings
    nstarsmin = cfg.getint("Processing", "nstarsmin")

    
    # Extract abbrevs for TLE files
    abbrevs, tlefiles = [], []
    for key, value in cfg.items("Elements"):
        if "tlefile" in key:
            tlefiles.append(os.path.basename(value))
        elif "abbrev" in key:
            abbrevs.append(value)

    # Read astrometric catalog
    acat = calibration.AstrometricCatalog(cfg.getfloat("Astrometry", "maximum_magnitude"))
            
    # Start processing loop
    solved = False
    while True:
        # Get files
        fitsfnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits")))
        procfnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits.png")))
        fnames = [fname for fname in fitsfnames if f"{fname}.png" not in procfnames]

        # Loop over files
        for fname in fnames:
            # Find stars
            scat = calibration.generate_star_catalog(fname)

            # Plate solve
            if not solved and scat.nstars > nstarsmin:
                print(colored(f"Computing astrometric calibration for {fname}", "yellow"))
                wref, tref = calibration.plate_solve(fname, cfg)

                # Mark as solved
                if wref is not None:
                    solved = True
            
            # Calibrate
            w, rmsx, rmsy, nused = calibration.calibrate(fname, cfg, acat, scat, wref, tref)
            output = f"{os.path.basename(fname)},{w.wcs.crval[0]:.6f},{w.wcs.crval[1]:.6f},{rmsx:.3f},{rmsy:.3f},{nused}/{scat.nstars}"

            print(scat.flag)
            
            print(colored(output, "green"))

            # Stars available and used
            #nused = np.sum(pix_catalog.flag == 1)
            #nstars = pix_catalog.nstars
#            nstars = starcatalog.nstars
            
            # Write output
#            screenoutput = "%s %10.6f %10.6f %4d/%4d %5.1f %5.1f %6.2f +- %6.2f" % (
#                os.path.basename(ff.fname), ff.crval[0], ff.crval[1], nused, nstars,
#                3600.0 * ff.crres[0], 3600.0 * ff.crres[1], np.mean(
#                    ff.zavg), np.std(ff.zavg))
#
#            if ff.is_calibrated():
#                color = "green"
#            else:
#                color = "red"
#            print(colored(screenoutput, color))
            
        # Sleep
        try:
            sys.exit()
            print("File queue empty, waiting for new files...\r", end = "")
            time.sleep(10)
        except KeyboardInterrupt:
            sys.exit()
