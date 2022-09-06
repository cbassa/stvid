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

from stvid.fourframe import FourFrame
from stvid.fourframe import Observation

from stvid.stars import pixel_catalog
from stvid.stars import store_calibration
from stvid.stars import generate_star_catalog
from stvid.astrometry import calibrate_from_reference
from stvid.astrometry import generate_reference_with_anet

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
    conf_parser.add_argument("-r",
                             "--reprocess",
                             help="Remove processed files and start from scratch.",
                             action="store_true")
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
            tlefiles.append(value)
        elif "abbrev" in key:
            abbrevs.append(value)

    # Remove processed files
    if args.reprocess:
        for f_pattern in ["test.fits","*.png","*.cat","*.cal","*.csv","*.dat"]:
            for files in glob.glob(os.path.join(args.file_dir, f_pattern)):
                os.remove(files)

    # Start processing loop
    while True:
        # Get files
        fitsfnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits")))
        procfnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits.png")))
        fnames = [fname for fname in fitsfnames if f"{fname}.png" not in procfnames]

        # Loop over files
        for fname in fnames:
            # Create reference calibration file in single threaded environment
            calfname = os.path.join(args.file_dir, "test.fits")
            if not os.path.exists(calfname):
                solved = False
                # Loop over files to find a suitable calibration file
                for fname in fnames:
                    # Was this file already tried?
                    if not os.path.exists(f"{fname}.cat"):
                        # Generate star catalog
                        pix_catalog = generate_star_catalog(fname)

                        # Solve
                        if pix_catalog.nstars > nstarsmin:
                            print(colored(f"Computing astrometric calibration for {fname}", "yellow"))
                            solved = generate_reference_with_anet(fname, "", calfname)

                        # Break when solved
                        if solved:
                            break
            else:
                # test.fits exists, so calibration has been solved
                solved = True

            # Exit on failure to solve
            if not solved:
                break

            # Generate star catalog
            if not os.path.exists(f"{fname}.cat"):
                pix_catalog = generate_star_catalog(fname)
            else:
                pix_catalog = pixel_catalog(f"{fname}.cat")

            # Calibrate from reference
            calibrate_from_reference(fname, calfname, pix_catalog)

            # Store calibration
            store_calibration(pix_catalog, f"{fname}.cal")

            # Read Fourframe
            ff = FourFrame(fname)

            # Find stars
#            starcatalog = ff.generate_star_catalog()

            # Calibrate
#            wref = ff.find_calibration(cfg)
            
            # Stars available and used
            nused = np.sum(pix_catalog.flag == 1)
            nstars = pix_catalog.nstars
            
            # Write output
            screenoutput = "%s %10.6f %10.6f %4d/%4d %5.1f %5.1f %6.2f +- %6.2f" % (
                os.path.basename(ff.fname), ff.crval[0], ff.crval[1], nused, nstars,
                3600.0 * ff.crres[0], 3600.0 * ff.crres[1], np.mean(
                    ff.zavg), np.std(ff.zavg))

            if ff.is_calibrated():
                color = "green"
            else:
                color = "red"
            print(colored(screenoutput, color))
            
            # Generate predictions
            predictions = ff.generate_satellite_predictions(cfg)

            # Find tracks
            if ff.is_calibrated():
                tracks = ff.find_tracks_by_hough3d(cfg)
            else:
                tracks = []

            # Identify tracks
            satno = 90000
            for t in tracks:
                is_identified = t.identify(predictions, satno, "22 500A", None, cfg, abbrevs, tlefiles)
                if not is_identified:
                    satno += 1

            # Format observations
            obs = []
            for t in tracks:
                # Add to observation
                obs.append(Observation(ff, t.tmid, t.x0, t.y0, ff.site_id,
                                       t.satno, t.cospar, t.catalogname))

            # Write observations
            for o in obs:
                iod_line = o.to_iod_line()

                # Open file
                outfname = f"{ff.froot}_{o.satno:05d}_{o.catalogname}.dat"
                with open(outfname, "w") as fp:
                    fp.write(f"{iod_line}\n")
                if o.catalogname == "catalog":
                    color = "grey"
                elif o.catalogname == "classfd":
                    color = "blue"
                elif o.catalogname == "unid":
                    color = "magenta"
                print(colored(iod_line, color))

            # Generate plots
            ff.diagnostic_plot(predictions, None, None, cfg)
            for track, o in zip(tracks, obs):
                ff.diagnostic_plot(predictions, track, o, cfg)


        # Sleep
        try:
            print("File queue empty, waiting for new files...\r", end = "")
            time.sleep(10)
        except KeyboardInterrupt:
            sys.exit()
