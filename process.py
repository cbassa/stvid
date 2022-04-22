#!/usr/bin/env python3
from __future__ import print_function
import glob
import numpy as np
from stvid.stio import FourFrame
from stvid.stars import generate_star_catalog
from stvid.stars import store_calibration
from stvid.stars import pixel_catalog
from stvid.astrometry import calibrate_from_reference
from stvid.astrometry import is_calibrated
from stvid.astrometry import generate_reference_with_anet
from stvid.satellite import generate_satellite_predictions
from stvid.satellite import find_hough3d_lines
from stvid.extract import extract_tracks
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import EarthLocation
from astropy.time import Time # For getting IERS table in single-thread session
import multiprocessing as mp
import warnings
import configparser
import argparse
import os
from termcolor import colored
import time
import shutil
import sys

""" 
process.py - Utility to process stvid/acquire.py FITS images to detect and 
   extract satellite positions and create IODs.

Terminal output Color Codes:
    GREEN:   Calibrated file
    RED:     Not calibrated file
    YELLOW:  Computing astrometric calibration
    BLUE:    Classified satellite
    GREY:    Catalog satellite
    MAGENTA: UNID satellite
"""

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def process_loop(fname):
    """
    Thread to process satobs FourFrame FITS files in a multi-thread compatible manner
    """

    # Generate star catalog
    if not os.path.exists(fname + ".cat"):
        pix_catalog = generate_star_catalog(fname)
    else:
        pix_catalog = pixel_catalog(fname + ".cat")

    # Calibrate from reference
    calibrate_from_reference(fname, "test.fits", pix_catalog)

    # Store calibration
    store_calibration(pix_catalog, fname + ".cal")

    # Generate satellite predictions
    generate_satellite_predictions(fname)

    # Detect lines with 3D Hough transform
    ids = find_hough3d_lines(fname, nhoughmin, houghrmin)

    # Get properties
    ff = FourFrame(fname)

    # Extract tracks
    if is_calibrated(ff):
        screenoutput_idents = extract_tracks(fname, trkrmin, drdtmin, drdtmax, trksig, ntrkmin, root_dir, results_dir, tle_dir)
    else:
        screenoutput_idents = None 

    # Stars available and used
    nused = np.sum(pix_catalog.flag == 1)
    nstars = pix_catalog.nstars

    # Write output
    screenoutput = "%s %10.6f %10.6f %4d/%4d %5.1f %5.1f %6.2f +- %6.2f" % (
        ff.fname, ff.crval[0], ff.crval[1], nused, nstars,
        3600.0 * ff.crres[0], 3600.0 * ff.crres[1], np.mean(
            ff.zavg), np.std(ff.zavg))

    if is_calibrated(ff):
        screenoutput = colored(screenoutput, "green")
    else:
        screenoutput = colored(screenoutput, "red")

    imgstat_output = ("%s,%.8lf,%.6f,%.6f,%.3f,%.3f,%.3f," + "%.3f,%d,%d\n") % (
        (ff.fname, ff.mjd, ff.crval[0], ff.crval[1],
            3600 * ff.crres[0], 3600 * ff.crres[1], np.mean(
                ff.zavg), np.std(ff.zavg), nstars, nused))

    # Move processed files
    shutil.move(fname, os.path.join(processed_dir, fname))
    shutil.move(fname + ".png", os.path.join(processed_dir, fname + ".png"))
    shutil.move(fname + ".id", os.path.join(processed_dir, fname + ".id"))
    shutil.move(fname + ".cat", os.path.join(processed_dir, fname + ".cat"))
    shutil.move(fname + ".cal", os.path.join(processed_dir, fname + ".cal"))

    return (screenoutput, imgstat_output, screenoutput_idents)

if __name__ == "__main__":
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description='Process captured' +
                                          ' video frames.')
    conf_parser.add_argument("-c",
                             "--conf_file",
                             help="Specify configuration file. If no file" +
                             " is specified 'configuration.ini' is used.",
                             metavar="FILE")
    conf_parser.add_argument("-d",
                             "--directory",
                             help="Specify directory of observations. If no" +
                             " directory is specified parent will be used.",
                             metavar='DIR',
                             dest='file_dir',
                             default=".")

    args = conf_parser.parse_args()

    # Process commandline options and parse configuration
    cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    conf_file = args.conf_file if args.conf_file else "configuration.ini"
    result = cfg.read([conf_file])

    if not result:
        print("Could not read config file: %s\nExiting..." % conf_file)
        sys.exit()

    # Set warnings
    warnings.filterwarnings("ignore", category=UserWarning, append=True)
    warnings.simplefilter("ignore", AstropyWarning)

    # Set location
    loc = EarthLocation(lat=cfg.getfloat('Common', 'observer_lat') * u.deg,
                        lon=cfg.getfloat('Common', 'observer_lon') * u.deg,
                        height=cfg.getfloat('Common', 'observer_height') * u.m)

    # Extract settings
    drdtmin = cfg.getfloat('Processing', 'drdtmin')
    drdtmax = cfg.getfloat('Processing', 'drdtmax')
    trksig = cfg.getfloat('Processing', 'trksig')
    trkrmin = cfg.getfloat('Processing', 'trkrmin')
    ntrkmin = cfg.getint('Processing', 'ntrkmin')
    houghrmin = cfg.getfloat('Processing', 'houghrmin')
    nhoughmin = cfg.getint('Processing', 'nhoughmin')
    nstarsmin = cfg.getint('Processing', 'nstarsmin')
    
    # Move to processing directory
    os.chdir(args.file_dir)

    # Force single-threaded IERS table update, if needed
    t = Time.now()
    _ = t.ut1

    # Statistics file
    fstat = open("imgstat.csv", "w")
    fstat.write("fname,mjd,ra,de,rmsx,rmsy,mean,std,nstars,nused\n")

    # Directory logic
    file_dir = os.path.abspath(args.file_dir.rstrip("/"))
    root_dir = os.path.split(file_dir)[0]
    tle_dir = cfg.get('Common', 'tle_path')
    if cfg.has_option('Common', 'results_path'):
        results_dir = os.path.join(cfg.get('Common', 'results_path'),
                                   os.path.split(root_dir)[-1])
    else:
        results_dir = root_dir
    processed_dir = os.path.join(file_dir, "processed")

    # Create output dirs
    if not os.path.exists(os.path.join(results_dir, "classfd")):
        os.makedirs(os.path.join(results_dir, "classfd"))
    if not os.path.exists(os.path.join(results_dir, "catalog")):
        os.makedirs(os.path.join(results_dir, "catalog"))
    if not os.path.exists(os.path.join(results_dir, "unid")):
        os.makedirs(os.path.join(results_dir, "unid"))
    if not os.path.exists(os.path.join(results_dir, "not_seen")):
        os.makedirs(os.path.join(results_dir, "not_seen"))
    if not os.path.exists(processed_dir):
        os.makedirs(processed_dir)

    cpu_count = os.cpu_count()
    if (cpu_count > 1):
        print("Processing with {} threads".format(cpu_count))

    # Forever loop
    while True:
        # Get files
        fnames = sorted(glob.glob("2*.fits"))

        # Create reference calibration file in single threaded environment
        if not os.path.exists("test.fits"):
            solved = False
            # Loop over files to find a suitable calibration file
            for fname in fnames:
                # Generate star catalog
                pix_catalog = generate_star_catalog(fname)

                # Solve
                if pix_catalog.nstars > nstarsmin:
                    print(colored("Computing astrometric calibration for %s" % fname, "yellow"))
                    solved = generate_reference_with_anet(fname, "")

                # Break when solved
                if solved:
                    break

        p = mp.Pool(processes=cpu_count)

        try:
            # Loop over files in parallel, print output every batch size of cpu_count
            satobs_chunks = chunks(fnames,cpu_count)
            for chk in satobs_chunks:
                for result in p.map(process_loop, chk):
                    (screenoutput, fileoutput, screenoutput_idents) = result    

                    fstat.write(fileoutput)
                    print(screenoutput)

                    if screenoutput_idents is not None:
                        for [outfilename, iod_line, color] in screenoutput_idents:
                            print(colored(iod_line,color))
                            # Write iodline
                            with open(outfilename, "a") as fp:
                                fp.write("%s\n" % iod_line)

            p.close()
            p.join()
        except KeyboardInterrupt:
            fstat.close()
            p.close()
            sys.exit()

        # Sleep
        try:
            print("File queue empty, waiting for new files...\r", end = '')
            time.sleep(1)
        except KeyboardInterrupt:
            fstat.close()
            sys.exit()

    # Close files
    fstat.close()
