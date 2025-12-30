#!/usr/bin/env python3
import os
import sys
import glob
import time
import json
import datetime

import argparse
import configparser

import numpy as np

import warnings

from termcolor import colored

import multiprocessing as mp

from stvid import calibration
from stvid.fourframe import FourFrame
from stvid.fourframe import Observation
from stvid.fourframe import AstrometricCatalog

from astropy.utils.exceptions import AstropyWarning

def number_to_letter(n):
    # 
    if n == 0:
        return ""
    x = (n - 1) % 24
    letters = "ABCDEFGHJKLMNPQRSTUVWXYZ"
    rest = (n - 1) // 24
    if rest == 0:
        return letters[x]
    return number_to_letter(rest) + letters[x]

def cospar(nfd, number):
    t = datetime.datetime.strptime(nfd[:19], "%Y-%m-%dT%H:%M:%S")
    year = int(t.strftime("%y"))
    doy = int(t.strftime("%j")) + 500
    letter = number_to_letter(number)
    return f"{year:02d} {doy:03d}{letter:3s}"

def chunk_list(l, n):
    o = []
    for i in range(0, len(l), n):
        o.append(l[i:i + n])
    return o

def process_loop(fname):
    """
    Thread to process satobs FourFrame FITS files in a multi-thread compatible manner
    """

    # File root
    froot = os.path.splitext(fname)[0]
    
    # Find stars
    if not os.path.exists(f"{froot}_stars.cat"):
        scat = calibration.generate_star_catalog(fname)
    else:
        scat = calibration.read_star_catalog(fname)

    # Calibrate
    screenoutput = None
    if not os.path.exists(f"{froot}_calib.wcs"):
        w, rmsx, rmsy, nused, is_calibrated = calibration.calibrate(fname, cfg, acat, scat, wref, tref)

        # Attempt plate solve
#        if not is_calibrated and scat.nstars > nstarsmin:
#            print(colored(f"Computing astrometric calibration for {fname}", "yellow"))
#            wtmp, ttmp = calibration.plate_solve(fname, cfg, calfname)
            
#            # Retry calibration
#            if wtmp is not None:
#                wref, tref = wtmp, ttmp
#                w, rmsx, rmsy, nused, is_calibrated = calibration.calibrate(fname, cfg, acat, scat, wref, tref)

        # Log output
        output = f"{os.path.basename(fname)} {w.wcs.crval[0]:10.6f} {w.wcs.crval[1]:10.6f} {rmsx:6.2f} {rmsy:6.2f} {nused}/{scat.nstars}"
        if is_calibrated:
            color = "green"
        else:
            color = "red"
        screenoutput = colored(output, color)

    # Skip if png exists
    if os.path.exists(f"{froot}_0.png"):
        return
        
    # Read Fourframe
    ff = FourFrame(fname, cfg)
        
    # Generate predictions
    predictions = ff.generate_satellite_predictions(cfg)

    # Find tracks
    if ff.is_calibrated():
        tracks = ff.find_tracks_by_hough3d(cfg)
    else:
        tracks = []

    # Output dictionary
    output_dict = {"site_id": ff.site_id,
                   "latitude": ff.lat,
                   "longitude": ff.lon,
                   "height": ff.height,
                   "observer": ff.observer,
                   "start": ff.nfd,
                   "exptime": ff.texp,
                   "ra": ff.ra0,
                   "dec": ff.dec0,
                   "sx": ff.sx,
                   "sy": ff.sy,
                   "wx": ff.wx,
                   "wy": ff.wy}
    
    # Loop over tracks
    ident_dicts = []
    obs = []
    satno = 90000
    number = 1
    for t in tracks:
        # Identify
        ident, is_identified = t.identify(predictions, satno, cospar(ff.nfd, number), None, cfg, abbrevs, tlefiles)
        if not is_identified:
            satno += 1
            number += 1

        # Identification dictionary
        ident_dict = {"satno": ident.satno,
                      "cospar": ident.cospar,
                      "tlefile": ident.tlefile,
                      "catalogname": ident.catalogname}

        # Save track
        t.save(f"{ff.froot}_{ident.satno:05d}_{ident.catalogname}.csv", ff)

        # Measure single position
        m = t.measure_single_position(ff)
        iod_line = m.to_iod_line(ff, ident)

        # Measure multiple position
        ms = t.measure_multiple_positions(ff)
        iod_lines = [mt.to_iod_line(ff, ident) for mt in ms]

        # Add to dictionary
        single_measurement = {"time": m.t.isot,
                              "ra": float(m.ra),
                              "dec": float(m.dec),
                              "drxdt": float(m.drxdt),
                              "drydt": float(m.drydt)}
        multiple_measurements = [{"time": mt.t.isot,
                                  "ra": float(mt.ra),
                                  "dec": float(mt.dec),
                                  "drxdt": float(mt.drxdt),
                                  "drydt": float(mt.drydt)} for mt in ms]
#        ident_dict["measurement"] = single_measurement
#        ident_dict["iod_line"] = iod_line
        ident_dict["measurements"] = multiple_measurements
        ident_dict["iod_lines"] = [line for line in iod_lines]
        ident_dicts.append(ident_dict)

        # Store observation
        obs.append(Observation(ident.satno, ident.catalogname, iod_line, iod_lines))

        
    # Store output
    if ident_dicts is not []:
        output_dict["satellites"] = ident_dicts

        with open(f"{ff.froot}_data.json", "w") as fp:
            json.dump(output_dict, fp)
        
    # Write observations
    screenoutput_idents = []
    for o in obs:
        # Open file
        outfname = f"{ff.froot}_{o.satno:05d}_{o.catalogname}.dat"
        with open(outfname, "w") as fp:
            fp.write(f"{o.iod_line}\n")
        outfname = f"{ff.froot}_{o.satno:05d}_{o.catalogname}_m.dat"
        with open(outfname, "w") as fp:
            for iod_line in o.iod_lines:
                fp.write(f"{iod_line}\n")
        if o.catalogname == "classfd":
            color = "blue"
        elif o.catalogname == "unid":
            color = "magenta"
        else:
            color = "grey"
        screenoutput_idents.append(colored(o.iod_line, color))

    # Generate plots
    ff.diagnostic_plot(predictions, None, None, cfg)
    for track, o in zip(tracks, obs):
        ff.diagnostic_plot(predictions, track, o, cfg)

    # Clean up
    del ff
    for t in tracks:
        del t
    for o in obs:
        del o
        
    return (screenoutput, screenoutput_idents)

if __name__ == "__main__":
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description="Process captured" +
                                          " video frames.")
    conf_parser.add_argument("-c",
                             "--conf_file",
                             help="Specify configuration file. If no file" +
                             " is specified 'configuration.ini' is used.",
                             action="append",
                             nargs="?",
                             metavar="FILE")
    conf_parser.add_argument("-d",
                             "--directory",
                             help="Specify directory of observations. If no" +
                             " directory is specified parent will be used.",
                             metavar="DIR",
                             dest="file_dir",
                             default=".")
    conf_parser.add_argument("-b",
                             "--batch",
                             help="Batch process observations, exit when done.",
                             action="store_true")
    conf_parser.add_argument("-r",
                             "--reprocess",
                             help="Remove processed files and start from scratch.",
                             action="store_true")
    conf_parser.add_argument("-C",
                             "--cpu_count",
                             help="Number of threads to use (overrides value from configuration).",
                             type=int, default=None)
    conf_parser.add_argument("-w",
                             "--wait",
                             help="Delay before processing new files (seconds, default: 10).",
                             type=int, default=10)
    args = conf_parser.parse_args()
    
    # Read configuration file
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ":"))
    conf_file = args.conf_file if args.conf_file else "configuration.ini"
    result = cfg.read(conf_file)
    if not result:
        print("Could not read config file: %s\nExiting..." % conf_file)
        sys.exit()

    # Set warnings
    warnings.filterwarnings("ignore", category=UserWarning, append=True)
    warnings.simplefilter("ignore", AstropyWarning)
    
    # Observer settings
    nstarsmin = cfg.getint("Astrometry", "min_stars")

    # Extract abbrevs for TLE files
    abbrevs, tlefiles = [], []
    for key, value in cfg.items("Elements"):
        if "tlefile" in key:
            tlefiles.append(os.path.basename(value))
        elif "abbrev" in key:
            abbrevs.append(value)

    # Remove processed files
    if args.reprocess:
        for f_pattern in ["test.fits", "2*_*.png", "2*_stars.cat", "2*_calib.wcs",
                          "2*_*.csv", "2*_*.dat", "2*_data.json"]:
            for files in glob.glob(os.path.join(args.file_dir, f_pattern)):
                os.remove(files)

    # Read astrometric catalog
    acat = AstrometricCatalog(cfg.getfloat("Astrometry", "max_magnitude"))

    # Start calibration loop
    while True:
        # Get files without star catalogs
        fitsfnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits")))
        froots = [os.path.splitext(fitsname)[0] for fitsname in fitsfnames]
        fnames = [f"{froot}.fits" for froot in froots if not os.path.exists(f"{froot}_stars.cat")]

        # Create reference calibration file
        calfname = os.path.join(args.file_dir, "test.fits")
        if not os.path.exists(calfname):
            solved = False
            wref = None

            # Loop over files to find a suitable calibration file
            for fname in fnames:
                # Generate star catalog
                scat = calibration.generate_star_catalog(fname)

                # Solve
                if scat.nstars > nstarsmin:
                    print(colored(f"Computing astrometric calibration for {fname}", "yellow"))
                    wref, tref = calibration.plate_solve(fname, cfg, calfname)
                    if wref is not None:
                        solved = True

                # Break when solved
                if solved:
                    break
        else:
            # test.fits exists, so calibration has been solved
            solved = True

            # Read calibration
            wref, tref = calibration.read_calibration(calfname)

        # Break when solved
        if solved:
            print("Calibration succeeded!")
            break
            
        try:
            if(args.batch):
                sys.exit()
            print("File queue empty, waiting for new files...\r", end = "")
            time.sleep(args.wait)
        except KeyboardInterrupt:
            sys.exit()

    # Get number of CPUs for multiprocessing
    if not args.cpu_count:
        if cfg.has_option("LineDetection", "cpu_count"):
            cpu_count = cfg.getint("LineDetection", "cpu_count")
        else:
            cpu_count = mp.cpu_count()
    else:
        cpu_count = args.cpu_count
    print(f"Processing with {cpu_count} threads")
        
    # Processing loop
    while True:
        # Get unprocessed files
        fitsfnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits")))
        froots = [os.path.splitext(fitsname)[0] for fitsname in fitsfnames]
        fnames = [f"{froot}.fits" for froot in froots if not os.path.exists(f"{froot}_0.png")][:100]

        # Process files
        p = mp.Pool(processes=cpu_count)

        try:
            chunks = chunk_list(fnames, cpu_count)
            for chunk in chunks:
                for result in p.map(process_loop, chunk):
                    (screenoutput, screenoutput_idents) = result

                    if screenoutput is not None:
                        print(screenoutput)
                    for screenoutput_ident in screenoutput_idents:
                        print(screenoutput_ident)

            p.close()
            p.join()
        except KeyboardInterrupt:
            p.close()
            p.join()

        # Sleep
        try:
            if(args.batch):
                sys.exit()
            print("File queue empty, waiting for new files...\r", end = "")
            time.sleep(args.wait)
        except KeyboardInterrupt:
            sys.exit()
