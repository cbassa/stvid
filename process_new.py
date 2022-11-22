#!/usr/bin/env python3
import os
import sys
import glob
import time
import yaml

import argparse
import configparser

import numpy as np

import warnings

from termcolor import colored

from stvid import calibration
from stvid.fourframe import FourFrame
from stvid.fourframe import Observation
from stvid.fourframe import AstrometricCatalog

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
    acat = AstrometricCatalog(cfg.getfloat("Astrometry", "maximum_magnitude"))
            
    # Start processing loop
    while True:
       # Get files
        fnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits")))

        # Create reference calibration file
        calfname = os.path.join(args.file_dir, "test.fits")
        if not os.path.exists(calfname):
            solved = False

            # Loop over files to find a suitable calibration file
            for fname in fnames:
                # Was this file already tried?
                if not os.path.exists(os.path.join(args.file_dir, f"{fname}.cat")):
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

        # Loop over files
        for fname in fnames:
            # File root
            froot = os.path.splitext(fname)[0]
            
            # Find stars
            if not os.path.exists(f"{froot}_stars.cat"):
                scat = calibration.generate_star_catalog(fname)
            else:
                scat = calibration.read_star_catalog(fname)

            # Calibrate
            if not os.path.exists(f"{froot}_calib.wcs"):
                w, rmsx, rmsy, nused, is_calibrated = calibration.calibrate(fname, cfg, acat, scat, wref, tref)

                # Attempt plate solve
#                if not is_calibrated and scat.nstars > nstarsmin:
#                    print(colored(f"Computing astrometric calibration for {fname}", "yellow"))
#                    wtmp, ttmp = calibration.plate_solve(fname, cfg, calfname)
                    
#                    # Retry calibration
#                    if wtmp is not None:
#                        wref, tref = wtmp, ttmp
#                        w, rmsx, rmsy, nused, is_calibrated = calibration.calibrate(fname, cfg, acat, scat, wref, tref)

                # Log output
                output = f"{os.path.basename(fname)},{w.wcs.crval[0]:.6f},{w.wcs.crval[1]:.6f},{rmsx:.3f},{rmsy:.3f},{nused}/{scat.nstars}"
                if is_calibrated:
                    color = "green"
                else:
                    color = "red"
                print(colored(output, color))

            # Skip if png exists
            if os.path.exists(f"{froot}_0.png"):
                continue
                
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
                           "observer": ff.observer}
            
            # Loop over tracks
            ident_dicts = []
            obs = []
            satno = 90000
            for t in tracks:
                # Identify
                ident, is_identified = t.identify(predictions, satno, "22 500A", None, cfg, abbrevs, tlefiles)
                if not is_identified:
                    satno += 1

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
                ident_dict["measurement"] = single_measurement
                ident_dict["iod_line"] = iod_line
                ident_dict["measurements"] = multiple_measurements
                ident_dict["iod_lines"] = [line for line in iod_lines]
                ident_dicts.append(ident_dict)

                # Store observation
                obs.append(Observation(ident.satno, ident.catalogname, iod_line, iod_lines))

                
            # Store output
            if ident_dicts is not []:
                output_dict["satellites"] = ident_dicts

                with open(f"{ff.froot}_data.yaml", "w") as fp:
                    yaml.dump(output_dict, fp, sort_keys=False)
                
            # Write observations
            for o in obs:
                # Open file
                outfname = f"{ff.froot}_{o.satno:05d}_{o.catalogname}.dat"
                with open(outfname, "w") as fp:
                    fp.write(f"{o.iod_line}\n")
                outfname = f"{ff.froot}_{o.satno:05d}_{o.catalogname}_m.dat"
                with open(outfname, "w") as fp:
                    for iod_line in o.iod_lines:
                        fp.write(f"{iod_line}\n")
                if o.catalogname == "catalog":
                    color = "grey"
                elif o.catalogname == "classfd":
                    color = "blue"
                elif o.catalogname == "unid":
                    color = "magenta"
                print(colored(o.iod_line, color))

            # Generate plots
            ff.diagnostic_plot(predictions, None, None, cfg)
            for track, o in zip(tracks, obs):
                ff.diagnostic_plot(predictions, track, o, cfg)

#        # Get files
#        fitsfnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits")))
#        procfnames = sorted(glob.glob(os.path.join(args.file_dir, "2*.fits.png")))
#        fnames = [fname for fname in fitsfnames if f"{fname}.png" not in procfnames]

        # Sleep
        try:
            sys.exit()
            print("File queue empty, waiting for new files...\r", end = "")
            time.sleep(10)
        except KeyboardInterrupt:
            sys.exit()
