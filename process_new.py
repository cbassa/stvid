#!/usr/bin/env python3
import os
import glob

import configparser

from stvid.fourframe import FourFrame
from stvid.fourframe import Observation

if __name__ == "__main__":
    # Read configuration file
    config_file = "config_new.ini"
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ":"))
    result = cfg.read([config_file])

    # Observer settings
    site_id = cfg.getint("Observer", "cospar")

    # Extract colors for TLE files
    colors, abbrevs, tlefiles, catalognames = [], [], [], []
    for key, value in cfg.items("Elements"):
        if "tlefile" in key:
            tlefiles.append(value)
        elif "color" in key:
            colors.append(value)
        elif "name" in key:
            catalognames.append(value)
        elif "abbrev" in key:
            abbrevs.append(value)
    color_detected = cfg.get("LineDetection", "color")

    # Colormap
    cmap = cfg.get("DiagnosticPlot", "colormap")

    fname = "/data3/satobs/test/185300/processed/2022-03-24T18:53:20.708.fits"
    fnames = sorted(glob.glob("/data3/satobs/test/185300/processed/2*.fits"))
    fnames = sorted(glob.glob("/data/projects/stvid/test/2*.fits"))
    fnames = sorted(glob.glob("/data3/satobs/test/asi174mm/2*.fits"))
    #    fname = "/data3/satobs/test/2022-04-02T21:35:17.038.fits"

    for fname in fnames:
        print(fname)
        ff = FourFrame(fname)

        # Output file root
        froot = os.path.splitext(fname)[0]

        # Generate predictions
        predictions = ff.generate_satellite_predictions(cfg)

        # Detect 3D lines
        tracks = ff.find_tracks_by_hough3d(cfg)

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
            obs.append(Observation(ff, t.tmid, t.x0, t.y0, site_id,
                                   t.satno, t.cospar, t.catalogname))

        # Write observations
        for o in obs:
            iod_line = o.to_iod_line()

            # Open file
            outfname = f"{froot}_{o.satno:05d}_{o.catalogname}.dat"
            with open(outfname, "w") as fp:
                fp.write(f"{iod_line}\n")
            print(iod_line, o.catalogname)

        # Generate plots
        ff.diagnostic_plot(predictions, None, None, cfg)
        for track, o in zip(tracks, obs):
            ff.diagnostic_plot(predictions, track, o, cfg)
