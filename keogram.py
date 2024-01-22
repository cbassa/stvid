#!/usr/bin/env python3
from __future__ import print_function
import glob
import numpy as np
from stvid.stio import fourframe
from stvid.config import add_argument_conf_file, load_config
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from astropy.time import Time


def generate_keogram(path):
    # Get files
    fnames = sorted(glob.glob(os.path.join(path, "processed/2*.fits")))

    # Allocate arrays
    nx = len(fnames)
    ny = cfg.getint('Camera', 'camera_y')
    ixmid = cfg.getint('Camera', 'camera_x') // 2
    keogram = np.zeros(nx * ny).reshape(ny, nx)
    mjds = np.zeros(nx)

    # Get data
    for i, fname in enumerate(fnames):
        if i % 10 == 0:
            print(i, fname)

        # Read file
        ff = fourframe(fname)

        # Extract data
        keogram[:, i] = ff.zavg[:, ixmid]
        mjds[i] = ff.mjd

    # Save npy arrays
    np.save(os.path.join(path, "mjds"), mjds)
    np.save(os.path.join(path, "keogram"), keogram)


if __name__ == "__main__":
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description='Process captured' +
                                          ' video frames.')
    conf_parser = add_argument_conf_file(conf_parser)
    conf_parser.add_argument("-d",
                             "--directory",
                             help="Specify directory of observations. If no" +
                             " directory is specified parent will be used.",
                             metavar='DIR',
                             dest='file_dir',
                             default=".")

    args = conf_parser.parse_args()
    cfg = load_config(args.conf_files)

    # Generate keogram is it does not exist
    if not os.path.exists(os.path.join(args.file_dir, "keogram.npy")):
        generate_keogram(args.file_dir)

    # Load data
    keogram = np.load(os.path.join(args.file_dir, "keogram.npy"))
    mjds = np.load(os.path.join(args.file_dir, "mjds.npy"))

    # Time limits
    tmin = mdates.date2num(Time(np.min(mjds), format="mjd").datetime)
    tmax = mdates.date2num(Time(np.max(mjds), format="mjd").datetime)

    # Plot keogram
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.imshow(np.log10(keogram),
              aspect="auto",
              origin="lower",
              cmap="magma",
              extent=[tmin, tmax, 0, 1])

    ax.axes.get_yaxis().set_visible(False)
    ax.xaxis_date()
    date_format = mdates.DateFormatter("%F\n%H:%M:%S")
    ax.xaxis.set_major_formatter(date_format)
    fig.autofmt_xdate(rotation=0, ha="center")

    plt.tight_layout()
    plt.savefig(os.path.join(args.file_dir, "keogram.png"))
