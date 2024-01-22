#!/usr/bin/env python3
from __future__ import print_function
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
import argparse
import os

from stvid.config import add_argument_conf_file, load_config


if __name__ == "__main__":
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description='Plot image statistics')
    conf_parser = add_argument_conf_file(conf_parser)
    conf_parser.add_argument("-i",
                             "--input",
                             help="Specify file to be processed. If no file" +
                             " is specified ./imgstat.csv will be used.",
                             metavar='FILE',
                             default="./imgstat.csv")
    conf_parser.add_argument("-d",
                             "--directory",
                             help="Specify directory of observations. If no" +
                             " directory is specified parent will be used.",
                             metavar='DIR',
                             dest='file_dir',
                             default=".")
    conf_parser.add_argument(
        "-o",
        "--output",
        help="Specify output file. Default is 'imgstat.png'",
        metavar='FILE',
        default="./imgstat.png")

    args = conf_parser.parse_args()
    cfg = load_config(args.conf_files)

    # Move to processing directory
    os.chdir(args.file_dir)

    table = ascii.read(args.input, format="csv")

    t = Time(table['mjd'], format="mjd", scale="utc")
    tplot = mdates.date2num(t.datetime)
    pos = SkyCoord(ra=table['ra'], dec=table['de'], frame="icrs", unit="deg")

    # Set location
    loc = EarthLocation(lat=cfg.getfloat('Common', 'observer_lat') * u.deg,
                        lon=cfg.getfloat('Common', 'observer_lon') * u.deg,
                        height=cfg.getfloat('Common', 'observer_height') * u.m)

    pa = pos.transform_to(AltAz(obstime=t, location=loc))

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,
                                             1,
                                             figsize=(20, 10),
                                             sharex=True)

    date_format = mdates.DateFormatter("%F\n%H:%M:%S")
    fig.autofmt_xdate(rotation=0, ha="center")

    ax1.plot(tplot, table['mean'], label='Brightness')
    ax1.plot(tplot, table['std'], label='Variation')
    ax1.set_ylabel("Brightness (ADU)")

    ax1.xaxis_date()
    ax1.xaxis.set_major_formatter(date_format)
    ax1.legend()

    ax2.plot(tplot, pa.az.degree)
    ax2.xaxis_date()
    ax2.xaxis.set_major_formatter(date_format)
    ax2.set_ylabel("Azimuth (deg)")

    ax3.plot(tplot, pa.alt.degree)
    ax3.xaxis_date()
    ax3.xaxis.set_major_formatter(date_format)
    ax3.set_ylabel("Altitude (deg)")

    ax4.plot(tplot, table['rmsx'], label="RA")
    ax4.plot(tplot, table['rmsy'], label="Dec")
    ax4.xaxis_date()
    ax4.xaxis.set_major_formatter(date_format)
    ax4.set_xlabel("Time (UTC)")
    ax4.set_ylabel("Residuals (arcseconds)")
    ax4.legend()

    plt.tight_layout()
    plt.savefig(args.output)
