#!/usr/bin/env python3
from __future__ import print_function
import os
import configparser
import argparse
import subprocess

if __name__ == '__main__':
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description='Plate solve FITS file' +
                                                      ' and add WCS on header')
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

    path = args.file_dir
    extension = 'fits'
    files = sorted((f for f in os.listdir(path) if f.endswith('.' + extension)),
                   key=lambda x: os.path.getctime(os.path.join(path, x)))

    if files:
        file_for_astrometry = os.path.join(path, files[0])
        print("Found " + file_for_astrometry + " for astrometric solving.")

        sex_config = cfg.get('Astrometry', 'sex_config')
        low_app = cfg.get('Astrometry', 'low_app')
        high_app = cfg.get('Astrometry', 'high_app')

        # Format command
        command = "solve-field -O -y -u app -L %s -H %s --downsample 2 " % (low_app, high_app) + \
                  "--use-sextractor --sextractor-config %s --x-column X_IMAGE " % sex_config + \
                  "--y-column Y_IMAGE  --sort-column MAG_AUTO --sort-ascending " + \
                  "--no-plots -T -N %s/test.fits %s" % (path, file_for_astrometry)

        # Run sextractor
        subprocess.run(command, shell=True, stderr=subprocess.STDOUT)

    else:
        print("No fits file found for astrometric solving.")
