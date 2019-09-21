#!/usr/bin/env python3
from __future__ import print_function
import configparser
import argparse
from spacetrack import SpaceTrackClient
from shutil import copyfile
import datetime
from io import BytesIO
from zipfile import ZipFile
from urllib.request import urlopen
import re
import os

if __name__ == '__main__':
    # Read commandline options
    conf_parser = argparse.ArgumentParser(description='Update TLEs from' +
                                                      ' online sources')
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify configuration file. If no file" +
                             " is specified 'configuration.ini' is used.",
                             metavar="FILE")

    args = conf_parser.parse_args()

    # Process commandline options and parse configuration
    cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    if args.conf_file:
        cfg.read([args.conf_file])
    else:
        cfg.read('configuration.ini')

    tle_path = cfg.get('Common', 'tle_path')

    now = datetime.datetime.utcnow()
    time = now.strftime("%Y%m%d_%H%M%S")

    # Get Space Track TLEs
    catalog_tle = os.path.join(tle_path, 'catalog.tle')
    st = SpaceTrackClient(identity=cfg.get('Credentials', 'st-username'),
                          password=cfg.get('Credentials', 'st-password'))

    data = st.tle_latest(iter_lines=True, epoch='>now-30',
                         ordinal=1, format='3le')

    with open(catalog_tle, 'w') as fp:
        for line in data:
            # Fix missing leading zeros
            line = re.sub("^1     ", "1 0000", line)
            line = re.sub("^2     ", "2 0000", line)
            line = re.sub("^1    ", "1 000", line)
            line = re.sub("^2    ", "2 000", line)
            line = re.sub("^1   ", "1 00", line)
            line = re.sub("^2   ", "2 00", line)
            line = re.sub("^1  ", "1 0", line)
            line = re.sub("^2  ", "2 0", line)
            fp.write(line + '\n')

    copyfile(catalog_tle, os.path.join(tle_path, time + '_catalog.txt'))

    # Get classified TLEs
    resp = urlopen("http://www.prismnet.com/~mmccants/tles/classfd.zip")
    zipfile = ZipFile(BytesIO(resp.read()))
    zipfile.extractall(path=tle_path)
    classfd_tle = os.path.join(tle_path, 'classfd.tle')

    content = ''
    outsize = 0
    with open(classfd_tle, 'rb') as infile:
        content = infile.read()
    with open(classfd_tle, 'wb') as output:
        for line in content.splitlines():
            outsize += len(line) + 1
            output.write(line + b'\n')

    copyfile(classfd_tle, os.path.join(tle_path, time + '_classfd.txt'))

    # Get int TLEs
    resp = urlopen("http://www.prismnet.com/~mmccants/tles/inttles.zip")
    zipfile = ZipFile(BytesIO(resp.read()))
    zipfile.extractall(path=tle_path)
    int_tle = os.path.join(tle_path, 'inttles.tle')

    content = ''
    outsize = 0
    with open(int_tle, 'rb') as infile:
        content = infile.read()
    with open(int_tle, 'wb') as output:
        for line in content.splitlines():
            outsize += len(line) + 1
            output.write(line + b'\n')

    copyfile(int_tle, os.path.join(tle_path, time + '_inttles.txt'))

    # Create bulk catalog
    catalogs = [catalog_tle, classfd_tle]
    with open(os.path.join(tle_path, 'bulk.tle'), 'w') as outfile:
        for fname in catalogs:
            with open(fname) as infile:
                outfile.write(infile.read())
