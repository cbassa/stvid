#!/usr/bin/env python3
import re
import os
import datetime
import argparse

from io import BytesIO
from shutil import copyfile
from zipfile import ZipFile
from urllib.request import urlopen

from spacetrack import SpaceTrackClient
from stvid.config import add_argument_conf_file, load_config


if __name__ == '__main__':
    conf_parser = argparse.ArgumentParser(description="Update TLEs from" +
                                                      " online sources")
    conf_parser = add_argument_conf_file(conf_parser)
    args = conf_parser.parse_args()

    cfg = load_config(args.conf_files)

    # Create TLE locatio
    tle_path = cfg.get("Elements", "tlepath")
    if not os.path.exists(tle_path):
        os.makedirs(tle_path)
    
    now = datetime.datetime.utcnow()
    time = now.strftime("%Y%m%d_%H%M%S")

    print("Get Space Track TLEs")
    catalog_tle = os.path.join(tle_path, "catalog.tle")
    st = SpaceTrackClient(identity=cfg.get("Credentials", "st-username"),
                          password=cfg.get("Credentials", "st-password"))

    data = st.tle_latest(iter_lines=True, epoch=">now-30",
                         ordinal=1, format="3le")

    with open(catalog_tle, "w") as fp:
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
            fp.write(line + "\n")

    copyfile(catalog_tle, os.path.join(tle_path, time + "_catalog.txt"))

    print("Get classified TLEs")
    resp = urlopen("https://mmccants.org/tles/classfd.zip")
    zipfile = ZipFile(BytesIO(resp.read()))
    zipfile.extractall(path=tle_path)
    classfd_tle = os.path.join(tle_path, "classfd.tle")

    content = ""
    outsize = 0
    with open(classfd_tle, "rb") as infile:
        content = infile.read()
    with open(classfd_tle, "wb") as output:
        for line in content.splitlines():
            outsize += len(line) + 1
            output.write(line + b"\n")

    copyfile(classfd_tle, os.path.join(tle_path, time + "_classfd.txt"))

    print("Get integrated TLEs")
    resp = urlopen("https://mmccants.org/tles/inttles.zip")
    zipfile = ZipFile(BytesIO(resp.read()))
    zipfile.extractall(path=tle_path)
    int_tle = os.path.join(tle_path, "inttles.tle")

    content = ""
    outsize = 0
    with open(int_tle, "rb") as infile:
        content = infile.read()
    with open(int_tle, "wb") as output:
        for line in content.splitlines():
            outsize += len(line) + 1
            output.write(line + b"\n")

    copyfile(int_tle, os.path.join(tle_path, time + "_inttles.txt"))

    print("Get supplemental Starlink TLEs")
    resp = urlopen("https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle")
    starlink_tle = os.path.join(tle_path, "starlink.tle")
    lines = resp.read().decode("utf-8")
    with open(starlink_tle, "w") as fp:
        fp.write(lines)

    copyfile(starlink_tle, os.path.join(tle_path, time + "_starlink.txt"))
        
    print("Get supplemental OneWeb TLEs")
    resp = urlopen("https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=oneweb&FORMAT=tle")
    oneweb_tle = os.path.join(tle_path, "oneweb.tle")
    lines = resp.read().decode("utf-8")
    with open(oneweb_tle, "w") as fp:
        fp.write(lines)
 
    copyfile(oneweb_tle, os.path.join(tle_path, time + "_oneweb.txt"))
    
    print("Create bulk catalog")
    catalogs = [catalog_tle, classfd_tle]
    with open(os.path.join(tle_path, "bulk.tle"), "w") as outfile:
        for fname in catalogs:
            with open(fname) as infile:
                outfile.write(infile.read())
