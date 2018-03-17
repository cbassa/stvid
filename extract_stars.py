#!/usr/bin/env python3
from __future__ import print_function
import os
import subprocess
import glob

# Get files
files=sorted(glob.glob("2*.fits"))

# Get sextractor location
env=os.getenv("ST_DATADIR")

# Loop over files
for file in files:
    # Format command
    command="sextractor %s -c %s/sextractor/default.sex"%(file,env)

    # Run sextractor
    output=subprocess.check_output(command,shell=True,stderr=subprocess.STDOUT)

    # Rename file
    if os.path.exists("test.cat"):
        os.rename("test.cat",file+".cat")
    else:
        print("test.cat does not exist")
