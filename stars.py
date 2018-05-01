#!/usr/bin/env python
from __future__ import print_function
import os
import subprocess
import numpy as np

class pixel_catalog:
    """Pixel catalog"""

    def __init__(self,fname):
        d=np.loadtxt(fname)
        if len(d.shape)==2:
            self.x=d[:,0]
            self.y=d[:,1]
            self.mag=d[:,2]
            self.ra=np.empty_like(self.x)
            self.dec=np.empty_like(self.x)
            self.imag=np.empty_like(self.x)
            self.flag=np.zeros_like(self.x)
            self.nstars=len(self.mag)
        else:
            self.x=None
            self.y=None
            self.mag=None
            self.ra=None
            self.dec=None
            self.imag=None
            self.flag=None
            self.nstars=0
            
def generate_star_catalog(fname):
    # Skip if file already exists
    if not os.path.exists(fname+".cat"):
        # Get sextractor location
        env=os.getenv("ST_DATADIR")

        # Format command
        command="sextractor %s -c %s/sextractor/default.sex"%(fname,env)

        # Run sextractor
        output=subprocess.check_output(command,shell=True,stderr=subprocess.STDOUT)

        # Rename file
        if os.path.exists("test.cat"):
            os.rename("test.cat",fname+".cat")

    return pixel_catalog(fname+".cat")
