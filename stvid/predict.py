#!/usr/bin/env python3
import os
import numpy as np
import subprocess

from astropy.io import ascii
from astropy.coordinates import SkyCoord

class Prediction:
    """Prediction class"""

    def __init__(self, satno, mjd, ra, dec, x, y, state, tlefile, age):
        self.satno = satno
        self.age = age
        self.mjd = mjd
        self.t = 86400 * (self.mjd - self.mjd[0])
        self.texp = self.t[-1] - self.t[0]
        self.ra = ra
        self.dec = dec
        self.x = x
        self.y = y
        self.state = state
        self.tlefile = tlefile

def generate_satellite_predictions(ff, cfg):
    # Output file name
    outfname = f"{ff.fname}.csv"

    # Run predictions
    if not os.path.exists(outfname):
        # Extract parameters
        nfd = ff.nfd
        texp = ff.texp
        nmjd = int(np.ceil(texp))
        ra0, de0 = ff.crval[0], ff.crval[1]
        radius = np.sqrt(ff.wx * ff.wx + ff.wy * ff.wy)
        lat = cfg.getfloat("Observer", "latitude")
        lon = cfg.getfloat("Observer", "longitude")
        height = cfg.getfloat("Observer", "height")
    
        # Format command
        command = f"satpredict -t {nfd} -l {texp} -n {nmjd} -L {lon} -B {lat} -H {height} -o {outfname} -R {ra0} -D {de0} -r {radius}"
        for key, value in cfg.items("Elements"):
            if "tlefile" in key:
                command += f" -c {value}"

        # Run command
        output = subprocess.check_output(command,
                                         shell=True,
                                         stderr=subprocess.STDOUT)

    # Read results
    d = ascii.read(outfname, format="csv")

    # Compute frame coordinates
    p = SkyCoord(ra=d["ra"], dec=d["dec"], unit="deg", frame="icrs")
    x, y = p.to_pixel(ff.w, 0)
    
    # Loop over satnos
    satnos = np.unique(d["satno"])
    predictions = []
    for satno in satnos:
        c = d["satno"] == satno
        tlefile = np.unique(d["tlefile"][c])[0]
        age = np.unique(np.asarray(d["age"])[c])[0]
        p = Prediction(satno, np.asarray(d["mjd"])[c], np.asarray(d["ra"])[c], np.asarray(d["dec"])[c], x[c], y[c], np.array(d["state"])[c], tlefile, age)
        predictions.append(p)
        
    return predictions
