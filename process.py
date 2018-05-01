#!/usr/bin/env python
from __future__ import print_function
import glob
import numpy as np
from stio import fourframe
from stars import pixel_catalog, generate_star_catalog
from astrometry import calibrate_from_reference
import astropy.units as u
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import EarthLocation, AltAz
import warnings

if __name__ == "__main__":
    # Set warnings
    warnings.filterwarnings("ignore", category=UserWarning, append=True)
    warnings.simplefilter("ignore", AstropyWarning)

    # Set location
    loc=EarthLocation(lat=52.8344*u.deg,lon=6.3785*u.deg,height=10*u.m)

    # Get files
    files=sorted(glob.glob("2*.fits"))

    # Statistics file
    fstat=open("imgstat.csv","w")
    fstat.write("fname,mjd,ra,de,rmsx,rmsy,mean,std,nstars,nused\n")
    
    # Loop over files
    for fname in files:
        # Generate star catalog
        pix_catalog=generate_star_catalog(fname)

        # Calibrate astrometry
        calibrate_from_reference(fname,"test.fits",pix_catalog)

        # Stars available and used
        nused=np.sum(pix_catalog.flag==1)
        nstars=pix_catalog.nstars

        # Get properties
        ff=fourframe(fname)

        print("%s,%.8lf,%.6f,%.6f,%.3f,%.3f,%.3f,%.3f,%d,%d"%(ff.fname,ff.mjd,ff.crval[0],ff.crval[1],3600*ff.crres[0],3600*ff.crres[1],np.mean(ff.zavg),np.std(ff.zavg),nstars,nused))
        fstat.write("%s,%.8lf,%.6f,%.6f,%.3f,%.3f,%.3f,%.3f,%d,%d\n"%(ff.fname,ff.mjd,ff.crval[0],ff.crval[1],3600*ff.crres[0],3600*ff.crres[1],np.mean(ff.zavg),np.std(ff.zavg),nstars,nused))


    fstat.close()
