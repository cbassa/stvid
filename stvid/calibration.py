#!/usr/bin/env python3
import os

import subprocess

import numpy as np

import astropy.units as u
from astropy import wcs
from astropy.time import Time
from astropy.io import fits
from astropy.coordinates import SkyCoord, ICRS, FK5

class AstrometricCatalog:
    """AstrometricCatalog class"""

    def __init__(self, maxmag):
        # Filename
        fname = os.path.normpath(os.path.join(
            os.path.dirname(__file__),
            "..",
            "data/tyc2.fits"))
        
        # Read catalog
        hdu = fits.open(fname)
        ra = hdu[1].data.field("RA")
        dec = hdu[1].data.field("DEC")
        mag = hdu[1].data.field("MAG_VT")
        hdu.close()
        
        # Select stars
        c = mag < maxmag

        self.ra = ra[c]
        self.dec = dec[c]
        self.mag = mag[c]
        self.nstars = len(self.ra)

class StarCatalog:
    """StarCatalog class"""

    def __init__(self, fname):
        # Load catalog
        d = np.loadtxt(fname)
        if len(d.shape) == 2:
            self.x = d[:, 0]
            self.y = d[:, 1]
            self.mag = d[:, 2]
            self.ra = np.empty_like(self.x)
            self.dec = np.empty_like(self.x)
            self.imag = np.empty_like(self.x)
            self.flag = np.zeros_like(self.x).astype("bool")
            self.nstars = len(self.mag)
        else:
            self.x = None
            self.y = None
            self.mag = None
            self.ra = None
            self.dec = None
            self.imag = None
            self.flag = None
            self.nstars = 0


def generate_star_catalog(fname):
    # Source-extractor configuration file
    path = os.path.normpath(os.path.join(os.path.dirname(__file__),
                                         "..",
                                         "source-extractor"))
    conffname = os.path.join(path, "default.sex")

    # Output catalog name
    froot = os.path.splitext(fname)[0]
    outfname = f"{froot}_stars.cat"

    # Skip if file already exists
    if not os.path.exists(outfname):
        # Format command
        command = f"sextractor {fname} -c {conffname} -CATALOG_NAME {outfname}"

        # Add sextractor config path to environment
        env = dict(os.environ)
        env["SEXTRACTOR_CFG"] = path
        
        # Run sextractor
        output = subprocess.check_output(command, shell=True, env=env,
                                         stderr=subprocess.STDOUT)
            
    return StarCatalog(outfname)


def read_star_catalog(fname):
    # Output catalog name
    froot = os.path.splitext(fname)[0]
    outfname = f"{froot}_stars.cat"

    return StarCatalog(outfname)


def plate_solve(fname, cfg, store_as_fname=None):
    # Arguments for solve-field
    if cfg.has_option("Astrometry", "solve-field_args"):
        cmd_args = cfg.get("Astrometry", "solve-field_args")
    else:
        cmd_args = ""

    # File root
    froot = os.path.splitext(fname)[0]
        
    # Generate command
    command = f"solve-field {cmd_args} {fname}"

    # Run command
    try:
        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)

        # Read header
        hdu = fits.open(f"{froot}.new")
        hdr = hdu[0].header
        data = hdu[0].data
        hdu[0].header["NAXIS"] = 2
        w = wcs.WCS(hdu[0].header)
        t = Time(hdu[0].header["MJD-OBS"], format="mjd", scale="utc")
        hdu.close()
    except OSError:
        w, t = None, None

    # Solved
    if w is not None:
        solved = True
    else:
        solved = False
        
    # Remove temporary files
    extensions = [".new", ".axy", "-objs.png", ".rdls", ".solved", "-indx.xyls",
                  ".match", ".corr", ".wcs", "-indx.png", "-ngc.png"]
    for extension in extensions:
        try:
            os.remove(f"{froot}{extension}")
        except OSError:
            pass

    # Store calibrated file
    if solved and store_as_fname is not None:
        hdu = fits.PrimaryHDU(header=hdr, data=data)
        hdu.writeto(store_as_fname, overwrite=True, output_verify="ignore")
        
    return w, t

def read_calibration(fname):
    # Read header
    hdu = fits.open(fname)
    hdu[0].header["NAXIS"] = 2
    w = wcs.WCS(hdu[0].header)
    t = Time(hdu[0].header["MJD-OBS"], format="mjd", scale="utc")
    hdu.close()

    return w, t

def calibrate(fname, cfg, astcat, pixcat, wref, tref):
    # Read FITS
    hdu = fits.open(fname)
    data = hdu[0].data
    header = hdu[0].header
    hdu.close()
    t = Time(header["MJD-OBS"], format="mjd", scale="utc")

    # Check for sidereal tracking
    try:
        tracked = bool(hdu[0].header['TRACKED'])
    except KeyError:
        tracked = False

    # Copy WCS
    w = wref.deepcopy()
        
    # Apply tracking
    if not tracked:
        # Reference position
        pref = SkyCoord(ra=wref.wcs.crval[0],
                        dec=wref.wcs.crval[1],
                        unit="deg",
                        frame="icrs").transform_to(FK5(equinox=tref))

        # RA shift
        dra = (t.sidereal_time("mean", "greenwich") - tref.sidereal_time("mean", "greenwich"))

        # Updated position
        p = FK5(ra=pref.ra + dra, dec=pref.dec, equinox=t).transform_to(ICRS)
        w.wcs.crval = np.array([p.ra.degree, p.dec.degree])

    # Exit on empty star catalog
    if pixcat.nstars == 0:
        return w, 0, 0, 0, False
        
    # Image size
    nx, ny = header["NAXIS1"], header["NAXIS2"]

    # Sky coordinates for astrometric standards
    p = SkyCoord(ra=astcat.ra, dec=astcat.dec, unit="deg", frame="fk5")
    
    # Pixel coordinates
    x, y = pixcat.x, pixcat.y

    # Settings
    rmin = 120 
    order = 1
    niter = 10
    rmsx, rmsy = 0, 0
    for i in range(niter):
        # Compute sky coordinates of detected stars
        pc = SkyCoord.from_pixel(x, y, w, 1)

        # Match
        idx, r, _ = pc.match_to_catalog_sky(p)
        c = r < rmin * u.arcsec
        nstars_used = np.sum(c)
        
        # Refit
        if nstars_used > 4:
            w = fit_wcs(x[c], y[c], p[idx[c]].ra.degree,
                        p[idx[c]].dec.degree, nx // 2, ny // 2, order)

            # Compute residuals
            rx, ry = residuals(x[c], y[c], p[idx[c]].ra.degree, p[idx[c]].dec.degree, w)

            r = np.sqrt(rx**2 + ry**2)        
            rms = np.sqrt(np.sum(r**2) / r.size)
            rmsx, rmsy = np.std(rx), np.std(ry)
            rmin = 2 * rms

    # Store flags
    pixcat.flag = c
            
    # CTYPE
    if order > 1:
        ctype1, ctype2 = "RA---TAN-SIP", "DEC--TAN-SIP"
    else:
        ctype1, ctype2 = "RA---TAN", "DEC--TAN"
             
    # Keywords to add
    whdr = {"CRPIX1": w.wcs.crpix[0], "CRPIX2": w.wcs.crpix[1],
            "CRVAL1": w.wcs.crval[0], "CRVAL2": w.wcs.crval[1],
            "CD1_1": w.wcs.cd[0, 0], "CD1_2": w.wcs.cd[0, 1],
            "CD2_1": w.wcs.cd[1, 0], "CD2_2": w.wcs.cd[1, 1],
            "CTYPE1": ctype1, "CTYPE2": ctype2, "CUNIT1": "DEG",
            "CUNIT2": "DEG", "CRRES1": rmsx / 3600, "CRRES2": rmsy / 3600}

    # Add SIP
    if order > 1:
        whdr["A_ORDER"] = order
        whdr["B_ORDER"] = order
        for i in range(order):
            for j in range(order):
                key = f"A_{i}_{j}"
                value = w.sip.a[i, j]
                whdr[key] = value
                key = f"B_{i}_{j}"
                value = w.sip.b[i, j]
                whdr[key] = value

    # Add keywords
    hdr = hdu[0].header
    for k, v in whdr.items():
        hdr[k] = v

    # Write FITS
    hdu = fits.PrimaryHDU(header=hdr, data=data)
    hdu.writeto(fname, overwrite=True, output_verify="ignore")

    # Is calibrated?
    sx = 3600 * np.sqrt(w.wcs.cd[0, 0] ** 2 + w.wcs.cd[1, 0] ** 2)
    sy = 3600 * np.sqrt(w.wcs.cd[0, 1] ** 2 + w.wcs.cd[1, 1] ** 2)
    if (rmsx < 1e-3) | (rmsy < 1e-3) | (rmsx / sx > 2) | (rmsy / sy > 2):
        is_calibrated = False
    else:
        is_calibrated = True

    # Log calibration
    froot = os.path.splitext(fname)[0]
    outfname = f"{froot}_calib.wcs"
    with open(outfname, "w") as fp:
        pass
        
    return w, rmsx, rmsy, nstars_used, is_calibrated

def solve_linear_equation(a, b):
    q, r = np.linalg.qr(a)
    y = np.dot(q.T, b)
    x = np.linalg.solve(r, y)
    return x

def fit_wcs(x, y, ra, dec, x0, y0, order):
    #ra0, dec0 = np.mean(ra), np.mean(dec) # Change with mean angle
    ra0, dec0 = ra[0], dec[0]
    dx, dy = x - x0, y - y0

    ixs, iys = np.meshgrid(np.arange(order+1), np.arange(order+1))
    c = ixs + iys <= order
    ix, iy = ixs[c], iys[c]
    a = np.array([dx**ix[i] * dy**iy[i] for i in range(len(iy))]).T

    for k in range(5):
        w = wcs.WCS(naxis=2)
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cd = [[1.0, 0.0], [0.0, 1.0]]
        w.wcs.crval = [ra0, dec0]
        w.wcs.crpix = [0.0, 0.0]

        rx, ry = w.wcs_world2pix(np.stack((ra, dec), axis=-1), 0).T

        ax = solve_linear_equation(a, rx)
        ay = solve_linear_equation(a, ry)

        ra0, dec0 = w.wcs_pix2world(([[ax[0], ay[0]]]), 0)[0]

    cd = np.array([[ax[1], ax[order+1]], [ay[1], ay[order+1]]])
    cdinv = np.linalg.inv(cd)

    axm = np.zeros_like(ixs).astype("float32")
    aym = np.zeros_like(ixs).astype("float32")
    for i in range(len(ix)):
        if ix[i] + iy[i]>=2:
            p = np.matmul(cdinv, np.array([ax[i], ay[i]]))
            axm[iy[i], ix[i]] = p[0]
            aym[iy[i], ix[i]] = p[1]    

    w = wcs.WCS(naxis=2)
    w.wcs.cd = cd
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [x0, y0]
    if order > 1:
        w.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
        w.sip = wcs.Sip(axm.T, aym.T, None, None, w.wcs.crpix)
    else:
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        
    return w

def residuals(xcen, ycen, ra, dec, w):
    pas = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")

    pstars = SkyCoord.from_pixel(xcen, ycen, w, 1)

    r = pas.separation(pstars)
    pa = pas.position_angle(pstars)
    rx, ry = r * np.sin(pa), r * np.cos(pa)

    return rx.to(u.arcsec).value, ry.to(u.arcsec).value
