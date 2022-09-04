#!/usr/bin/env python3
import os

import subprocess

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import astropy.units as u
from astropy import wcs
from astropy.time import Time
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord, FK5

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
            self.flag = np.zeros_like(self.x)
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


class Prediction:
    """Prediction class"""

    def __init__(self, satno, cospar, mjd, ra, dec, x, y, rx, ry, state, tlefile, age):
        self.satno = satno
        self.cospar = cospar
        self.age = age
        self.mjd = mjd
        self.t = 86400 * (self.mjd - self.mjd[0])
        self.texp = self.t[-1] - self.t[0]
        self.ra = ra
        self.dec = dec
        self.x = x
        self.y = y
        self.rx = rx
        self.ry = ry
        self.state = state
        self.tlefile = tlefile

    def position_and_velocity(self, t, dt):
        # Skip if predicted track is too short to fit a 3rd order polynomial
        if len(self.t) < 3:
            return np.nan, np.nan, np.nan, np.nan, np.nan

        # Create polynomials
        px = np.polyfit(self.t, self.rx, 2)
        py = np.polyfit(self.t, self.ry, 2)

        # Derivatives
        dpx = np.polyder(px)
        dpy = np.polyder(py)

        # Evaluate
        rx0, ry0 = np.polyval(px, t), np.polyval(py, t)
        drxdt, drydt = np.polyval(dpx, t), np.polyval(dpy, t)
        drdt = np.sqrt(drxdt**2 + drydt**2)
        pa = np.mod(np.arctan2(-drxdt, drydt), 2 * np.pi)
        dr = drdt * dt

        return rx0, ry0, drdt, pa, dr

    def predicted_track(self, tmin, tmid, tmax, ax, color):
        # Skip if predicted track is too short to fit a 3rd order polynomial
        if len(self.t) < 3:
            return np.nan, np.nan, np.nan, np.nan

        # Create polynomials
        px = np.polyfit(self.t, self.x, 2)
        py = np.polyfit(self.t, self.y, 2)

        # Derivatives
        dpx = np.polyder(px)
        dpy = np.polyder(py)

        # Evaluate
        x0, y0 = np.polyval(px, tmid), np.polyval(py, tmid)
        dxdt, dydt = np.polyval(dpx, tmid), np.polyval(dpy, tmid)

        # Extrema
        xmin = x0 + dxdt * (tmin - tmid)
        xmax = x0 + dxdt * (tmax - tmid)
        ymin = y0 + dydt * (tmin - tmid)
        ymax = y0 + dydt * (tmax - tmid)

        ax.plot([xmin, xmax], [ymin, ymax], color=color, linestyle="-")
        ax.plot(x0, y0, color=color, marker="o", markerfacecolor="none")

    def residual(self, t0, rx0, ry0):
        # Skip if predicted track is too short to fit a 3rd order polynomial
        if len(self.t) < 3:
            return np.nan, np.nan

        # Create polynomials
        px = np.polyfit(self.t, self.rx, 2)
        py = np.polyfit(self.t, self.ry, 2)

        # Derivatives
        dpx = np.polyder(px)
        dpy = np.polyder(py)

        # Evaluate
        rx, ry = np.polyval(px, t0), np.polyval(py, t0)
        drxdt, drydt = np.polyval(dpx, t0), np.polyval(dpy, t0)
        drdt = np.sqrt(drxdt**2 + drydt**2)
        pa = np.arctan2(-drxdt, drydt)

        # Compute cross-track, in-track residual
        drx, dry = rx0 - rx, ry0 - ry
        ca, sa = np.cos(pa), np.sin(pa)
        rm = ca * drx - sa * dry
        wm = sa * drx + ca * dry
        dtm = rm / drdt

        return dtm, wm


class Track:
    """Track class"""

    def __init__(self, t, x, y, z, ra, dec, rx, ry):
        self.x = x
        self.y = y
        self.t = t
        self.z = z
        self.ra = ra
        self.dec = dec
        self.rx = rx
        self.ry = ry
        self.n = len(x)
        self.satno = None
        self.cospar = None
        self.tlefile = None

        # Compute mean position
        self.tmin, self.tmax = np.min(self.t), np.max(self.t)
        self.tmid = 0.5 * (self.tmax + self.tmin)

        # Position and velocity on the sky
        px = np.polyfit(self.t - self.tmid, self.rx, 1)
        py = np.polyfit(self.t - self.tmid, self.ry, 1)
        self.rx0 = px[-1]
        self.ry0 = py[-1]
        self.drxdt = px[-2]
        self.drydt = py[-2]
        self.rxmin = self.rx0 + self.drxdt * (self.tmin - self.tmid)
        self.rxmax = self.rx0 + self.drxdt * (self.tmax - self.tmid)
        self.rymin = self.ry0 + self.drydt * (self.tmin - self.tmid)
        self.rymax = self.ry0 + self.drydt * (self.tmax - self.tmid)
        self.drdt = np.sqrt(self.drxdt**2 + self.drydt**2)
        self.pa = np.mod(np.arctan2(-self.drxdt, self.drydt), 2 * np.pi)
        self.dr = self.drdt * (self.tmax - self.tmin)

        # Position and velocity on the image
        self.px = np.polyfit(self.t - self.tmid, self.x, 1)
        self.py = np.polyfit(self.t - self.tmid, self.y, 1)
        self.x0 = self.px[-1]
        self.y0 = self.py[-1]
        self.dxdt = self.px[-2]
        self.dydt = self.py[-2]
        self.xp = np.polyval(self.px, self.t - self.tmid)
        self.yp = np.polyval(self.py, self.t - self.tmid)
        self.xmin = self.x0 + self.dxdt * (self.tmin - self.tmid)
        self.xmax = self.x0 + self.dxdt * (self.tmax - self.tmid)
        self.ymin = self.y0 + self.dydt * (self.tmin - self.tmid)
        self.ymax = self.y0 + self.dydt * (self.tmax - self.tmid)
        self.r = np.sqrt((self.xmax - self.xmin) ** 2 + (self.ymax - self.ymin) ** 2)

    def identify(self, predictions, satno, cospar, tlefile, cfg, abbrevs, tlefiles):
        # Identification settings
        rm_max = cfg.getfloat("Identification", "max_off_track_offset_deg")
        dtm_max = cfg.getfloat("Identification", "max_along_track_offset_s")
        dpa_max = cfg.getfloat("Identification", "max_direction_difference_deg")
        fdr_max = cfg.getfloat("Identification", "max_velocity_difference_percent")

        # Loop over predictions
        is_identified = False
        for p in predictions:
            # Compute identification constraints
            rx0, ry0, drdt, pa, dr = p.position_and_velocity(
                self.tmid, self.tmax - self.tmin
            )
            dtm, rm = p.residual(self.tmid, self.rx0, self.ry0)
            dpa = angle_difference(self.pa, pa) * 180 / np.pi
            fdr = (dr / self.dr - 1) * 100
            if (
                (np.abs(dtm) < dtm_max)
                & (np.abs(rm) < rm_max)
                & (np.abs(dpa) < dpa_max)
                & (np.abs(fdr) < fdr_max)
            ):
                satno = p.satno
                cospar = p.cospar
                tlefile = p.tlefile
                is_identified = True

        self.satno = satno
        self.cospar = cospar
        self.catalogname = "unid"

        # Get catalog abbreviation
        for abbrev, tfile in zip(abbrevs, tlefiles):
            if tfile == tlefile:
                self.catalogname = abbrev

        return is_identified

    def match_to_prediction(self, p, dt, w):
        # Return if predicted track is too short to fit a 3rd order polynomial
        if len(p.t) < 3:
            return False

        # Create polynomials
        px = np.polyfit(p.t, p.x, 2)
        py = np.polyfit(p.t, p.y, 2)

        # Compute extrema
        xmin = np.polyval(px, self.tmin)
        xmax = np.polyval(px, self.tmax)
        ymin = np.polyval(py, self.tmin)
        ymax = np.polyval(py, self.tmax)

        # Check if observed track endpoints match with prediction
        pmin = inside_selection_area(
            self.tmin,
            self.tmax,
            self.x0,
            self.y0,
            self.dxdt,
            self.dydt,
            xmin,
            ymin,
            dt,
            w,
        )
        pmax = inside_selection_area(
            self.tmin,
            self.tmax,
            self.x0,
            self.y0,
            self.dxdt,
            self.dydt,
            xmax,
            ymax,
            dt,
            w,
        )

        return pmin & pmax


class Observation:
    """Satellite observation"""

    def __init__(self, ff, t, x, y, site_id, satno, cospar, catalogname):
        self.t = t
        self.mjd = ff.mjd + self.t / 86400
        self.nfd = Time(self.mjd, format="mjd", scale="utc").isot
        self.x = x
        self.y = y
        self.site_id = site_id
        self.satno = satno
        self.cospar = cospar
        self.catalogname = catalogname

        # Compute coordinates
        p = SkyCoord.from_pixel(self.x, self.y, ff.w, 1)

        # Correct for motion of stationary camera
        if ff.tracked == False:
            t = Time(self.mjd, format="mjd")
            tmid = Time(ff.mjd, format="mjd") + 0.5 * ff.texp * u.s
            p = correct_stationary_coordinates(tmid, t, p, direction=-1)
        
        self.ra = p.ra.degree
        self.dec = p.dec.degree

    def to_iod_line(self):
        pstr = format_position(self.ra, self.dec)
        tstr = (
            self.nfd.replace("-", "").replace("T", "").replace(":", "").replace(".", "")
        )
        iod_line = "%05d %-9s %04d G %s 17 25 %s 37 S" % (
            self.satno,
            self.cospar,
            self.site_id,
            tstr,
            pstr,
        )
        return iod_line


class FourFrame:
    """Four Frame class"""

    def __init__(self, fname=None):
        if fname is None:
            # Initialize empty fourframe
            self.nx = 0
            self.ny = 0
            self.nz = 0
            self.mjd = -1
            self.nfd = None
            self.zavg = None
            self.zstd = None
            self.zmax = None
            self.znum = None
            self.dt = None
            self.site_id = 0
            self.observer = None
            self.texp = 0.0
            self.fname = None
            self.froot = None
            self.crpix = np.array([0.0, 0.0])
            self.crval = np.array([0.0, 0.0])
            self.cd = np.array([[1.0, 0.0], [0.0, 1.0]])
            self.ctype = ["RA---TAN", "DEC--TAN"]
            self.cunit = np.array(["deg", "deg"])
            self.crres = np.array([0.0, 0.0])
            self.ra0 = None
            self.dec0 = None
            self.tracked = False
            self.lat = None
            self.lon = None
            self.height = None
            self.has_location = False
        else:
            # Read FITS file
            hdu = fits.open(fname)

            # Read header
            header = hdu[0].header

            # Read image planes
            self.zavg, self.zstd, self.zmax, self.znum = hdu[0].data

            # Generate sigma frame
            self.zsig = (self.zmax - self.zavg) / (self.zstd + 1e-9)

            # Frame properties
            self.ny, self.nx = self.zavg.shape
            self.nz = header["NFRAMES"]

            # Read frame time oselfsets
            self.dt = np.array([header["DT%04d" % i] for i in range(self.nz)])

            # Read header
            self.mjd = header["MJD-OBS"]
            self.nfd = header["DATE-OBS"]
            self.site_id = header["COSPAR"]
            self.observer = header["OBSERVER"]
            self.texp = header["EXPTIME"]
            self.fname = fname
            self.froot = os.path.splitext(fname)[0]

            # Location info
            keys = ["SITELONG", "SITELAT", "ELEVATIO"]
            if np.all([key in header for key in keys]):
                self.lon = header["SITELONG"]
                self.lat = header["SITELAT"]
                self.height = header["ELEVATIO"]
                self.has_location = True
            else:
                self.lon = None
                self.lat = None
                self.height = None
                self.has_location = False
            
            # Astrometry keywords
            self.crpix = np.array([header["CRPIX1"], header["CRPIX2"]])
            self.crval = np.array([header["CRVAL1"], header["CRVAL2"]])
            self.cd = np.array(
                [
                    [header["CD1_1"], header["CD1_2"]],
                    [header["CD2_1"], header["CD2_2"]],
                ]
            )
            self.ctype = [header["CTYPE1"], header["CTYPE2"]]
            self.cunit = [header["CUNIT1"], header["CUNIT2"]]
            self.crres = np.array([header["CRRES1"], header["CRRES2"]])
            self.ra0 = self.crval[0]
            self.dec0 = self.crval[1]

            # Check for sidereal tracking
            try:
                self.tracked = bool(header["TRACKED"])
            except KeyError:
                self.tracked = False

            hdu.close()

        # Compute image properties
        self.sx = np.sqrt(self.cd[0, 0] ** 2 + self.cd[1, 0] ** 2)
        self.sy = np.sqrt(self.cd[0, 1] ** 2 + self.cd[1, 1] ** 2)
        self.wx = self.nx * self.sx
        self.wy = self.ny * self.sy
        self.zmaxmin = np.mean(self.zmax) - 2.0 * np.std(self.zmax)
        self.zmaxmax = np.mean(self.zmax) + 6.0 * np.std(self.zmax)
        self.zavgmin = np.mean(self.zavg) - 2.0 * np.std(self.zavg)
        self.zavgmax = np.mean(self.zavg) + 6.0 * np.std(self.zavg)
        self.zstdmin = np.mean(self.zstd) - 2.0 * np.std(self.zstd)
        self.zstdmax = np.mean(self.zstd) + 6.0 * np.std(self.zstd)
        self.znummin = 0
        self.znummax = self.nz
        self.zsigmin = 0
        self.zsigmax = 10

        # Setup WCS
        self.w = wcs.WCS(naxis=2)
        self.w.wcs.crpix = self.crpix
        self.w.wcs.crval = self.crval
        self.w.wcs.cd = self.cd
        self.w.wcs.ctype = self.ctype
        self.w.wcs.set_pv([(2, 1, 45.0)])

    def generate_satellite_predictions(self, cfg):
        # Output file name
        outfname = f"{self.froot}_predict.csv"

        # Run predictions
        if not os.path.exists(outfname):
            # Extract parameters
            nfd = self.nfd
            texp = self.texp
            nmjd = int(np.ceil(texp))
            ra0, de0 = self.crval[0], self.crval[1]
            radius = np.sqrt(self.wx * self.wx + self.wy * self.wy)

            # Use FITS location if available, config otherwise
            if self.has_location:
                lat = self.lat
                lon = self.lon
                height = self.height
            else:
                lat = cfg.getfloat("Observer", "latitude")
                lon = cfg.getfloat("Observer", "longitude")
                height = cfg.getfloat("Observer", "height")

            # Format command
            command = f"satpredict -t {nfd} -l {texp} -n {nmjd} -L {lon} -B {lat} -H {height}"
            command = command + f" -o {outfname} -R {ra0} -D {de0} -r {radius}"
            for key, value in cfg.items("Elements"):
                if "tlefile" in key:
                    command += f" -c {value}"
            # Run command
            output = subprocess.check_output(
                command, shell=True, stderr=subprocess.STDOUT
            )

        # Read results
        d = ascii.read(outfname, format="csv")

        # Coordinates (correct for image drift)
        p = SkyCoord(ra=d["ra"], dec=d["dec"], unit="deg", frame=FK5(equinox="J2000"))
        if self.tracked == False:
            t = Time(d["mjd"], format="mjd")
            tmid = Time(self.mjd, format="mjd") + 0.5 * self.texp * u.s
            p = correct_stationary_coordinates(tmid, t, p, direction=1)

        # Compute pixel positions
        x, y = p.to_pixel(self.w, 0)

        # Compute angular offsets
        rx, ry = deproject(self.ra0, self.dec0, p.ra.degree, p.dec.degree)

        # Loop over satnos
        satnos = np.unique(d["satno"])
        tlefiles = np.unique(d["tlefile"])
        predictions = []
        for tlefile in tlefiles:
            for satno in satnos:
                c = (d["satno"] == satno) & (d["tlefile"] == tlefile)
                if np.sum(c) == 0:
                    continue
                age = np.unique(np.asarray(d["age"])[c])[0]
                cospar = str(np.unique(np.asarray(d["cospar"])[c])[0])

                # Fix COSPAR designation for IOD format
                if len(cospar) > 1:
                    if cospar[2] != " ":
                        cospar = f"{cospar[0:2]} {cospar[2:]}"
                p = Prediction(
                    satno,
                    cospar,
                    np.asarray(d["mjd"])[c],
                    np.asarray(d["ra"])[c],
                    np.asarray(d["dec"])[c],
                    x[c],
                    y[c],
                    rx[c],
                    ry[c],
                    np.array(d["state"])[c],
                    tlefile,
                    age,
                )
                predictions.append(p)

        return predictions

    def in_frame(self, x, y):
        if (x >= 0) & (x <= self.nx) & (y >= 0) & (y <= self.ny):
            return True
        else:
            return False

    def find_tracks_by_hough3d(self, cfg):
        # Config settings
        sigma = cfg.getfloat("LineDetection", "trksig")
        trkrmin = cfg.getfloat("LineDetection", "trkrmin")
        ntrkmin = cfg.getfloat("LineDetection", "ntrkmin")

        # Find significant pixels
        c = self.zsig > sigma
        xm, ym = np.meshgrid(np.arange(self.nx), np.arange(self.ny))
        x, y = np.ravel(xm[c]), np.ravel(ym[c])
        znum = np.ravel(self.znum[c]).astype("int")
        zmax = np.ravel(self.zmax[c])
        sig = np.ravel(self.zsig[c])
        t = np.array([self.dt[i] for i in znum])

        # Compute ra/dec
        p = SkyCoord.from_pixel(x, y, self.w, 0)
        ra, dec = p.ra.degree, p.dec.degree

        # Compute angular offsets
        rx, ry = deproject(self.ra0, self.dec0, ra, dec)

        # Skip if not enough points
        if len(t) < ntrkmin:
            return []

        # Save points to file
        fname = f"{self.froot}_threshold.csv"
        with open(fname, "w") as fp:
            for i in range(len(t)):
                fp.write(f"{x[i]:f},{y[i]:f},{znum[i]:f}\n")

        # Run 3D Hough line-finding algorithm
        command = f"hough3dlines -dx {trkrmin} -minvotes {ntrkmin} -raw {fname}"

        try:
            output = subprocess.check_output(
                command, shell=True, stderr=subprocess.STDOUT
            )
        except Exception:
            return []

        # Decode output
        fname = f"{self.froot}_hough.csv"
        with open(fname, "w") as fp:
            fp.write("ax,ay,az,bx,by,bz,n\n")
            tracks = []
            for line in output.decode("utf-8").splitlines()[2:]:
                # lines.append(ThreeDLine(line, self.nx, self.ny, self.nz))
                ax, ay, az, bx, by, bz, n = decode_line(line)

                # Write result
                fp.write(f"{ax},{ay},{az},{bx},{by},{bz},{n}\n")

                # Select points
                f = (znum - az) / bz
                xr = ax + f * bx
                yr = ay + f * by
                r = np.sqrt((x - xr) ** 2 + (y - yr) ** 2)
                c = r < trkrmin

                # Number of selected points and unique times
                nsel = np.sum(c)
                nt = len(np.unique(t[c]))

                if (nsel > 0) & (nt > 1):
                    tracks.append(
                        Track(t[c], x[c], y[c], zmax[c], ra[c], dec[c], rx[c], ry[c])
                    )

        return tracks

    def generate_star_catalog(self):
        # Source-extractor configuration file
        conffname = os.path.normpath(os.path.join(os.path.dirname(__file__),
                                                  "..",
                                                  "source-extractor/default.sex"))

        # Output catalog name
        outfname = f"{self.froot}.cat"

        # Skip if file already exists
        if not os.path.exists(outfname):
            # Format command
            command = f"sextractor {self.fname} -c {conffname} -CATALOG_NAME {outfname}"

            # Run sextractor
            output = subprocess.check_output(command, shell=True,
                                             stderr=subprocess.STDOUT)
            
        return StarCatalog(outfname)

    def find_calibration(self, cfg):
        # Arguments for solve-field
        if cfg.has_option("Astrometry", "solve-field_args"):
            cmd_args = cfg.get("Astrometry", "solve-field_args")
        else:
            cmd_args = ""

        # Generate command
        command = f"solve-field {cmd_args} {self.fname}"

        # Run command
        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)

            # Did the astrometry succeed?
            solved = os.path.exists(f"{self.froot}.solved")

            # Read header
            hdu = fits.open(f"{self.froot}.wcs")
            hdu[0].header["NAXIS"] = 2
            w = wcs.WCS(hdu[0].header)
            hdu.close()
        except OSError:
            solved = False
            w = None

        # Remove temporary files
        extensions = [".new", ".axy", "-objs.png", ".wcs", ".rdls", ".solved", "-indx.xyls", ".match", ".corr"]
        for extension in extensions:
            try:
                os.remove(f"{self.froot}{extension}")
            except OSError:
                pass

        return w
        
    
    def is_calibrated(self):
        if (3600.0 * self.crres[0] < 1e-3) | \
           (3600.0 * self.crres[1] < 1e-3) | \
           (self.crres[0] / self.sx > 2.0) | \
           (self.crres[1] / self.sy > 2.0):
            return False
        else:
            return True

    
    def diagnostic_plot(self, predictions, track, obs, cfg):
        # Get info
        if track is not None:
            iod_line = obs.to_iod_line()
            satno = obs.satno
            catalogname = obs.catalogname
            outfname = f"{self.froot}_{satno:05d}_{catalogname}.png"
        else:
            iod_line = ""
            outfname = f"{self.fname}.png"

        # Configuration parameters
        color_detected = cfg.get("LineDetection", "color")
        cmap = cfg.get("DiagnosticPlot", "colormap")
        colors, abbrevs, tlefiles, catalognames = [], [], [], []
        for key, value in cfg.items("Elements"):
            if "tlefile" in key:
                tlefiles.append(os.path.basename(value))
            elif "color" in key:
                colors.append(value)
            elif "name" in key:
                catalognames.append(value)
            elif "abbrev" in key:
                abbrevs.append(value)

        # Create plot
        fig, ax = plt.subplots(figsize=(12, 10), dpi=75)

        # Set color indicating accurate astrometry
        if self.is_calibrated():
            color = "k"
        else:
            color = "r"
        
        ax.set_title(
            f'UT Date: {self.nfd} COSPAR ID: {self.site_id}\nR.A.: {self.crval[0]:10.6f} ({3600 * self.crres[0]:.1f}") Decl.: {self.crval[1]:10.6f} ({3600 * self.crres[1]:.1f}")\nFOV: {self.wx:.2f}$^\circ$x{self.wy:.2f}$^\circ$ Scale: {3600 * self.sx:.2f}"x{3600 * self.sy:.2f}" pix$^{{-1}}$\n\n{iod_line}',
            fontdict={"fontsize": 14, "horizontalalignment": "left"},
            loc="left", color=color,
        )

        ax.imshow(
            self.zmax,
            origin="lower",
            interpolation="none",
            vmin=self.zmaxmin,
            vmax=self.zmaxmax,
            cmap=cmap,
        )

        for p in predictions:
            self.plot_prediction(p, ax, tlefiles, colors, dt=0)

        if track is not None:
            ax.plot(track.xp, track.yp, color=color_detected, linestyle="-")
            ax.plot(
                track.x0,
                track.y0,
                color=color_detected,
                marker="o",
                markerfacecolor="none",
            )
            ax.text(
                track.x0,
                track.y0,
                f" {track.satno:05d}",
                color=color_detected,
                ha="center",
                in_layout=False,
            )

        ax.set_xlim(0, self.nx)
        ax.set_ylim(0, self.ny)
        ax.set_xlabel("x (pixel)")
        ax.set_ylabel("y (pixel)")
        #        ax.xaxis.set_ticklabels([])
        #        ax.yaxis.set_ticklabels([])

        # Create legend handles
        handles = []
        for catalogname, color in zip(catalognames, colors):
            handles.append(
                mlines.Line2D([], [], color=color, marker="", label=catalogname)
            )
        handles.append(
            mlines.Line2D(
                [],
                [],
                color=color_detected,
                marker="o",
                markerfacecolor="none",
                label="Detected",
            )
        )
        for state, linestyle in zip(
            ["Sunlit", "Penumbra", "Eclipsed"], ["solid", "dashed", "dotted"]
        ):
            handles.append(
                mlines.Line2D(
                    [], [], color="k", linestyle=linestyle, marker="", label=state
                )
            )
        ax.legend(
            handles=handles,
            ncol=7,
            bbox_to_anchor=(0.5, -0.1),
            loc="center",
            frameon=True,
        )

        plt.tight_layout()
        plt.savefig(outfname, bbox_inches="tight", pad_inches=0.5)
        plt.close()

    def plot_prediction(self, p, ax, tlefiles, colors, dt=2.0, w=10.0):
        color = "k"
        for tlefile, temp_color in zip(tlefiles, colors):
            if tlefile in p.tlefile:
                color = temp_color
        marker = "."

        # Get direction of motion
        dx, dy = p.x[-1] - p.x[0], p.y[-1] - p.y[0]
        theta = np.arctan2(dy, dx)
        sa, ca = np.sin(theta), np.cos(theta)

        # Start of trail
        dx = np.zeros(2)
        dy = np.array([w, -w])
        xs = ca * dx - sa * dy + p.x[0]
        ys = sa * dx + ca * dy + p.y[0]

        # Trail
        dx, dy = p.x - p.x[0], p.y - p.y[0]
        r = np.sqrt(dx**2 + dy**2)
        dx, dy = r, -w * np.ones_like(r)
        xt = ca * dx - sa * dy + p.x[0]
        yt = sa * dx + ca * dy + p.y[0]

        ax.plot(xs, ys, color=color)
        ax.plot(xs, ys, color=color, marker=marker)
        if theta < 0:
            ha = "left"
        else:
            ha = "right"

        if self.in_frame(xs[0], ys[0]):
            ax.text(
                xs[0], ys[0], f" {p.satno:05d} ", color=color, ha=ha, in_layout=False
            )

        for state, linestyle in zip(
            ["sunlit", "umbra", "eclipsed"], ["solid", "dashed", "dotted"]
        ):
            c = correct_bool_state(p.state == state)
            ax.plot(
                np.ma.masked_where(~c, xt),
                np.ma.masked_where(~c, yt),
                color=color,
                linestyle=linestyle,
            )

        return


def decode_line(line):
    p = line.split(" ")
    ax = float(p[0])
    ay = float(p[1])
    az = float(p[2])
    bx = float(p[3])
    by = float(p[4])
    bz = float(p[5])
    n = int(p[6])

    return ax, ay, az, bx, by, bz, n


# IOD position format 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
def format_position(ra, de):
    ram = 60.0 * ra / 15.0
    rah = int(np.floor(ram / 60.0))
    ram -= 60.0 * rah

    des = np.sign(de)
    dem = 60.0 * np.abs(de)
    ded = int(np.floor(dem / 60.0))
    dem -= 60.0 * ded

    if des == -1:
        sign = "-"
    else:
        sign = "+"

    return ("%02d%06.3f%c%02d%05.2f" % (rah, ram, sign, ded, dem)).replace(".", "")


# Inside selection
def inside_selection_area(tmin, tmax, x0, y0, dxdt, dydt, x, y, dt=2.0, w=10.0):
    dx, dy = x - x0, y - y0
    ang = -np.arctan2(dy, dx)
    r = np.sqrt(dx**2 + dy**2)
    drdt = r / (tmax - tmin)
    sa, ca = np.sin(ang), np.cos(ang)
    tmid = 0.5 * (tmin + tmax)

    xmid = x0 + dxdt * tmid
    ymid = y0 + dydt * tmid

    dx, dy = x0 - xmid, y0 - ymid
    rm = ca * dx - sa * dy
    wm = sa * dx + ca * dy
    dtm = rm / drdt

    print(">> ", wm, dtm)

    if (np.abs(wm) < w) & (np.abs(dtm) < dt):
        return True
    else:
        return False


# Angular offsets from spherical angles
def deproject(l0, b0, l, b):
    lt = l * np.pi / 180
    bt = b * np.pi / 180
    l0t = l0 * np.pi / 180
    b0t = b0 * np.pi / 180

    # To vector
    r = np.array([np.cos(lt) * np.cos(bt), np.sin(lt) * np.cos(bt), np.sin(bt)])

    # Rotation matrices
    cl, sl = np.cos(l0t), np.sin(l0t)
    Rl = np.array([[cl, sl, 0], [-sl, cl, 0], [0, 0, 1]])
    cb, sb = np.cos(b0t), np.sin(b0t)
    Rb = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])

    # Apply rotations
    r = Rl.dot(r)
    r = Rb.dot(r)

    # Back to angles
    radius = np.arccos(r[0])
    position_angle = np.arctan2(r[1], r[2])

    # To offsets
    dl, db = radius * np.sin(position_angle), radius * np.cos(position_angle)

    return dl * 180 / np.pi, db * 180 / np.pi


def angle_difference(ang1, ang2):
    x1, y1 = np.cos(ang1), np.sin(ang1)
    x2, y2 = np.cos(ang2), np.sin(ang2)

    return np.arccos(x1 * x2 + y1 * y2)


def correct_bool_state(c):
    # Return on no changes
    if np.all(c):
        return c
    if np.all(~c):
        return c

    # Find indices of first false to true flip
    idx = np.argwhere(c)[0].squeeze()

    # Decrement index and keep in range
    idx = max(0, idx - 1)

    # Flip bool at that index
    c[idx] = True

    return c


def correct_stationary_coordinates(tmid, t, p, direction=1):
    # Compute LST
    t.delta_ut1_utc = 0
    tmid.delta_ut1_utc = 0
    h = t.sidereal_time("mean", longitude=0.0).degree
    hmid = tmid.sidereal_time("mean", longitude=0.0).degree

    dra = direction * (hmid - h)
    
    # Transform to epoch of date
    pt = p.transform_to(FK5(equinox=tmid))

    # Apply correction
    ra = pt.ra.degree + dra
    dec = pt.dec.degree

    # Transform to J2000
    pt = SkyCoord(ra=ra, dec=dec, unit="deg", frame="fk5", equinox=tmid).transform_to(FK5(equinox="J2000"))
    
    return pt
