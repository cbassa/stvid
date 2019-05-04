#!/usr/bin/env python
from __future__ import print_function
import os
import glob
import shutil
from stvid.stio import fourframe, satid, observation
from stvid.astrometry import is_calibrated
import numpy as np
import ppgplot as ppg
from scipy import optimize, ndimage
from termcolor import colored
import datetime

# Gaussian model
def model(a, nx, ny):
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    dx, dy = (x-a[0])/a[2], (y-a[1])/a[2]
    arg = -0.5*(dx**2+dy**2)
    return a[3]*np.exp(arg)+a[4]

# Residual function
def residual(a, img):
    ny, nx = img.shape
    mod = model(a, nx, ny)
    return (img-mod).ravel()

# Find peak
def peakfind(img, w=1.0):
    # Find approximate location
    ny, nx = img.shape
    i = np.argmax(img)
    y0 = int(i/nx)
    x0 = i-y0*nx

    # Image properties
    imgavg = np.mean(img)
    imgstd = np.std(img)

    # Estimate
    a = np.array([x0, y0, w, img[y0, x0]-imgavg, imgavg])
    q, cov_q, infodict, mesg, ier = optimize.leastsq(residual,
                                                     a,
                                                     args=(img),
                                                     full_output=1)

    # Extract
    xc, yc, w = q[0], q[1], q[2]

    # Significance
    sigma = (a[3]-imgavg)/(imgstd+1e-9)

    return xc, yc, w, sigma

# Plot selection
def plot_selection(id, x0, y0, dt=2.0, w=10.0):
    dx, dy = id.x1-id.x0, id.y1-id.y0
    ang = np.arctan2(dy, dx)
    r = np.sqrt(dx**2+dy**2)
    drdt = r/(id.t1-id.t0)
    sa, ca = np.sin(ang), np.cos(ang)

    dx = np.array([-dt, -dt, dt, dt, -dt])*drdt
    dy = np.array([w, -w, -w, w, w])
    x = ca*dx-sa*dy+x0
    y = sa*dx+ca*dy+y0

    ppg.pgsci(7)
    ppg.pgline(x, y)

    return

# Check if point is inside selection
def inside_selection(ident, tmid, x0, y0, dt=2.0, w=10.0):
    dx, dy = ident.x1-ident.x0, ident.y1-ident.y0
    ang = -np.arctan2(dy, dx)
    r = np.sqrt(dx**2+dy**2)
    drdt = r/(ident.t1-ident.t0)
    sa, ca = np.sin(ang), np.cos(ang)

    xmid = ident.x0+ident.dxdt*tmid
    ymid = ident.y0+ident.dydt*tmid
    
    dx, dy = x0-xmid, y0-ymid
    rm = ca*dx-sa*dy
    wm = sa*dx+ca*dy
    dtm = rm/drdt

    if (abs(wm) < w) & (abs(dtm) < dt):
        return True
    else:
        return False

# Get COSPAR ID
def get_cospar(norad, nfd):
    f = open(os.path.join(os.getenv("ST_DATADIR"), "data/desig.txt"))
    lines = f.readlines()
    f.close()
    
    try:
        cospar = ([line for line in lines if str(norad) in line])[0].split()[1]
    except IndexError:
        t = datetime.datetime.strptime(nfd[:-4], "%Y-%m-%dT%H:%M:%S")
        doy = int(t.strftime("%y%j"))+500
        cospar = "%sA"%doy

    return "%2s %s" % (cospar[0:2], cospar[2:])

# IOD position format 2: RA/DEC = HHMMmmm+DDMMmm MX   (MX in minutes of arc)
def format_position(ra, de):
    ram = 60.0*ra/15.0
    rah = int(np.floor(ram/60.0))
    ram -= 60.0*rah

    des = np.sign(de)
    dem = 60.0*np.abs(de)
    ded = int(np.floor(dem/60.0))
    dem -= 60.0*ded

    if des == -1:
        sign = "-"
    else:
        sign = "+"

    return ("%02d%06.3f%c%02d%05.2f"
            % (rah, ram, sign, ded, dem)).replace(".", "")

# Format IOD line
def format_iod_line(norad, cospar, site_id, t, ra, de):
    pstr = format_position(ra, de)
    tstr = t.replace("-", "") \
            .replace("T", "") \
            .replace(":", "") \
            .replace(".", "")

    return "%05d %-9s %04d G %s 17 25 %s 37 S" % (norad,
                                                  cospar,
                                                  site_id,
                                                  tstr,
                                                  pstr)

def store_results(ident, fname, path, iod_line):
    # Find destination
    if ident.catalog.find("classfd.tle") > 0:
        outfname = os.path.join(path, "classfd/classfd.dat")
        dest = os.path.join(path, "classfd")
        color = "blue"
    elif ident.catalog.find("inttles.tle") > 0:
        outfname = os.path.join(path, "classfd/classfd.dat")
        dest = os.path.join(path, "classfd")
        color = "blue"
    elif ident.catalog.find("catalog.tle") > 0:
        outfname = os.path.join(path, "catalog/catalog.dat")
        dest = os.path.join(path, "catalog")
        color = "grey"
    else:
        dest = os.path.join(path, "unid")
        outfname = os.path.join(path, "unid/unid.dat")
        color = "magenta"

    # Print iod line
    print(colored(iod_line, color))
        
    # Copy files
    shutil.copy2(fname, dest)
    shutil.copy2(fname + ".cat", dest)
    shutil.copy2(fname + ".id", dest)
    shutil.copy2(fname + ".png", dest)
    try:
        shutil.move(fname.replace(".fits", "_%05d.png"%ident.norad), dest)
    except:
        pass
        
    # Write iodline
    fp = open(outfname, "a")
    fp.write("%s\n" % iod_line)
    fp.close()
    
    return

def plot_header(fname, ff, iod_line):
    # ppgplot arrays
    heat_l = np.array([0.0, 0.2, 0.4, 0.6, 1.0])
    heat_r = np.array([0.0, 0.5, 1.0, 1.0, 1.0])
    heat_g = np.array([0.0, 0.0, 0.5, 1.0, 1.0])
    heat_b = np.array([0.0, 0.0, 0.0, 0.3, 1.0])
    
    # Plot
    ppg.pgopen(fname)
    ppg.pgpap(0.0, 1.0)
    ppg.pgsvp(0.1, 0.95, 0.1, 0.8)

    ppg.pgsch(0.8)
    ppg.pgmtxt("T", 6.0, 0.0, 0.0,
               "UT Date: %.23s  COSPAR ID: %04d"
               % (ff.nfd, ff.site_id))
    if is_calibrated(ff):
        ppg.pgsci(1)
    else:
        ppg.pgsci(2)
    ppg.pgmtxt("T", 4.8, 0.0, 0.0,
               "R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')"
               % (ff.crval[0],
                  3600.0*ff.crres[0],
                  ff.crval[1],
                  3600.0*ff.crres[1]))
    ppg.pgsci(1)
    ppg.pgmtxt("T", 3.6, 0.0, 0.0,
               ("FoV: %.2f\\(2218)x%.2f\\(2218) "
                "Scale: %.2f''x%.2f'' pix\\u-1\\d")
               % (ff.wx, ff.wy, 3600.0*ff.sx, 3600.0*ff.sy))
    ppg.pgmtxt("T", 2.4, 0.0, 0.0,
               "Stat: %5.1f+-%.1f (%.1f-%.1f)"
               % (np.mean(ff.zmax),
                  np.std(ff.zmax),
                  ff.zmaxmin,
                  ff.zmaxmax))
    ppg.pgmtxt("T", 0.3, 0.0, 0.0, iod_line)

    ppg.pgsch(1.0)
    ppg.pgwnad(0.0, ff.nx, 0.0, ff.ny)
    ppg.pglab("x (pix)", "y (pix)", " ")
    ppg.pgctab(heat_l, heat_r, heat_g, heat_b, 5, 1.0, 0.5)


# Extract tracks
def extract_tracks(fname, trkrmin, drdtmin, trksig, ntrkmin, path):
    # Read four frame
    ff = fourframe(fname)

    # Skip saturated frames
    if np.sum(ff.zavg > 240.0)/float(ff.nx*ff.ny) > 0.95:
        return

    # Read satelite IDs
    try:
        f = open(fname+".id")
    except OSError:
        print("Cannot open", fname+".id")
    else:
        lines = f.readlines()
        f.close()

    tr = np.array([-0.5, 1.0, 0.0, -0.5, 0.0, 1.0])

    # Parse identifications
    idents = [satid(line) for line in lines]
    
    # Identify unknowns
    for ident0 in idents:
        if ident0.catalog == "unidentified":
            for ident1 in idents:
                if ident1.catalog == "unidentified":
                    continue

                # Find matches
                p1 = inside_selection(ident1, ident0.t0, ident0.x0, ident0.y0)
                p2 = inside_selection(ident1, ident0.t1, ident0.x1, ident0.y1)

                # Match found
                if p1 and p2:
                    # Copy info
                    ident0.norad = ident1.norad
                    ident0.catalog = ident1.catalog
                    ident0.state = ident1.state
                    ident1.state = "remove"
                    break
    
    # Loop over identifications
    for ident in idents:
        # Skip superseded unknowns
        if ident.state == "remove":
            continue
        
        # Skip slow moving objects
        drdt = np.sqrt(ident.dxdt**2+ident.dydt**2)
        if drdt < drdtmin:
            continue

        # Extract significant pixels along a track
        x, y, t, sig = ff.significant_pixels_along_track(trksig,
                                                         ident.x0,
                                                         ident.y0,
                                                         ident.dxdt,
                                                         ident.dydt,
                                                         trkrmin)

        # Fit tracks
        if len(t) > ntrkmin:
            # Get times
            tmin = np.min(t)
            tmax = np.max(t)
            tmid = 0.5*(tmax+tmin)
            mjd = ff.mjd+tmid/86400.0

            # Skip if no variance in time
            if np.std(t-tmid) == 0.0:
                continue

            # Very simple polynomial fit; no weighting, no cleaning
            px = np.polyfit(t-tmid, x, 1)
            py = np.polyfit(t-tmid, y, 1)

            # Extract results
            x0, y0 = px[1], py[1]
            dxdt, dydt = px[0], py[0]
            xmin = x0+dxdt*(tmin-tmid)
            ymin = y0+dydt*(tmin-tmid)
            xmax = x0+dxdt*(tmax-tmid)
            ymax = y0+dydt*(tmax-tmid)

            cospar = get_cospar(ident.norad, ff.nfd)
            obs = observation(ff, mjd, x0, y0)
            iod_line = "%s" % format_iod_line(ident.norad,
                                              cospar,
                                              ff.site_id,
                                              obs.nfd,
                                              obs.ra,
                                              obs.de)

            # Create diagnostic plot
            plot_header(fname.replace(".fits", "_%05d.png/png" % ident.norad), ff, iod_line)

            ppg.pgimag(ff.zmax, ff.nx, ff.ny, 0, ff.nx-1,
                       0, ff.ny-1, ff.zmaxmax, ff.zmaxmin, tr)
            ppg.pgbox("BCTSNI", 0., 0, "BCTSNI", 0., 0)
            ppg.pgstbg(1)
            
            ppg.pgsci(0)
            if ident.catalog.find("classfd.tle") > 0:
                ppg.pgsci(4)
            elif ident.catalog.find("inttles.tle") > 0:
                ppg.pgsci(3)

            ppg.pgpt(np.array([x0]), np.array([y0]), 4)
            ppg.pgmove(xmin, ymin)
            ppg.pgdraw(xmax, ymax)
            ppg.pgsch(0.65)
            ppg.pgtext(np.array([x0]), np.array([y0]), " %05d" % ident.norad)
            ppg.pgsch(1.0)
            ppg.pgsci(1)
            
            ppg.pgend()

            # Store results
            store_results(ident, fname, path, iod_line)

        elif ident.catalog.find("classfd.tle") > 0:
            # Track and stack
            t = np.linspace(0.0, ff.texp)
            x, y = ident.x0+ident.dxdt*t, ident.y0+ident.dydt*t
            c = (x > 0) & (x < ff.nx) & (y > 0) & (y < ff.ny)

            # Skip if no points selected
            if np.sum(c) == 0:
                continue

            # Compute track
            tmid = np.mean(t[c])
            mjd = ff.mjd+tmid/86400.0
            xmid = ident.x0+ident.dxdt*tmid
            ymid = ident.y0+ident.dydt*tmid
            ztrk = ndimage.gaussian_filter(ff.track(ident.dxdt, ident.dydt, tmid),
                                           1.0)
            vmin = np.mean(ztrk)-2.0*np.std(ztrk)
            vmax = np.mean(ztrk)+6.0*np.std(ztrk)

            # Select region
            xmin = int(xmid-100)
            xmax = int(xmid+100)
            ymin = int(ymid-100)
            ymax = int(ymid+100)
            if xmin < 0:
                xmin = 0
            if ymin < 0:
                ymin = 0
            if xmax > ff.nx:
                xmax = ff.nx-1
            if ymax > ff.ny:
                ymax = ff.ny-1

            # Find peak
            x0, y0, w, sigma = peakfind(ztrk[ymin:ymax, xmin:xmax])
            x0 += xmin
            y0 += ymin

            # Skip if peak is not significant
            if sigma < trksig:
                continue

            # Skip if point is outside selection area
            if inside_selection(ident, tmid, x0, y0) is False:
                continue

            # Format IOD line
            cospar = get_cospar(ident.norad, ff.nfd)
            obs = observation(ff, mjd, x0, y0)
            iod_line = "%s" % format_iod_line(ident.norad,
                                              cospar,
                                              ff.site_id,
                                              obs.nfd,
                                              obs.ra,
                                              obs.de)


            # Create diagnostic plot
            plot_header(fname.replace(".fits", "_%05d.png/png" % ident.norad), ff, iod_line)


            ppg.pgimag(ztrk, ff.nx, ff.ny, 0, ff.nx-1,
                       0, ff.ny-1, vmax, vmin, tr)
            ppg.pgbox("BCTSNI", 0., 0, "BCTSNI", 0., 0)
            ppg.pgstbg(1)

            plot_selection(ident, xmid, ymid)

            ppg.pgsci(0)
            if ident.catalog.find("classfd.tle") > 0:
                ppg.pgsci(4)
            elif ident.catalog.find("inttles.tle") > 0:
                ppg.pgsci(3)
            ppg.pgpt(np.array([ident.x0]), np.array([ident.y0]), 17)
            ppg.pgmove(ident.x0, ident.y0)
            ppg.pgdraw(ident.x1, ident.y1)
            ppg.pgpt(np.array([x0]), np.array([y0]), 4)
            ppg.pgsch(0.65)
            ppg.pgtext(np.array([ident.x0]),
                       np.array([ident.y0]),
                       " %05d" % ident.norad)
            ppg.pgsch(1.0)
            ppg.pgsci(1)

            ppg.pgend()

            # Store results
            store_results(ident, fname, path, iod_line)
