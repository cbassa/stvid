#!/usr/bin/env python3
import os
import glob

import configparser

from stvid.fourframe import FourFrame
from stvid.fourframe import Observation

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms

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
    
def plot_prediction(p, ax, tlefiles, colors, dt=2.0, w=10.0):
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

    
    if ff.in_frame(xs[0], ys[0]):
        ax.text(xs[0], ys[0], f" {p.satno:05d} ", color=color, ha=ha, in_layout=False)
    
    for state, linestyle in zip(["sunlit", "umbra", "eclipsed"], ["solid", "dashed", "dotted"]):
        c = correct_bool_state(p.state == state)
        ax.plot(np.ma.masked_where(~c, xt), np.ma.masked_where(~c, yt), color=color, linestyle=linestyle)

    return


def deproject(l0, b0, l, b):
    lt = l * np.pi / 180
    bt = b * np.pi / 180
    l0t = l0 * np.pi / 180
    b0t = b0 * np.pi / 180
    
    # To vector
    r = np.array([np.cos(lt) * np.cos(bt),
                  np.sin(lt) * np.cos(bt),
                  np.sin(bt)])

    # Rotation matrices
    cl, sl = np.cos(l0t), np.sin(l0t)
    Rl = np.array([[cl, sl, 0],
                    [-sl, cl, 0],
                    [0, 0, 1]])
    cb, sb = np.cos(b0t), np.sin(b0t)
    Rb = np.array([[cb, 0, sb],
                     [0, 1, 0],
                     [-sb, 0, cb]])    

    # Apply rotations
    r = Rl.dot(r)
    r = Rb.dot(r)

    # Back to angles
    radius = np.arccos(r[0])
    position_angle = np.arctan2(r[1], r[2]) 

    # To offsets
    dl, db = radius * np.sin(position_angle), radius * np.cos(position_angle)

    return dl * 180 / np.pi, db * 180 / np.pi


def residuals(o, p):
    ra0, dec0 = np.mean(o.ra), np.mean(o.dec)

    # Compute offsets
    rxo, ryo = deproject(ra0, dec0, o.ra, o.dec)
    rxp, ryp = deproject(ra0, dec0, p.ra, p.dec)

    # Compute polynomials for prediction
    if len(p.t) > 3:
        px = np.polyfit(p.t, rxp, 2)
        py = np.polyfit(p.t, ryp, 2)
    elif len(p.t) > 2:
        px = np.polyfit(p.t, rxp, 1)
        py = np.polyfit(p.t, ryp, 1)
    else:
        return np.nan

   
    drx = rxo - np.polyval(px, o.t)
    dry = ryo - np.polyval(py, o.t)

    r = np.sqrt(drx**2+dry**2)

    return np.sqrt(np.sum(r**2) / len(r))

def angle_difference(ang1, ang2):
    x1, y1 = np.cos(ang1), np.sin(ang1)
    x2, y2 = np.cos(ang2), np.sin(ang2)

    return np.arccos(x1 * x2 + y1 * y2)

if __name__ == "__main__":
    # Read configuration file
    config_file = "config_new.ini"
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ":"))
    result = cfg.read([config_file])

    # Observer settings
    site_id = cfg.getint("Observer", "cospar")
    
    # Extract colors for TLE files
    colors, abbrevs, tlefiles, catalognames = [], [], [], []
    for key, value in cfg.items("Elements"):
        if "tlefile" in key:
            tlefiles.append(value)
        elif "color" in key:
            colors.append(value)
        elif "name" in key:
            catalognames.append(value)
        elif "abbrev" in key:
            abbrevs.append(value)
    color_detected = cfg.get("LineDetection", "color")

    # Colormap
    cmap = cfg.get("DiagnosticPlot", "colormap")

    # Identification settings
    rm_max = cfg.getfloat("Identification", "max_off_track_offset_deg")
    dtm_max = cfg.getfloat("Identification", "max_along_track_offset_s")
    dpa_max = cfg.getfloat("Identification", "max_direction_difference_deg")
    fdr_max = cfg.getfloat("Identification", "max_velocity_difference_percent")
    
    fname = "/data3/satobs/test/185300/processed/2022-03-24T18:53:20.708.fits"
    fnames = sorted(glob.glob("/data3/satobs/test/185300/processed/2*.fits"))
    fnames = sorted(glob.glob("/data3/satobs/test/asi174mm/2*.fits"))    
    #    fname = "/data3/satobs/test/2022-04-02T21:35:17.038.fits"

    
    for fname in fnames:
        print(fname)
        ff = FourFrame(fname)

        # Output file root
        froot = os.path.splitext(fname)[0]

        # Generate predictions
        predictions = ff.generate_satellite_predictions(cfg)

        # Detect tracks
        tracks = ff.find_tracks_by_hough3d(cfg)

        # Identify tracks and format observations
        obs = []
        for i, t in enumerate(tracks):
            # Default satno
            satno = 90000 + i
            cospar = "22 500A"
            tlefile = None
            for p in predictions:
                # Compute identification constraints
                rx0, ry0, drdt, pa, dr = p.position_and_velocity(t.tmid, t.tmax - t.tmin)
                dtm, rm = p.residual(t.tmid, t.rx0, t.ry0)
                dpa = angle_difference(t.pa, pa) * 180 / np.pi
                fdr = (dr / t.dr - 1) * 100
                if (np.abs(dtm) < dtm_max) & (np.abs(rm) < rm_max) & (np.abs(dpa) < dpa_max) & (np.abs(fdr) < fdr_max):
                    satno = p.satno
                    cospar = p.cospar
                    tlefile = p.tlefile

            t.satno = satno
            t.cospar = cospar
            t.catalogname = "unid"

            # Get catalog abbreviation
            for abbrev, tfile in zip(abbrevs, tlefiles):
                if tfile == tlefile:
                    t.catalogname = abbrev

            # Add to observation
            obs.append(Observation(ff, t.tmid, t.x0, t.y0, site_id, t.satno, t.cospar, t.catalogname))

        # Write observations
        for o in obs:
            iod_line = o.to_iod_line()

            # Open file
            outfname = f"{froot}_{o.satno:05d}_{o.catalogname}.dat"
            with open(outfname, "w") as fp:
                fp.write(f"{iod_line}\n")
            print(iod_line, o.catalogname)

        # Generate plots
        for i in range(len(tracks) + 1):
            if i < len(tracks):
                track = tracks[i]
                iod_line = obs[i].to_iod_line()
                satno = obs[i].satno
                catalogname = obs[i].catalogname
                outfname = f"{froot}_{satno:05d}_{catalogname}.png"
            else:
                track = None
                iod_line = ""
                outfname = f"{fname}.png"
            
            fig, ax = plt.subplots(figsize=(12, 10), dpi=75)

            ax.set_title(f"UT Date: {ff.nfd} COSPAR ID: {ff.site_id}\nR.A.: {ff.crval[0]:10.6f} ({3600 * ff.crres[0]:.1f}\") Decl.: {ff.crval[1]:10.6f} ({3600 * ff.crres[1]:.1f}\")\nFOV: {ff.wx:.2f}$^\circ$x{ff.wy:.2f}$^\circ$ Scale: {3600 * ff.sx:.2f}\"x{3600 * ff.sy:.2f}\" pix$^{{-1}}$\n\n{iod_line}", fontdict={"fontsize": 14, "horizontalalignment": "left"}, loc="left")
    
            ax.imshow(ff.zmax, origin="lower", interpolation="none", vmin=ff.zmaxmin, vmax=ff.zmaxmax,
                      cmap=cmap)

            for p in predictions:
                plot_prediction(p, ax, tlefiles, colors, dt=0)

            if track is not None:
                ax.plot(track.xp, track.yp, color=color_detected, linestyle="-")
                ax.plot(track.x0, track.y0, color=color_detected, marker="o", markerfacecolor="none")
                ax.text(track.x0, track.y0, f" {track.satno:05d}", color=color_detected, ha="center", in_layout=False)
            
            ax.set_xlim(0, ff.nx)
            ax.set_ylim(0, ff.ny)
            ax.set_xlabel("x (pixel)")
            ax.set_ylabel("y (pixel)")
            #        ax.xaxis.set_ticklabels([])
            #        ax.yaxis.set_ticklabels([])


            # Create legend handles
            handles = []
            for catalogname, color in zip(catalognames, colors):
                handles.append(mlines.Line2D([], [], color=color, marker="", label=catalogname))
            handles.append(mlines.Line2D([], [], color=color_detected, marker="o", markerfacecolor="none", label="Detected"))
            for state, linestyle in zip(["Sunlit", "Penumbra", "Eclipsed"], ["solid", "dashed", "dotted"]):
                handles.append(mlines.Line2D([], [], color="k", linestyle=linestyle, marker="", label=state))
            ax.legend(handles=handles, ncol=7, bbox_to_anchor=(0.5, -0.1), loc="center", frameon=True)
                    
            plt.tight_layout()
            #        plt.show()
            plt.savefig(outfname, bbox_inches="tight", pad_inches=0.5)
            plt.close()
