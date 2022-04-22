#!/usr/bin/env python3
import configparser

from stvid.stio import FourFrame
from stvid import predict

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

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
        ax.text(xs[0], ys[0], f" {p.satno:05d} ", color=color, ha=ha)
    
    for state, linestyle in zip(["sunlit", "umbra", "eclipsed"], ["solid", "dashed", "dotted"]):
        c = correct_bool_state(p.state == state)
        ax.plot(np.ma.masked_where(~c, xt), np.ma.masked_where(~c, yt), color=color, linestyle=linestyle)



    
    return

if __name__ == "__main__":
    config_file = "config.ini"

    fname = "/data3/satobs/test/185300/processed/2022-03-24T18:53:20.708.fits"
#    fname = "/data3/satobs/test/2022-04-02T21:35:17.038.fits"
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ":"))
    result = cfg.read([config_file])

    # Extract colors for TLE files
    colors, tlefiles, catalognames = [], [], []
    for key, value in cfg.items("Elements"):
        if "tlefile" in key:
            tlefiles.append(value)
        elif "color" in key:
            colors.append(value)
        elif "name" in key:
            catalognames.append(value)

    ff = FourFrame(fname)

    predictions = predict.generate_satellite_predictions(ff, cfg)

    fig, ax = plt.subplots(figsize=(15, 10), dpi=75)

    ax.set_title(f"UT Date: {ff.nfd} COSPAR ID: {ff.site_id}\nR.A.: {ff.crval[0]:10.6f} ({3600 * ff.crres[0]:.1f}\") Decl.: {ff.crval[1]:10.6f} ({3600 * ff.crres[1]:.1f}\")\nFOV: {ff.wx:.2f}$^\circ$x{ff.wy:.2f}$^\circ$ Scale: {3600 * ff.sx:.2f}\"x{3600 * ff.sy:.2f}\" pix$^{{-1}}$", fontdict={"fontsize": 14, "horizontalalignment": "left"}, loc="left")
    
    ax.imshow(ff.zmax, origin="lower", interpolation="none", vmin=ff.zmaxmin, vmax=ff.zmaxmax,
              cmap="viridis")

    for p in predictions:
        plot_prediction(p, ax, tlefiles, colors, dt=0)

    ax.set_xlim(0, ff.nx)
    ax.set_ylim(0, ff.ny)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])

    # Create legend handles
    handles = []
    for catalogname, color in zip(catalognames, colors):
        handles.append(mlines.Line2D([], [], color=color, marker="", label=catalogname))
    for state, linestyle in zip(["Sunlit", "Penumbra", "Eclipsed"], ["solid", "dashed", "dotted"]):
        handles.append(mlines.Line2D([], [], color="k", linestyle=linestyle, marker="", label=state))
    ax.legend(handles=handles, ncol=6, bbox_to_anchor=(0.5, -0.02), loc="center", frameon=False)
                    
    plt.tight_layout()
    plt.savefig("test.png", bbox_inches="tight")
