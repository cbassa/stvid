#!/usr/bin/env python
from __future__ import print_function
import subprocess
from stvid.stio import fourframe
from stvid.stio import satid


def generate_satellite_predictions(fname):
    # Format command
    command = "satid %s %s.png/png" % (fname, fname)

    # Run command
    output = subprocess.check_output(command,
                                     shell=True,
                                     stderr=subprocess.STDOUT)

    return output


def find_hough3d_lines(fname, ntrkmin=20, dr=8):
    # Read four frame
    ff = fourframe(fname)

    # Mask frame
    ff.mask(10, 10, 5, 5)

    # Compute selection mask
    x, y, z, t, sig = ff.selection_mask(5.0, 40.0)

    # Skip if not enough points
    if len(t) < 2:
        return []

    # Save points to temporary file
    with open("/tmp/hough.dat", "w") as f:
        for i in range(len(t)):
            f.write("%f,%f,%f\n" % (x[i], y[i], z[i]))
    f.close()

    # Run 3D Hough line-finding algorithm
    command = "hough3dlines -dx %d -minvotes %d %s" % (dr, ntrkmin,
                                                       "/tmp/hough.dat")
    try:
        output = subprocess.check_output(command,
                                         shell=True,
                                         stderr=subprocess.STDOUT)
    except Exception:
        return []

    # Clean output (a bit cluncky)
    cleaned_output = output.decode("utf-8").replace("npoints=", "")
    cleaned_output = cleaned_output.replace(", a=(", " ")
    cleaned_output = cleaned_output.replace("), b=(", " ")
    cleaned_output = cleaned_output.replace(")", "")
    cleaned_output = cleaned_output.replace(",", " ")

    # Generate identifications
    lines = []
    s = cleaned_output.split()
    for i in range(len(s) // 7):
        # Extract points
        x0, y0, z0 = float(s[1 + 7 * i]), float(s[2 + 7 * i]), float(s[3 +
                                                                       7 * i])
        dx, dy, dz = float(s[4 + 7 * i]), float(s[5 + 7 * i]), float(s[6 +
                                                                       7 * i])

        # Reconstruct start and end points
        xmin = x0 - z0 * dx / (dz + 1e-9)
        xmax = x0 + (ff.nz - z0) * dx / (dz + 1e-9)
        ymin = y0 - z0 * dy / (dz + 1e-9)
        ymax = y0 + (ff.nz - z0) * dy / (dz + 1e-9)

        # Output line
        line = "%.23s %8.3f %8.3f %8.3f %8.3f %8.5f  %s unidentified sunlit\n" %\
               (ff.nfd, xmin, ymin, xmax, ymax, ff.texp, 99999-i)
        lines.append(line)

    # Store identifications
    fp = open(fname + ".id", "a")
    for line in lines:
        fp.write(line)
    fp.close()

    return [satid(line) for line in lines]
