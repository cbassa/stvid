#!/usr/bin/env python
from __future__ import print_function
import os
import subprocess
from stvid.stio import fourframe
from stvid.stio import satid
from stvid.stio import observation
from tempfile import mkstemp

def generate_satellite_predictions(fname):
    # Format command
    command = "satid %s %s.png/png"%(fname, fname)

    # Run command
    output = subprocess.check_output(command, shell=True,
                                     stderr=subprocess.STDOUT)

    return

def find_hough3d_lines(fname,ntrkmin=20,dr=8):
    # Read four frame
    ff=fourframe(fname)

    # Mask frame
    ff.mask(10,10,5,5)

    # Compute selection mask
    x,y,z,t,sig=ff.selection_mask(5.0,40.0)
    
    # Save points to temporary file
    f,fpath=mkstemp()
    with open(fpath,"w") as f:
        for i in range(len(t)):
            f.write("%f,%f,%f\n"%(x[i],y[i],z[i]))
    f.close()

    # Run 3D Hough line-finding algorithm
    command="hough3dlines -dx %d -minvotes %d %s"%(dr,ntrkmin,fpath)
    output=subprocess.check_output(command,shell=True,stderr=subprocess.STDOUT)

    # Remove file
    os.remove(fpath)

    # Clean output (a bit cluncky)
    cleaned_output = output.replace("npoints=", "")
    cleaned_output = cleaned_output.replace(", a=(", " ")
    cleaned_output = cleaned_output.replace("), b=("," ")
    cleaned_output = cleaned_output.replace(")", "")
    cleaned_output = cleaned_output.replace(",", " ")
    print(cleaned_output)
    # Generate identifications
    lines=[]
    s=cleaned_output.split()
    for i in range(len(s)//7):
        x0,y0,z0=float(s[1+7*i]),float(s[2+7*i]),float(s[3+7*i])
        dx,dy,dz=float(s[4+7*i]),float(s[5+7*i]),float(s[6+7*i])
        xmin=x0-z0*dx/dz
        xmax=x0+(ff.nz-z0)*dx/dz
        ymin=y0-z0*dy/dz
        ymax=y0+(ff.nz-z0)*dy/dz
        line="%.23s %8.3f %8.3f %8.3f %8.3f %8.5f  %s unidentified sunlit"%(ff.nfd,xmin,ymin,xmax,ymax,ff.texp,99999)
        lines.append(line)

    return [satid(line) for line in lines]

