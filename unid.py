#!/usr/bin/env python
import glob,os
import numpy as np
from stio import fourframe
import ppgplot
from tempfile import mkstemp
import subprocess

def find_tracks(fname,trksig,ntrkmin):
    # Read four frame
    ff=fourframe(fname)

    # Skip saturated frames
    if np.sum(ff.zavg>240.0)/float(ff.nx*ff.ny)>0.95:
        return

    # Extract significant pixels
    x,y,t,sig=ff.significant_pixels(trksig)

    # Mask image
    zsig=(ff.zmax-ff.zavg)/(ff.zstd+1e-9)
    c=zsig<trksig
    zsig[c]=0.0
    
    # Save points
    f,fpath=mkstemp()
    with open(fpath,"w") as f:
        for i in range(len(t)):
            f.write("%f,%f,%f\n"%(x[i],y[i],t[i]))
    f.close()
    
    # Run 3D Hough line-finding algorithm
    command="/home/bassa/temp/hough3d/hough3d-code/hough3dlines -dx 4 -minvotes %d %s"%(ntrkmin,fpath)
    output=str(subprocess.check_output(command,shell=True,stderr=subprocess.STDOUT),'utf-8')

    # Remove file
    os.remove(fpath)

    # Skip if no output
    if len(output)==0:
        return
    
    # Unpack ouput
    s=output.split()
    x0,y0,t0=float(s[1]),float(s[2]),float(s[3])
    dx,dy,dt=float(s[4]),float(s[5]),float(s[6])
    print(output)
    
    # ppgplot arrays
    tr=np.array([-0.5,1.0,0.0,-0.5,0.0,1.0])
    heat_l=np.array([0.0,0.2,0.4,0.6,1.0])
    heat_r=np.array([0.0,0.5,1.0,1.0,1.0])
    heat_g=np.array([0.0,0.0,0.5,1.0,1.0])
    heat_b=np.array([0.0,0.0,0.0,0.3,1.0])

    ppgplot.pgopen("/xs")
    ppgplot.pgpap(0.0,1.0)
    ppgplot.pgsvp(0.1,0.95,0.1,0.8)
    ppgplot.pgsfs(2)
            
    ppgplot.pgsch(0.8)
    ppgplot.pgmtxt("T",6.0,0.0,0.0,"UT Date: %.23s  COSPAR ID: %04d"%(ff.nfd,ff.site_id))
    if (3600.0*ff.crres[0]<1e-3) | (3600.0*ff.crres[1]<1e-3) | (ff.crres[0]/ff.sx>2.0) | (ff.crres[1]/ff.sy>2.0):
        ppgplot.pgsci(2)
    else:
        ppgplot.pgsci(1)
    ppgplot.pgmtxt("T",4.8,0.0,0.0,"R.A.: %10.5f (%4.1f'') Decl.: %10.5f (%4.1f'')"%(ff.crval[0],3600.0*ff.crres[0],ff.crval[1],3600.0*ff.crres[1]))
    ppgplot.pgsci(1)
    ppgplot.pgmtxt("T",3.6,0.0,0.0,"FoV: %.2f\\(2218)x%.2f\\(2218) Scale: %.2f''x%.2f'' pix\\u-1\\d"%(ff.wx,ff.wy,3600.0*ff.sx,3600.0*ff.sy))
    ppgplot.pgmtxt("T",2.4,0.0,0.0,"Stat: %5.1f+-%.1f (%.1f-%.1f)"%(np.mean(ff.zmax),np.std(ff.zmax),ff.zmaxmin,ff.zmaxmax))
        
    ppgplot.pgsch(1.0)
    ppgplot.pgwnad(0.0,ff.nx,0.0,ff.ny)
    ppgplot.pglab("x (pix)","y (pix)"," ")
    ppgplot.pgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5)
            
    ppgplot.pgimag(zsig,ff.nx,ff.ny,0,ff.nx-1,0,ff.ny-1,10.0,0.0,tr)
    ppgplot.pgbox("BCTSNI",0.,0,"BCTSNI",0.,0)
    ppgplot.pgstbg(1)

    ppgplot.pgsci(4)
    t=np.linspace(0.0,ff.texp)
    x=x0+(t-t0)*dx/dt
    y=y0+(t-t0)*dy/dt
    ppgplot.pgline(x,y)
    ppgplot.pgpt(np.array([x0]),np.array([y0]),4)

    ppgplot.pgend()

    
if __name__ == '__main__':
    # Track selection sigma
    trksig=4.0

    # Minimum track points
    ntrkmin=10

    # Get files
    files=sorted(glob.glob("2*.fits"))
    #find_tracks("2018-03-14T19:08:43.568.fits",trksig,ntrkmin)
    find_tracks("2018-03-14T19:24:04.127.fits",trksig,ntrkmin)
    
    # Process files
#    for file in files:
#        find_tracks(file,trksig,ntrkmin)
