#!/usr/bin/env python

import sys,os,glob
from stio import fourframe,satid,observation
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ff=fourframe("2018-02-26T03:09:05.866.fits")

zsig=(ff.zmax-ff.zavg)/(ff.zstd+1e-9)
sigma=5.0

c=zsig>sigma
print(np.sum(c)/float(ff.nx*ff.ny))
xm,ym=np.meshgrid(np.arange(ff.nx),np.arange(ff.ny))
x,y=np.ravel(xm[c]),np.ravel(ym[c])
inum=np.ravel(ff.znum[c]).astype('int')
sig=np.ravel(zsig[c])
t=np.array([ff.dt[i] for i in inum])

plt.figure()
plt.scatter(x,y,c=t,s=1)
plt.colorbar()
plt.show()

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.scatter(x,y,t,c=t,s=1)
plt.show()
