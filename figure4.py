#!/usr/bin/env python
import matplotlib.pyplot as plt
from scipy.signal import periodogram
import numpy as np
import math
from scipy import signal
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import cartopy.crs as ccrs
import cartopy as cart
import matplotlib as mpl
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font', size=18)


comps = ['Radial to Pressure', 'Vertical to HT of Pressure']
units = ['nm/s/Pa', 'nm/s/Pa']
fig = plt.figure(1, figsize=(12,12))

# We wanna make some maps!

f = open('Theresults','r')
rlatg, rlong, rcorrg, rratg =[],[],[],[] 
rlatb, rlonb, rcorrb, rratb =[],[],[],[] 

vlatg, vlong, vcorrg, vratg =[],[],[],[] 
vlatb, vlonb, vcorrb, vratb =[],[],[],[] 

for line in f:
    line = line.split(', ')
    if float(line[3]) >= 0.8:

        rlatg.append(float(line[1]))
        rlong.append(float(line[2]))
        rcorrg.append(float(line[3]))
        rratg.append(float(line[4])*10**9)
    else:
        rlatb.append(float(line[1]))
        rlonb.append(float(line[2]))
        rcorrb.append(float(line[3]))
        rratb.append(float(line[4])*10**9)

    if float(line[5]) >= 0.8:
        vlatg.append(float(line[1]))
        vlong.append(float(line[2]))
        vcorrg.append(float(line[5]))
        vratg.append(float(line[6])*10**9)
    else:
        vlatb.append(float(line[1]))
        vlonb.append(float(line[2]))
        vcorrb.append(float(line[5]))
        vratb.append(float(line[6])*10**9)

f.close()

ax = fig.add_subplot(2,2,1, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
im = ax.scatter(rlong, rlatg, c=rcorrg, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='o', vmin=-1, vmax=1, edgecolor='k')
im = ax.scatter(rlonb, rlatb, c=rcorrb, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='s', vmin=-1, vmax=1, edgecolor='k')
cbar = plt.colorbar(im, orientation='horizontal')
cbar.set_label('Correlation Radial') 
ax.set_title('(a)',loc='left')
ax = fig.add_subplot(2,2,2, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
im = ax.scatter(rlong, rlatg, c=rratg, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='o', vmin=0, vmax=20, edgecolor='k')
#im = ax.scatter(rlonb, rlatb, c=rratb, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='s', vmin=0, vmax=10, edgecolor='k')
cbar = plt.colorbar(im, orientation='horizontal')
cbar.set_label('Ratio Radial (nm/s/Pa)') 
ax.set_title('(b)',loc='left')
ax = fig.add_subplot(2,2,3, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
im = ax.scatter(vlong, vlatg, c=vcorrg, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='o', vmin=-1, vmax=1, edgecolor='k')
im = ax.scatter(vlonb, vlatb, c=vcorrb, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='s', vmin=-1, vmax=1, edgecolor='k')
cbar = plt.colorbar(im, orientation='horizontal')
cbar.set_label('Correlation Vertical') 
ax.set_title('(c)',loc='left')
ax = fig.add_subplot(2,2,4, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
im = ax.scatter(vlong, vlatg, c=vratg, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='o', vmin=0, vmax=20, edgecolor='k')
#im = ax.scatter(vlonb, vlatb, c=vratb, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='s', vmin=0, vmax=10, edgecolor='k')
cbar = plt.colorbar(im, orientation='horizontal')
cbar.set_label('Ratio Vertical (nm/s/Pa)') 
ax.set_title('(d)',loc='left')





plt.tight_layout()
plt.savefig('Figure4.PNG', format='PNG', dpi=400)
plt.savefig('Figure4.PDF', format='PDF', dpi=400)
plt.close('all') 

