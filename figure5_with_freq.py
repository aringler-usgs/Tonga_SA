#!/usr/bin/env python
import numpy as np
import crust1 
import cartopy.feature as cfeature
import cartopy as cart
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
from scipy.interpolate import griddata
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
model = crust1.crustModel()

import matplotlib as mpl
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font', size=18)

win = 0.5





def get_upo(lat, lon, depths):
    model_result = model.get_point(lat, lon)
    print(model_result)
    cvp, cvs, crho = [],[],[]
    for depth in depths:
        mdepth = 0.
        for layer in model_result:
            if layer == 'water':
                continue
            mdepth += model_result[layer][3]*1000.
            if mdepth >= depth:
                cvp.append(model_result[layer][0]*1000)
                cvs.append(model_result[layer][1]*1000)
                crho.append(model_result[layer][2]*1000)
                break
    cvp = np.array(cvp)
    cvs = np.array(cvs)
    crho = np.array(crho)
    return cvp, cvs, crho


# stations
stas = ['ANMO', 'FUNA', 'HKT', 'QSPA']
lats = [34.94591, -8.5259, 29.96478, -89.9289]
lons =[-106.4572, 179.1966, -95.83812, 144.4382]


fig = plt.figure(1, figsize=(18,12))
fig, ax = plt.subplots(1,4,figsize=(18,12))
ax = ax.flatten()
depths = np.linspace(0,5000., 100)

for idx, sta in enumerate(stas):
    vp, vs, rho = get_upo(lats[idx], lons[idx], depths)
    #try:
    #    vp = np.mean(vp[depths <= 3500.])
    #    vs = np.mean(vs[depths <= 3500.])
    #    rho =np.mean(rho[depths <= 3500.])
    #except:
    #    continue
    mu = rho*vs**2
    lam = rho*vp**2 - 2*mu
    # eqn 23 of sorrells
    c0=340.
    # # return in units of nm/s/Pa
    upo = c0*(lam + 2*mu)/(2*mu*(lam + mu))*10**9
    for freq in [1./25., 1./50., 1./75., 1./100., 1./125. ]:
    

        w= 2*np.pi*freq
        upow = (c0/(2*mu))*((lam + 2*mu)/(lam + mu) - w*(-depths)/c0)*np.exp(w*(-depths)/c0)*10**9


        ax[idx].plot(upow, depths/1000., linewidth=2, label=str(round(freq,3)) + ' Hz')
        ax[idx].axhspan(0,3.5,color='silver', alpha=0.3)
        ax[idx].set_ylim((min(depths/1000.), max(depths/1000.)))
        ax[idx].set_xlabel('Vertical Ratio ($nm/s/Pa$)')
        ax[idx].set_xlim((0,100))
        ax[idx].invert_yaxis()
        if idx == 0:
            ax[idx].set_ylabel('Depth (km)')
        #ax[idx].text(100, 4.8, sta)


ax[0].text(-5, -0.2, '(a) ANMO')
ax[1].text(-5, -0.2, '(b) FUNA')
ax[2].text(-5, -0.2, '(c) HKT')
ax[3].text(-5, -0.2, '(d) QSPA')
fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.12)
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=5)
plt.savefig('Figure5_w_freq.PNG', format='PNG', dpi=400)
plt.savefig('Figure5_w_freq.PDF', format='PDF', dpi=400)







