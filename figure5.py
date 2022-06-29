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
    #upo = (c0/(2*(lam+mu)))*10**9
    mval = upo
    mval = np.mean(upo[depths <= 3500.])
    ax[0].plot(vp/1000.,depths/1000., label=sta + ' $u_z /P_0 =$' + str(round(mval,2)) + ' nm/s/Pa', linewidth=2)
    ax[0].axhspan(0,3.5,color='silver', alpha=0.3)
    ax[0].set_ylabel('Depth (km)')
    ax[0].set_ylim((min(depths/1000.), max(depths/1000.)))
    ax[0].set_xlabel('$V_P$ ($km/s$)')
    ax[0].invert_yaxis()


    ax[1].plot(vs/1000.,depths/1000., linewidth=2)
    ax[1].axhspan(0,3.5,color='silver', alpha=0.3)
    ax[1].set_ylim((min(depths/1000.), max(depths/1000.)))
    ax[1].invert_yaxis()
    ax[1].set_xlabel('$V_S$ ($km/s$)')

    ax[2].plot(rho/1000., depths/1000., linewidth=2)
    ax[2].axhspan(0,3.5,color='silver', alpha=0.3)
    ax[2].set_ylim((min(depths/1000.), max(depths/1000.)))
    ax[2].invert_yaxis()
    ax[2].set_xlabel(r'$\rho$ ($Kg/m^3$)')

    ax[3].plot(upo, depths/1000., linewidth=2)
    ax[3].axhspan(0,3.5,color='silver', alpha=0.3)
    ax[3].set_ylim((min(depths/1000.), max(depths/1000.)))
    ax[3].set_xlabel('Vertical Ratio ($nm/s/Pa$)')
    ax[3].invert_yaxis()
    f = open(sta + '_data','w')
    f.write('depth, vs, vp, rho, upo\n')
    for trip in zip(depths, vp, vs, rho, upo):
        f.write(str(trip[0]/1000.) + ', ' + str(trip[1]/1000.) + ', ' + str(trip[2]/1000.) + ', ' + str(trip[3]/1000.)  + ', ' + str(trip[4]) + '\n')
    f.close()



ax[0].text(1.1, -0.3, '(a)')
ax[1].text(-0.2, -0.3, '(b)')
ax[2].text(0.7, -0.3, '(c)')
ax[3].text(-1, -0.3, '(d)')
#fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.12)
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=4)
plt.savefig('Figure5.PNG', format='PNG', dpi=400)
plt.savefig('Figure5.PDF', format='PDF', dpi=400)







# upos, vps, vss, rhos =[], [], [], []
# for idx, sta in enumerate(stas):
#     lat, lon = lats[idx], lons[idx]
#     for depth in depths:
#         upo, vp, vs, rho = get_upo(lat, lon, depth)
#         upos.append(upo)
#         vps.append(vp/1000.)
#         vss.append(vs)
#         rhos.append(rho)


#     plt.subplot(4,1,1)
#     plt.plot(vps, depths)
#     # plt.subplot(4,1,2)
#     # plt.plot(vss, depths)
#     # plt.subplot(4,1,3)
#     # plt.plot(rhos, depths)
#     # plt.subplot(4,1,4)
#     # plt.plot(upos, depths)

# plt.show()






