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
#mpl.rc('text', usetex=True)
mpl.rc('font', size=18)

win = 0.5


def get_upo(lat, lon, depths):
    model_result = model.get_point(lat, lon)
    cvp, cvs, crho = [],[],[]
    for depth in depths:
        mdepth = 0.
        for layer in model_result:
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







# Eruption Time
stime = UTCDateTime('2022-01-15T04:14:00')
# TUC
etime = UTCDateTime('2022-01-16T00:00:00')

client = Client("IRIS")

stas = []
inv = client.get_stations(network='IU,IC,CU,II', station='*', starttime=stime-60*60*win,
                          endtime=etime+60*60*win, level="response",
                          location='00,31', channel='LDO,LH*,LDI')


for net in inv:
    for sta in net:
        stas.append(sta.code)

s_lats, s_lons, s_stas = [], [], []
for sta in stas:
    # Get the data from IRIS pad time to account for filter ring
    try:
        inv = client.get_stations(network='IU,IC,CU,II', station=sta, starttime=stime-60*60*win,
                              endtime=etime+60*60*win, level="response",
                              location='00,31', channel='LDO,LH*,LDI')


        st = client.get_waveforms(network='IU,IC,CU,II', station=sta, location='00,31',
                              channel='LDO,LH*,LDI', starttime=stime-60*60*win,
                              endtime=etime + win*60*60)
    except:
        continue
    chans = []
    for tr in st:
        chans.append(tr.stats.channel)
    if ('LDO' in chans) or ('LDI' in chans):
        pass
    else:
        continue 

    for tr in st:
        if tr.stats.channel == 'LDO':
            chan_type='LDO'
            break
        if tr.stats.channel == 'LDI':
            chan_type = 'LDI'


    coors = inv.get_coordinates(st[0].id)
    s_lats.append(coors['latitude'])
    s_lons.append(coors['longitude'])
    s_stas.append(st[0].stats.station)

lats = np.arange(89.5, -90, -1)
lons = np.arange(-179.5, 180, 1)

depths = np.linspace(0,6000., 1000)
vals, lonv, latv =[], [], []
for lon in lons:
    for lat in lats:
        vp, vs, rho = get_upo(lat, lon, depths)
        try:
            vp = np.mean(vp[depths <= 5200.])
            vs = np.mean(vs[depths <= 5200.])
            rho =np.mean(rho[depths <= 5200.])
        except:
            continue
        mu = rho*vs**2
        lam = rho*vp**2 - 2*mu
        # eqn 23 of sorrells
        c0=330.
        # # return in units of nm/s/Pa
        upo = (c0*(lam + 2*mu)/(2*mu*(lam+mu)))*10**9
        #upo = (c0/(2*(lam+mu)))*10**9
        mval = upo
        #mval = np.mean(upo[depths <= 3500.])
        vals.append(mval)
        lonv.append(lon)
        latv.append(lat)
        #except:
        #    continue

print(min(vals))
print(max(vals))
fig = plt.figure(1, figsize=(16,12)) 
ax = fig.add_subplot(1,1,1, projection = ccrs.Robinson())

ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines(zorder=20)
im=ax.scatter(lonv, latv, c=vals, transform=ccrs.Geodetic(), zorder=30, vmin=0., vmax=35.)
ax.scatter(s_lons, s_lats, c='white', transform=ccrs.Geodetic(), zorder=30, marker='^', s=200, label='GSN Stations with Pressure', edgecolors='k')
for sta, lat, lon in zip(s_stas, s_lats, s_lons):
    ax.text(lon, lat, sta, transform=ccrs.Geodetic(), fontsize=14, zorder=35)
ax.scatter(175.4, -20.5, c='r',marker='*',s= 300, transform=ccrs.Geodetic(), zorder=30, vmin=0., vmax=10, label='Hunga Tonga Eruption')
ax.scatter(180-175.4, 20.5, c='r',marker='p',s= 230, transform=ccrs.Geodetic(), zorder=30, vmin=0., vmax=10, label='Eruption Antipode')
#cb_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
cbar = fig.colorbar(im,  orientation='horizontal')
cbar.set_label('Crust1.0 Vertical Pressure Ratio (nm/s/Pa)') 
fig.legend(ncol=3, fontsize=18, loc='lower center')

plt.savefig('Figure1.png', format='PNG', dpi=300)
plt.savefig('Figure1.pdf', format='PDF', dpi=200)