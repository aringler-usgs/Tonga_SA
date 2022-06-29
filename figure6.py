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
mpl.rc('font', size=14)



def get_upo(lat, lon, depths):
    model_result = model.get_point(lat, lon)
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



depths = np.linspace(0,10000., 1000)
fig = plt.figure(1, figsize=(12,9))
zrats, upos, zcorrs =[],[], []
rrats, wpos, rcorrs = [], [], []
zratsb, uposb, zcorrsb =[],[],[]
f = open('Theresults')
for line in f:
    line = line.split(', ')
    lat, lon = float(line[1]), float(line[2])
    rcorr = float(line[3])
    rrat = float(line[4])*10**9
    zcorr = float(line[5])
    zrat = float(line[6])*10**9
    vp, vs, rho = get_upo(lat, lon, depths)
    try:
        vp = np.mean(vp[depths <= 3500.])
        vs = np.mean(vs[depths <= 3500.])
        rho =np.mean(rho[depths <= 3500.])
    except:
        continue
    mu = rho*vs**2
    lam = rho*vp**2 - 2*mu
    # eqn 23 of sorrells
    c0=340.
    # # return in units of nm/s/Pa
    #upo = (c0/(2*(lam+mu)))*10**9
    upo = c0*(lam + 2*mu)/(2*mu*(lam+mu))*10**9
    print(line[0])
    mval = upo
    #mval = np.mean(upo[depths <= 7500.])
    upo = mval


    if upo >= 20:
        upo=20.
    #if zrat >=30:
    #    zrat = 30.
    #if wpo > 100:
    #    wpo=100.
    if rrat >=100:
        rrat = 100.
    if zrat >= 20:
        zrat = 20.
        #continue

    #if zrat > 20:
    if zcorr >= 0.8:
        print(line[0] + ' ' + str(zrat) + ' '  + str(upo))

    if (zcorr >= 0.8) and (upo < 20.):
        zrats.append(zrat)
        upos.append(upo)
        zcorrs.append(zcorr)
    else:
        zratsb.append(zrat)
        uposb.append(upo)
        zcorrsb.append(zcorr)
    if rcorr >= 0.8:
        rrats.append(rrat)
        rcorrs.append(rcorr)
    #    wpos.append(wpo)

pzv = np.polyfit(upos, zrats,1)
print(pzv)
pz = np.poly1d(pzv)
#prv= np.polyfit(wpos, rrats,1)
#pr = np.poly1d(prv)
#plt.subplot(1,2,1)
ec=plt.scatter(upos, zrats, c=zcorrs, marker='o', vmin=0, vmax=1., label='Used in Fit')
plt.scatter(uposb, zratsb, c=zcorrsb, marker='s', vmin=0, vmax=1., alpha=0.5, label='Not Used in Fit')
plt.plot(np.linspace(0,90,100), pz(np.linspace(0,90,100)), color='r', label= 'Slope=' + str(round(pzv[0],3)) )

plt.plot(np.linspace(0,90,100),np.linspace(0,90,100), label= 'Slope=1', color='k')
plt.xlabel('Vertical Ratio Crust 1.0 (nm/s/Pa)')
plt.ylabel('Vertical Ratio Hunga Tonga (nm/s/Pa)')
plt.xlim((0.,21.))
plt.ylim(0.,21.)
plt.legend(loc='lower right', ncol=2)
cbar = plt.colorbar(ec)
cbar.set_label('Correlation Vertial to Hilbert Transform of Pressure')
# plt.subplot(1,2,2)
# plt.scatter(wpos, rrats, c=rcorrs)
# plt.plot(np.linspace(0,30,100), pr(np.linspace(0,30,100)), color='r')
# plt.xlabel('Radial Ratio Crust 1.0 (nm/s/Pa)')
# plt.ylabel('Radial Ratio Hunga Tonga (nm/s/Pa)')
# plt.xlim((0.,100.))
# plt.ylim(0.,100.)
# #plt.tight_layout()
plt.savefig('Figure6.png', format='PNG', dpi=400)
plt.savefig('Figure6.pdf', format='PDF', dpi=400)
f.close()

