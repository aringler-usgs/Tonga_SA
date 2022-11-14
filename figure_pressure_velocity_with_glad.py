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



depths = np.linspace(0,7000., 1000)
fig = plt.figure(1, figsize=(16,10))
zrats, upos, zcorrs =[],[], []
rrats, wpos, rcorrs = [], [], []
zratsb, uposb, zcorrsb =[],[],[]
f = open('Theresults')
cnt = 0
stas =[]
therats =[]
crus =[]
zrats =[]
lats =[]
lons =[]
vels =[]
hvs = []
rrats =[]

f2 = open('gladvals','r')
for idx, line in enumerate(f):
    line = line.split(', ')

    rcorr = float(line[3])
    rrat = float(line[4])*10**9
    zcorr = float(line[5])
    zrat = float(line[6])*10**9
    rrat = float(line[4])*10**9
    lat = float(line[1])
    lon = float(line[2])
    vel = float(line[-1])
    if zcorr < 0.8:
        continue
    if zrat > 20:
        zrat = 20.
    cnt +=1
    cvp, cvs, crho = get_upo(lat, lon, depths)
    for idx2, dep in enumerate([5200]):
        vp = np.mean(cvp[(depths <= dep)])
        vs = np.mean(cvs[(depths <= dep)])
        rho =np.mean(crho[(depths <= dep)])
        print('Important')
        print(str(vp) + ' ' + str(vs) + ' ' + str(rho))
        blah = f2.readline()
        blah = blah.split(', ')
        vs = float(blah[4])*1000
        vp = float(blah[5])*1000
        rho = float(blah[3])*1000.
        print(str(vp) + ' ' + str(vs) + ' ' + str(rho))
        mu = rho*vs**2
        lam = rho*vp**2 - 2*mu
        # eqn 23 of sorrells
        c0=340.
        # # return in units of nm/s/Pa
        try:
            upo = (c0*(lam + 2*mu)/(2*mu*(lam+mu)))*10**9
            print(upo)
            #if upo > 20:
                #upo = 20.

            therat = np.abs(100*(upo-zrat)/zrat)
            if therat >= 100:
                therat = 100.

            stas.append(line[0])
            therats.append(therat)
            lats.append(lat)
            lons.append(lon)
            upos.append(upo)
            zrats.append(zrat)
            rrats.append(rrat)
            vels.append(vel)
            hvs.append(float(rrat)/float(zrat))
            #plt.plot(upo, cnt, '.', color='C' + str(idx2+1), marker='*')  
            #plt.plot(zrat,cnt, '.', color='C0')     
        except:
            continue     

print(therats)

ax = fig.add_subplot(1,1,1, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
ax.scatter(175.4, -20.5, c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
im = ax.scatter(lons, lats, c=therats, s = 100., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='o', vmin=0, vmax=50, edgecolor='k')
#im = ax.scatter(rlonb, rlatb, c=rcorrb, s = 50., transform=ccrs.Geodetic(), zorder=3, alpha=0.5, marker='s', vmin=-1, vmax=1, edgecolor='k')
cbar = plt.colorbar(im, orientation='horizontal')
#cbar.set_label('Vertical Ratio Hunga Tonga (nm/s/Pa)') 
cbar.set_label('Vertical Deviation from GLAD25 (%)') 



    #upo = mval

# therats = np.array(therats)
# stas = np.array(stas)
# upos = np.array(upos)
# zrats = np.array(zrats)
# idx = np.argsort(np.array(therats))
# stas = stas[idx]
# therats = therats[idx]
# upos = upos[idx]
# zrats = zrats[idx]


# plt.subplot(1,2,1)
# plt.plot(therats,range(len(stas)),'.', markersize=20)
# plt.yticks(range(len(stas)), stas, fontsize=12)
# plt.xlabel('Deviation from Crust1.0 (%)')
# plt.ylim((min(range(len(stas)))-0.5, max(range(len(stas)))+0.5))
# plt.subplot(1,2,2)
# plt.plot(upos,range(len(stas)), '.', 'C1', markersize=20)
# plt.plot(zrats, range(len(stas)), '.', 'C2', markersize=20)
# for idx, sta in enumerate(stas):
#     plt.plot([zrats[idx],upos[idx]], [idx,idx], 'C3')
# plt.yticks(range(len(stas)), stas, fontsize=12)

# plt.xlabel('Vertical Ratio (nm/s/Pa)')
# plt.ylim((min(range(len(stas)))-0.5, max(range(len(stas)))+0.5))



    #if upo >= 30:
    #    upo=30.
    #if zrat >=30:
    #    zrat = 30.
    #if wpo > 100:
    #    wpo=100.
#     if rrat >=100:
#         rrat = 100.
#     if zrat >= 100:
#         continue

#     #if zrat > 20:
#     # if zcorr >= 0.8:
#     #     print(line[0] + ' ' + str(zrat) + ' '  + str(upo))

#     if (zcorr >= 0.8) and (upo <40):
#         zrats.append(zrat)
#         upos.append(upo)
#         zcorrs.append(zcorr)
#     else:
#         zratsb.append(zrat)
#         uposb.append(upo)
#         zcorrsb.append(zcorr)
#     if rcorr >= 0.8:
#         rrats.append(rrat)
#         rcorrs.append(rcorr)
#     #    wpos.append(wpo)

# # pzv = np.polyfit(upos, zrats,1)
# print(pzv)
# pz = np.poly1d(pzv)
# #prv= np.polyfit(wpos, rrats,1)
# #pr = np.poly1d(prv)
# #plt.subplot(1,2,1)
# ec=plt.scatter(upos, zrats, c=zcorrs, marker='o', vmin=0, vmax=1.)
# plt.plot(np.linspace(0,90,100), pz(np.linspace(0,90,100)), color='r', label= 'Slope=' + str(round(pzv[0],3)) )
# plt.scatter(uposb, zratsb, c=zcorrsb, marker='s', vmin=0, vmax=1., alpha=0.5)
# plt.plot(np.linspace(0,90,100),np.linspace(0,90,100), label= 'Slope=1', color='k')
# plt.xlabel('Vertical Ratio Crust 1.0 (nm/s/Pa)')
# plt.ylabel('Vertical Ratio Hunga Tonga (nm/s/Pa)')
# plt.xlim((0.,81.))
# plt.ylim(0.,81.)
# plt.legend(loc='upper left')
# cbar = plt.colorbar(ec)
# cbar.set_label('Correlation Vertial to Hilbert Transform of Pressure')
# # plt.subplot(1,2,2)
# # plt.scatter(wpos, rrats, c=rcorrs)
# # plt.plot(np.linspace(0,30,100), pr(np.linspace(0,30,100)), color='r')
# # plt.xlabel('Radial Ratio Crust 1.0 (nm/s/Pa)')
# # plt.ylabel('Radial Ratio Hunga Tonga (nm/s/Pa)')
# # plt.xlim((0.,100.))
# # plt.ylim(0.,100.)
# # #plt.tight_layout()
plt.savefig('Figure102GLAD.png', format='PNG', dpi=400)
plt.savefig('Figure102GLAD.pdf', format='PDF', dpi=400)
# f.close()
