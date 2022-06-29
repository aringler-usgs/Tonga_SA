#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime, Stream, read
from obspy.geodetics.base import gps2dist_azimuth
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
from obspy.signal.cross_correlation import correlate, xcorr_max
from scipy import signal
import matplotlib as mpl
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font', size=18)




debug = False

# Tonga latitude and longitude (from USGS) #####################################


Tonga_Lat = -20.546
Tonga_Lon = -175.390


# Eruption Time
stime = UTCDateTime('2022-01-15T04:14:00')
# TUC
etime = UTCDateTime('2022-01-16T00:00:00')


# velocity of Lamb pulse (m/s)
Lamb_V = 315.


# delay time from predicted arrival of lamb Pulse
# size of window to do correlation
WD = 0.75
win = 1.0

# plotting window
Xmin = 0
Xmax = 600


fmin = 1./90.
fmax = 1./40.
permin = 1./fmax
permax = 1./fmin


client = Client("IRIS")

stas = []
inv = client.get_stations(network='IU,IC,CU,II', station='*', starttime=stime-60*60*win,
                          endtime=etime+60*60*win, level="response",
                          location='00,31', channel='LDO,LH*,LDI')


for net in inv:
    for sta in net:
        stas.append(sta.code)


# Adam's function to remove polynomial response#########################

def remove_polynomial(st, inv):
    # Simple function to remove the response of polynomial responses
    for tr in st:
        # Grab polynomial response
        resp = inv.get_response(tr.id, tr.stats.starttime).instrument_polynomial.coefficients
        # Go to floats to avoid type issues.
        tr.data = np.require(tr.data, dtype=np.float64)
        # Make a data vector of zeros
        data_temp = np.zeros_like(tr.data)
        # For each coefficient add the data to data_temp
        for idx, coef in enumerate(resp):
            data_temp += coef*tr.data**idx
        tr.data = data_temp

    return



fig =plt.figure(1,figsize=(12,12))
plotme = True
for idx, sta in enumerate(stas):


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


    if debug:
        print(st)

    st.detrend('linear')
    st.merge(fill_value=0)


    # Remove Response of Seismic data

    if st[0].stats.network == 'II':
        st.remove_response(inventory=inv)
    else:
        st.select(channel='LH*').remove_response(inventory=inv)
        remove_polynomial(st.select(channel='LD*'),inv)




    coors = inv.get_coordinates(st[0].id)
    s_lat = coors['latitude']
    s_lon = coors['longitude']

    # Now we get dist and the BAZ to the Event
    dist,BAZ, BAZ_INV = gps2dist_azimuth(s_lat,s_lon,Tonga_Lat,Tonga_Lon)

    #print(sta, 'Distance to Hunga is ', dist/1000 )

    # Calculate arrival of Lamb Pulse from eruption
    TD = dist/Lamb_V

    
    #nstime = st[0].stats.starttime + np.argmax(st.select(channel='LD*')[0].data)
    #print('Here is the starttime:' + str(st[0].stats.starttime))
    #print('Here is the new starttime' + str(nstime))


    # filter and trim data
    #st.filter('bandpass', freqmin=fmin, freqmax=fmax, zerophase=True)
    
    #st.trim(stime+TD+WD*60*60, stime+TD+(WD*60*60)+(win*60*60))
    # 4 hours window starting 10 minutes before predicted
    st2 = st.copy()
    st2.filter('bandpass',freqmin=fmin, freqmax=fmax, zerophase=True)

    st2.trim(stime+TD, stime+TD+4*60*60 )
    st.trim(stime +TD, stime+TD+4*60*60)
    nstime = st[0].stats.starttime + np.argmax(st2.select(channel='LD*')[0].data)
    st.trim(nstime -10*60, nstime-10*60 + (win*60*60))
    speed = dist / (nstime - stime)
    print(st[0].stats.station + ' ' + str(nstime) + ' ' + str(dist) + ' ' + str(speed))
    #############

    # Now we need to rotate everything to radial

    # first get station coordinates



    # Now we rotate the seismometer first from arbitary orientation to N and # -*- coding: utf-8 -*-
    st.rotate('->ZNE', inventory=inv)
    st.rotate('NE->RT',back_azimuth=BAZ)

    try:
        # Grab the Traces we want
        TR_P = st.select(channel='LD*')[0]
        TR_Z = st.select(channel='LHZ')[0]
        TR_R = st.select(channel='LHR')[0]
    except:
        continue



    
    plt.subplot(2,1,2)
    if plotme:
        plt.axvspan(permin,permax, color='C2')

    f, pz = signal.coherence(TR_P.data, TR_Z.data, 1, nperseg=1024)
    f, pz = f[1:], pz[1:]
    per =1./f
    valz = np.mean(pz[(per >= permin) & (per <= permax)])



    plt.semilogx(1./f, pz, alpha=0.1, color='C0')
    plt.subplot(2,1,1)
    if plotme:
        plt.axvspan(permin,permax, color='C2')
        plotme = False
    f, pr = signal.coherence(TR_P.data, TR_R.data, 1, nperseg=1024)
    f, pr = f[1:], pr[1:]
    per =1./f
    valr = np.mean(pr[(per >= permin) & (per <= permax)])
    plt.semilogx(1./f, pr, alpha=0.1, color='C0')
    try:
        if valr >= 0.7:
            if 'prs' in vars():

                prs = np.vstack((prs, pr))
            
            else:
                prs = pr
    except:
        continue
    try:
        if valz >= 0.7:
            if 'pzs' in vars():

                pzs = np.vstack((pzs, pz))
            
            else:
                pzs = pz
    except:
        continue

pz= np.mean(pzs, axis=0)
pr = np.mean(prs, axis=0)
plt.subplot(2,1,2)
plt.semilogx(1./f, pz, color='C1', label='Vertical Mean')
plt.xlim((2,1000))
plt.text(.9, 1, '(b)' )
plt.ylabel('Coherence ($\gamma^2$)')
plt.legend(loc=2)
plt.subplot(2,1,1)
plt.text(.9, 1, '(a)' )
plt.semilogx(1./f, pr, color='C1', label='Radial Mean')
plt.xlim((2,1000))
plt.xlabel('Period (s)')
plt.ylabel('Coherence ($\gamma^2$)')
plt.legend(loc=2)

plt.savefig('Figure2.png', format='PNG', dpi=400)
plt.savefig('Figure2.pdf', format='PDF', dpi=400)



