#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime, Stream, read
from obspy.geodetics.base import gps2dist_azimuth
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
from obspy.signal.cross_correlation import correlate, xcorr_max

import matplotlib as mpl
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font', size=18)




debug = True

# Tonga latitude and longitude (from USGS) #####################################


Tonga_Lat = -20.546
Tonga_Lon = -175.390


# Eruption Time
stime = UTCDateTime('2022-01-15T04:14:00')
# TUC
etime = UTCDateTime('2022-01-16T00:00:00')


# velocity of Lamb pulse (m/s)
Lamb_V = 315.

# Seismic station you want to look at
net = 'II'
sta = 'PFO'

# delay time from predicted arrival of lamb Pulse
# size of window to do correlation
WD = 0.75
win = 0.5

# plotting window
Xmin = 0
Xmax = 600


fmin = 1./90.
fmax = 1./40.



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

f = open('Theresults_Jordan','w')

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

    print(st)
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

    print(sta, 'Distance to Hunga is ', dist/1000 )

    # Calculate arrival of Lamb Pulse from eruption
    TD = dist/Lamb_V


    #st2 = st.copy()
    st.filter('bandpass',freqmin=fmin, freqmax=fmax, zerophase=True)
    st.trim(stime+TD, stime+TD+4*60*60 )
    #st2.trim(stime+TD, stime+TD+4*60*60 )
    
    nstime = st[0].stats.starttime + np.argmax(st.select(channel='LD*')[0].data)

    #st.filter('bandpass', freqmin=fmin, freqmax=fmax, zerophase=True)
    st.trim(nstime -10*60, nstime-10*60 + (win*60*60))
    #fig = plt.figure(1)
    #for idx, tr in enumerate(st):
    #    plt.subplot(4,1,idx+1)
    #    plt.plot(tr.times(), tr.data, label=tr.id)
    #    plt.legend()
    #plt.show()


    # filter and trim data
    
    # 4 hours window starting 10 minutes before predicted
    #st.trim(stime+TD - 10*60, stime+TD+4*60*60 -10*60)
    #nstime = st[0].stats.starttime + np.argmax(st.select(channel='LD*')[0].data)
    #st.trim(nstime -10*60, nstime-10*60 + (win*60*60))
    #st.trim(stime+TD+WD*60*60, stime+TD+(WD*60*60)+(win*60*60))


    # Now we need to rotate everything to radial

    # first get station coordinates



    # Now we rotate the seismometer first from arbitary orientation to N and # -*- coding: utf-8 -*-
    st.rotate('->ZNE', inventory=inv)
    bestcor, bestangle = -1.,0.
    st2 = st.copy()
    for angle in range(-10,11):
        st2 = st.copy()
        nang = (BAZ+angle) % 360.
        st2.rotate('NE->RT',back_azimuth=nang)

        try:
            # Grab the Traces we want
            TR_P = st2.select(channel='LD*')[0]
            TR_Z = st2.select(channel='LHZ')[0]
            TR_R = st2.select(channel='LHR')[0]
        except:
            continue
        # Now we do our Hilber Transforms
        # Take the Hilbert Transforms that we want
        # Keep in mind these will no longer be trave objects

        TR_RHT = np.imag(hilbert(TR_R.data))
        TR_PHT = np.imag(hilbert(TR_P.data))

        # Now we do our correlations

        cc_PR =  correlate(TR_R.data, TR_P.data,0)
        shift_PR, value_PR = xcorr_max(cc_PR)
        value_PR = np.round(value_PR, decimals=2)
        print( sta, 'Radial to Pressure Correlation is: ', value_PR)
        rat_PR = np.std(TR_R.data)/np.std(TR_P.data)
        if value_PR > bestcor:
            bestcor = value_PR
            bestangle = angle


        cc_ZPHT=  correlate(TR_Z.data, TR_PHT,0)
        shift_ZPHT, value_ZPHT = xcorr_max(cc_ZPHT)
        value_ZPHT = np.round(value_ZPHT, decimals=2)
        print( sta, 'Vertical to HT Pressure Correlation is: ', value_ZPHT)
        rat_ZPHT= np.std(TR_Z.data)/np.std(TR_PHT)




        cc_Ray =  correlate(TR_Z.data,TR_RHT,0)
        shift_Ray, value_Ray = xcorr_max(cc_Ray)
        value_Ray = np.round(value_Ray, decimals=2)
        print(sta, 'HT Radial to Vertical Correlation is ', value_Ray)
        rat_Ray= np.std(TR_Z.data)/np.std(TR_RHT)



    # sta, lat, lon, RP_corr, VHTpress_corr, HTV_corr
    f.write(sta +', ' + str(s_lat) + ', ' + str(s_lon) + ', ' + str(value_PR) +', '  
        + str(rat_PR) + ', ' + str(value_ZPHT) + ', ' + str(rat_ZPHT) + ', '+str(value_Ray) +', ' 
        + str(rat_Ray) + ', ' + str(bestcor) + ', ' + str(bestangle) + '\n')


f.close()
