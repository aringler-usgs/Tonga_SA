#!/usr/bin/env python

from obspy.clients.fdsn import Client
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import obspy
from obspy import read_inventory, UTCDateTime, read, Stream
from obspy.geodetics.base import gps2dist_azimuth

#Global Variables
# Distance from ANMO to Tonga (km)

debug = True

Event_time = UTCDateTime('2022-015T04:14:45')
# Rayleigh wave velocity and Acoustic velocity (km/s)
# Pressure taken from Junghyun Park, SMU estimate

# Tonga latitude and longitude (from USGS) #####################################


Tonga_Lat = -20.546
Tonga_Lon = -175.390

# velocity of Lamb pulse (m/s)
Lamb_V = 315.


# Eruption Time
stime = UTCDateTime('2022-01-15T00:00:00')
# TUC
etime = UTCDateTime('2022-01-16T00:00:00')


# Seismic station you want to look at
net = 'IU'
sta = 'GNI'

# delay time from predicted arrival of lamb Pulse
# size of window to do correlation

win = 0.5

# plotting window
Xmin = 0
Xmax = 1800

I_Filt_min = 1/100.
I_Filt_max = 1/10.

S_Filt_min = 1/90.
S_Filt_max = 1/40.



client = Client("IRIS")

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




# Get the data from IRIS pad time to account for filter ring

inv = client.get_stations(network='IU,IC,CU,II', station=sta, starttime=stime-60*60*win,
                          endtime=etime+60*60*win, level="response",
                          location='00,31', channel='LDO,LH*')


st = client.get_waveforms(network='IU,IC,CU,II', station=sta, location='00,31',
                          channel='LDO,LH*', starttime=stime-60*60*win,
                          endtime=etime + win*60*60)


if debug:
    print(st)

st.detrend('linear')
st.merge(fill_value=0)
#st.remove_response(inventory=inv)

# Remove Response of Seismic data
st.select(channel='LHZ')[0].remove_response(inventory=inv)
st.select(channel='LH1')[0].remove_response(inventory=inv)
st.select(channel='LH2')[0].remove_response(inventory=inv)

#if st.select(channel='LDI')[0].stats.network == 'II':
#    st.select(channel='LDI')[0].remove_response(inventory=inv)
# Remove Poynomial Response of Pressure Data (do for ANMO)
remove_polynomial(st.select(channel='LDO'),inv)

# Filter the data

st.select(channel='LHZ')[0].filter('bandpass',freqmin=S_Filt_min, freqmax = S_Filt_max ,corners=4)
# Convert seismic data to micrometers per second
st.select(channel='LHZ')[0].data *= 1000000


st.select(channel='LDO')[0].filter('bandpass',freqmin=I_Filt_min, freqmax = I_Filt_max ,corners=4)

# Trim to remove filter Ring
st.trim(stime,etime)


# Now we get the windows
# 1) 4 Hour search window

coors = inv.get_coordinates(net+'.'+sta+'.00.LHZ')
s_lat = coors['latitude']
s_lon = coors['longitude']

# Now we get dist and the BAZ to the Event
dist,BAZ, BAZ_INV = gps2dist_azimuth(s_lat,s_lon,Tonga_Lat,Tonga_Lon)

if debug:
    print(sta, 'Distance to Hunga is ', dist/1000 )

# Calculate arrival of Lamb Pulse from eruption
TD = dist/Lamb_V

st2 = st.copy()
st2.trim(Event_time+TD, Event_time+TD+4*60*60 )
nstime = st2[0].stats.starttime + np.argmax(st2.select(channel='LD*')[0].data)



# get decimal hour of eruption

Eruption_DH = Event_time.hour + ((Event_time.second/60.) + Event_time.minute)/60.
Eruption_tvec = [Eruption_DH,Eruption_DH]

# Lamb window
L_S = Eruption_DH + (TD/3600)

# Event window
E_S = L_S + np.argmax(st2.select(channel='LD*')[0].data)/3600. - 1./6.

######## Make the Figure ###########################################
fig = plt.figure(1, figsize=(12,16))

ax1 = plt.subplot(111)
ax1.plot(st.select(channel = 'LHZ')[0].times()/3600., st.select(channel = 'LHZ')[0].data,linewidth=1.5,c='k',alpha=0.7)
ax1.axvspan(Eruption_DH-0.03, Eruption_DH+0.03, color='k')
ax1.text(Eruption_DH + 0.1, 1.0*st.select(channel = 'LHZ')[0].max(), 'Hunga Eruption', fontsize=15)
# 4-Hours after lamp arrival
ax1.axvspan(L_S, L_S +4, color='C1', alpha = 0.3)
ax1.axvspan(E_S, E_S + win, color='C0', alpha = 0.3)


ax1.set_ylabel('Ground Velocity ($\mu$m/s)', fontsize=18, c='k')
ax1.tick_params(axis='y', colors='k', labelsize=15)
ax1.tick_params(axis='x', labelsize=15)

#ax1.set_ylim(-1.5, 1.8)
ax1.set_xlim(2.0,24.0)
ax1.set_xlabel('Hour of 15th January', fontsize=18)



ax2 = ax1.twinx()
ax2.plot(st.select(channel = 'LDO')[0].times()/3600., st.select(channel = 'LDO')[0].data,linewidth=1.0,c='C2',alpha=0.7)
ax2.set_ylabel('Pressure (Pa)', fontsize=18, c='C2' )
ax2.tick_params(axis='y', colors='C2', labelsize=15)

plt.show()
