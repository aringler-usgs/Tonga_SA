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

# Tonga latitude and longitude and event time (from USGS) #####################################

Tonga_Lat = -20.546
Tonga_Lon = -175.390

Event_time = UTCDateTime('2022-015T04:14:45')


# stime
stime = UTCDateTime('2022-01-15T00:00:00')
# TUC
etime = UTCDateTime('2022-01-16T00:00:00')


# velocity of Lamb pulse (m/s)
Lamb_V = 315.

# Seismic station you want to look at
net = 'IU'
sta = 'ANMO'

# delay time from predicted arrival of lamb Pulse
# size of window to do correlation

win = 0.5

# plotting window
Xmin = 0
Xmax = 60*60*win

Filt_min = 1/90.
Filt_max = 1/40.

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


# Remove Response of Seismic data
st.select(channel='LHZ')[0].remove_response(inventory=inv)
st.select(channel='LH1')[0].remove_response(inventory=inv)
st.select(channel='LH2')[0].remove_response(inventory=inv)
# Remove Poynomial Response of Pressure Data
remove_polynomial(st.select(channel='LDO'),inv)

st.filter('bandpass', freqmin=Filt_min, freqmax=Filt_max)

# Now we do the stuff to get the window

coors = inv.get_coordinates(net+'.'+sta+'.00.LHZ')
s_lat = coors['latitude']
s_lon = coors['longitude']

# Now we get dist and the BAZ to the Event
dist,BAZ, BAZ_INV = gps2dist_azimuth(s_lat,s_lon,Tonga_Lat,Tonga_Lon)
print(sta, 'Distance to Hunga is ', dist/1000 )

# Calculate arrival of Lamb Pulse from eruption
TD = dist/Lamb_V

# now we get the peak amplitude
st2 = st.copy()
st2.trim(Event_time+TD, Event_time+TD+4*60*60 )
nstime = st2[0].stats.starttime + np.argmax(st2.select(channel='LD*')[0].data)

# trim the original trace around the peak acoustic amplitude

st.trim(nstime -10*60, nstime-10*60 + (win*60*60))

# Now we need to rotate everything to radial
# Now we rotate the seismometer first from arbitary orientation to N and # -*- coding: utf-8 -*-
st.rotate('->ZNE', inventory=inv)
st.rotate('NE->RT',back_azimuth=BAZ)

# Grab the Traces we want
TR_P = st.select(channel='LDO')[0]
TR_Z = st.select(channel='LHZ')[0]
TR_R = st.select(channel='LHR')[0]

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


cc_ZPHT=  correlate(TR_Z.data, TR_PHT,0)
shift_ZPHT, value_ZPHT = xcorr_max(cc_ZPHT)
value_ZPHT = np.round(value_ZPHT, decimals=2)
print( sta, 'Vertical to HT Pressure Correlation is: ', value_ZPHT)


cc_Ray =  correlate(TR_Z.data,TR_RHT,0)
shift_Ray, value_Ray = xcorr_max(cc_Ray)
value_Ray = np.round(value_Ray, decimals=2)
print(sta, 'HT Radial to Vertical Correlation is ', value_Ray)

# Now we get our ratios

RMS_Sv, RMS_Sr, RMS_A = np.std(TR_Z.data*1E9), np.std(TR_R.data*1E9), np.std(TR_P.data)
Sv_A, Sr_A = RMS_Sv/RMS_A, RMS_Sr/RMS_A

print('Vertical Seismic/Acoustic (nm/s/Pa) is:', Sv_A)
print('Radial Seismic/Acoustic (nm/s/Pa) is:', Sr_A)



# Creating figure and moving subplots to have no space between them
fig = plt.figure(1, figsize=(20,16))
#plt.subplots_adjust(hspace=0.001)
#plt.subplots_adjust(vspace=0.1)

# Radial Vs. Pressure
ax1 = plt.subplot(311)
ax1.plot(TR_R.times(),TR_R.data*1000000.,linewidth=3,c='C1',alpha=1.0, label= sta + ' Radial: ' + str(value_PR))
ax1.set_ylabel('Ground\nVelocity ($\mu$m/s)', fontsize=24, c='C1')
ax1.tick_params(labelsize = 20, colors='C1')
plt.xlim(Xmin, Xmax)
#plt.xlabel('Time (s)', fontsize = 18)
plt.legend(loc=2)
plt.text(-50, 0.08, '(a)', fontsize=24 )


ax2 = ax1.twinx()
ax2.plot(TR_P.times(),TR_P.data,linewidth=2,c='C2',alpha=1.0)
ax2.set_ylabel('Pressure (Pa)', fontsize=24, c='C2', )
ax2.tick_params(axis='y', colors='C2', labelsize=20)
ax2.set_xlim(Xmin, Xmax)



# Vertical vs HT Pressure
ax3 = plt.subplot(312)
ax3.plot(TR_Z.times(),TR_Z.data*1000000.,linewidth=3,c='k',alpha=1.0, label= sta +' Vertical: ' + str(value_ZPHT))
ax3.set_ylabel('Ground\nVelocity ($\mu$m/s)', fontsize=24)
plt.tick_params(labelsize = 20)
#plt.text(Xmin-0.05,1.2,'(b)', fontsize=24)
plt.xlim(Xmin, Xmax)
plt.legend(loc=2)
plt.text(-50, 0.12, '(b)', fontsize=24 )

ax4 = ax3.twinx()
ax4.plot(TR_P.times(),TR_PHT,linewidth=2,c='C2',alpha=1.0, linestyle='dashed', label='Pressure HT)')
ax4.set_ylabel('HT Pressure (Pa)', fontsize=24, c='C2')
ax4.tick_params(axis='y', colors='C2', labelsize=20)
ax4.set_xlim(Xmin, Xmax)

# Vertical vs HT Radial
ax5 = plt.subplot(313)
ax5.plot(TR_Z.times(),TR_Z.data*1000000.,linewidth=3,c='k',alpha=1.0, label= sta + ' Vertical: ' + str(value_Ray))
ax5.plot(TR_R.times(),TR_RHT*1000000.,linewidth=3,c='C1',linestyle='dashed')
ax5.set_ylabel('Ground\nVelocity ($\mu$m/s)', fontsize=24)
plt.tick_params(labelsize = 20)
plt.xlabel('Time (s)', fontsize = 24)
#plt.text(Xmin-0.05,1.2,'(c)', fontsize=24)
plt.xlim(Xmin, Xmax)
#plt.xlabel('Time (s)', fontsize = 18)
plt.legend(loc=2)
plt.text(-50, 0.12, '(c)',fontsize=24 )

########################################################################

plt.show()
