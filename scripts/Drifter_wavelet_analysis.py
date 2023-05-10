#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 12:36:35 2023

@author: Ashleigh Womack

Wavelet Analysis of the Drifter
    - This yields information about the time series and frequencies together. By
      decomposing a time series into a time-frequency space, it allows for the 
      determination of the dominant modes of variability, and how these modes vary
      in time. 
      
First run script: wavelet_functions.py 
"""

#%%
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from waveletFunctions import wavelet, wave_signif
from matplotlib.dates import DayLocator,DateFormatter
        
#%% Read in the buoy drifter data 
BDIR = '../data/'
theBuoy = 'ISVP1'                                     
drifter = pd.read_csv(BDIR+theBuoy+'.csv', index_col='time', parse_dates=True)

#%% If NOT appling the high pass filter
#   - Use these variables for the wavelet analysis
#   First index is nan - therefore skip it. 

x = drifter['u (m/s)'][1:]                                               # zonal velocity
y = drifter['v (m/s)'][1:]                                               # meridional velocity
var_x = np.std(x, ddof=1) ** 2                                           # variance of x (standard deiviation)
var_y = np.std(y, ddof=1) ** 2                                           # variance of y
t = drifter.index[1:]
dt = 1                                                                   # sampling interval of drifter (in hours)

#%% Application of a Butterworth High-Pass Filter 
#   This is a 12th order filter used to remove all frequencies below the cutoff
#   frequency of 0.04 Hz which is approximately the daily frequency, allowing 
#   the sub-daily frequencies to become more significant. 

dt = 1                                                                   # sampling intervals in hours
Fs = 1/dt                                                                # sampling frequency (1/dt)
b, a = signal.butter(12, 0.04, 'high', Fs)                               # 12th order with 0.04 h-1 cutoff frequency            
w, h = signal.freqs(b, a)                                                # frequency response of a digital filter
x_filt = signal.filtfilt(b, a, drifter['u (m/s)'][1:], padtype=None)     # filter data series of zonal component (or speed)
y_filt = signal.filtfilt(b ,a, drifter['v (m/s)'][1:], padtype=None)     # filter data series of meridional component

x = x_filt                           
y = y_filt                                                           
var_x = np.std(x, ddof=1) ** 2                                           # varience of x
var_y = np.std(y, ddof=1) ** 2                                           # varience of y
t = drifter.index[1:]

#%% Load the variables for the Wavelet transform
NFFT = np.size(x)                                                        # FFT length
pad = 1                                                                  # pad the time series with zeroes (recommended)
dj = 0.25                                                                # this will do 4 sub-octaves per octavet
s0 = 2 * dt                                                              # the smallest scale of the wavelet
j1 = 7 / dj                                                              # does 7 powers-of-two with dj sub-octaves each
lag1 = 0.72                                                              # lag-1 autocorrelation for red noise background
print("lag1 = ", lag1)  
mother = 'MORLET'                                                        # choice of mother wavelet function

# Complete the wavelet transform:
wave, period, scale, coi = wavelet(x, dt, pad, dj, s0, j1, mother)       # change to y
power = (np.abs(wave)) ** 2                                              # compute wavelet power spectrum
wt_spectrum = (np.nansum(power, axis=1) /NFFT)                           # time-average over all times

# Significance levels:
signif = wave_signif(([var_x]), dt=dt, sigtest=0, scale=scale,           # change to var_y
    lag1=lag1, mother=mother)
sig95 = signif[:, np.newaxis].dot(np.ones(NFFT)[np.newaxis, :])          # expand signif --> (J+1)x(N) array
sig95 = power / sig95                                                    # where ratio > 1, power is significant

# Wavelet spectrum & significance levels:
dof = NFFT - scale                                                       # the -scale corrects for padding at edges
signif_wt = wave_signif(var_x, dt=dt, scale=scale, sigtest=1,            # change to var_y
    lag1=lag1, dof=dof, mother=mother)

# Repeat for the meridional component of the drifter - where it has been stated "change to y or var_y"

#%% Plotting Wavelet Power Spectrum and Wavelet Spectrum - Plots zonal component (change to meridional)
mpl.rcParams['font.size'] = 30
fig = plt.figure(figsize=(30, 23))
gs = GridSpec(3, 4, wspace=0.025, hspace=0.05)
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.95, wspace=0, hspace=0)

# Wave Power Spectrum 
plt1 = plt.subplot(gs[1, 0:3])
CS = plt.contourf(t, period, power, cmap=plt.cm.BuGn) # , vmin=0, vmax=2.4)                           
plt.xlabel('Date')
plt.ylabel('Period (hours)')
plt.title('Wavelet Power Spectrum', fontsize= 30)                        
plt.contour(t, period, sig95, [-99, 1], colors='k')                      # significance contour, levels at -99 (fake) and 1 (95% signif)  
# cone-of-influence (coi), anything "below" is dubious - edge effects           
plt.plot(t, coi[1:], 'r', linewidth=2.5)                                 # remove the 1 if sizes don't match                                                  
plt1.set_yscale('log', base=2, subs=None)
plt1.xaxis.set_major_locator(DayLocator(interval=10))
plt1.xaxis.set_major_formatter(DateFormatter('%d-%m'))
plt.ylim(2.5,256)
ax = plt.gca().yaxis
ax.set_major_formatter(ticker.ScalarFormatter())
plt1.ticklabel_format(axis='y', style='plain')
plt1.invert_yaxis()
plt.colorbar(CS)
# if you need to scale the colorbar
#from matplotlib.cm import ScalarMappable
#cb = plt.colorbar(ScalarMappable(norm=CS.norm,cmap=plt.cm.BuGn),
#                  ticks=[0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4],
#                  boundaries=np.arange(0, 2.5, 0.3))

# Wavelet Spectrum
mpl.rcParams['font.size'] = 30
plt2 = plt.subplot(gs[1, -1])
plt.plot(wt_spectrum, period, linewidth=3, color='mediumseagreen')      # wavelet power spectrum
plt.plot(signif_wt, period, '--', linewidth=3, color='black')           # wavelet significance level 
plt.xlabel('Power (variance)') 
plt.title('Wavelet Spectrum', fontsize= 30)                            
plt.xlim([0, 1.25 * np.max(wt_spectrum)])
plt2.set_yscale('log', base=2, subs=None)
plt.ylim(2.5,256)
plt.xlim(0,1.5)
ax = plt.gca().yaxis
ax.set_major_formatter(ticker.ScalarFormatter())
ax.set_major_formatter(plt.NullFormatter())
plt2.invert_yaxis()

# end of code
