#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:18:57 2021

@author: Ashleigh Womack

Wavelet Analysis of the Drifter
    - This yields information about the time series and frequencies together. By
      decomposing a time series into a time-frequency space, it allows for the 
      determination of the dominant modes of variability, and how these modes vary
      in time. 

First - compute the PSD (or filtered spectrum) of the drifter.

NB: it requires wavelet_functions.py 
"""

#%%
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from waveletFunctions import wavelet,wave_signif
from matplotlib.dates import DayLocator,DateFormatter
#%% Use data from the power spectral analysis - read in either the unfiltered data
## or the filtered data

#drifter = drifter_1
                                   # varience (standard deiviation)            
#%% Read in the buoy drifter data 
BDIR = '../data/'
theBuoy = 'ISVP1'                                     
drifter = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
x = drifter['u (m/s)']                                        # zonal velocity
y = drifter['v (m/s)']                                         # meridional velocity
x = drifter['U (m/s)']                                         # Speed
var_x = np.std(x, ddof=1) ** 2                                 # variance (standard deiviation)
var_y = np.std(y, ddof=1) ** 2                                 # variance of y
t = drifter.index
#%% Load the variables for the Wavelet transform
dt = 1                                                              # sampling intervals in hours
NFFT = np.size(x)                                                   # FFT length
pad = 1                                                             # pad the time series with zeroes (recommended)
dj = 0.25                                                           # this will do 4 sub-octaves per octavet
s0 = 2 * dt                                                         # the smallest scale of the wavelet
j1 = 7 / dj                                                         # does 7 powers-of-two with dj sub-octaves each
lag1 = 0.72                                                         # lag-1 autocorrelation for red noise background
print("lag1 = ", lag1)  
mother = 'MORLET'                                                   # choice of mother wavelet function

# Complete the wavelet transform:
wave, period, scale, coi = wavelet(x, dt, pad, dj, s0, j1, mother)  # change to y
power = (np.abs(wave)) ** 2                                         # compute wavelet power spectrum
wt_spectrum = (np.nansum(power, axis=1) /NFFT)                      # time-average over all times

# Significance levels:
signif = wave_signif(([var_x]), dt=dt, sigtest=0, scale=scale,      # change to var_y
    lag1=lag1, mother=mother)
sig95 = signif[:, np.newaxis].dot(np.ones(NFFT)[np.newaxis, :])     # expand signif --> (J+1)x(N) array
sig95 = power / sig95                                               # where ratio > 1, power is significant

# Wavelet spectrum & significance levels:
dof = NFFT - scale                                                  # the -scale corrects for padding at edges
signif_wt = wave_signif(var_x, dt=dt, scale=scale, sigtest=1,       # change to var_y
    lag1=lag1, dof=dof, mother=mother)

# Repeat for the meridional component of the drifter - where it has been stated "change to y or var_y"

#%% Plotting Wavelet Power Spectrum and Wavelet Spectrum
mpl.rcParams['font.size'] = 30
fig = plt.figure(figsize=(30, 23))
gs = GridSpec(3, 4, wspace=0.025, hspace=0.05)
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.95, wspace=0, hspace=0)

# Wave Power Spectrum
plt1 = plt.subplot(gs[1, 0:3])
CS = plt.contourf(t, period, power, cmap=plt.cm.BuGn, vmin=0, vmax=3.2)                            
plt.xlabel('Date')
plt.ylabel('Period (hours)')
plt.title('Wavelet Power Spectrum', fontsize= 30)                   # change to (V component)
plt.contour(t, period, sig95, [-99, 1], colors='k')                 # significance contour, levels at -99 (fake) and 1 (95% signif)            
plt.plot(t, coi, 'r', linewidth=2.5)                            # cone-of-influence, anything "below" is dubious                            
plt1.set_yscale('log', base=2, subs=None)
plt1.xaxis.set_major_locator(DayLocator(interval=10))
plt1.xaxis.set_major_formatter(DateFormatter('%d-%m'))
plt.ylim([np.min(period[0]), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(ticker.ScalarFormatter())
plt1.ticklabel_format(axis='y', style='plain')
plt1.invert_yaxis()
plt.colorbar()
# if you need to specify colorbar limits
#from matplotlib.cm import ScalarMappable
#cb = plt.colorbar(ScalarMappable(norm=CS.norm,cmap=plt.cm.BuGn),
#                  ticks=[0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2],
#                  boundaries=np.arange(0, 3.6, 0.4))

# Wavelet Spectrum
mpl.rcParams['font.size'] = 30
plt2 = plt.subplot(gs[1, -1])
plt.plot(wt_spectrum, period, linewidth=3, color='mediumseagreen')   # wavelet power spectrum
plt.plot(signif_wt, period, '--', linewidth=3, color='black')        # wavelet significance level 
plt.xlabel('Power (variance)')
plt.title('Wavelet Spectrum', fontsize= 30)                          # change to (V component)
plt.xlim([0, 1.25 * np.max(wt_spectrum)])
plt2.set_yscale('log', base=2, subs=None)
plt.ylim([np.min(period[0]), np.max(period)])
plt.xlim(0,1.5)
ax = plt.gca().yaxis
ax.set_major_formatter(ticker.ScalarFormatter())
ax.set_major_formatter(plt.NullFormatter())
plt2.invert_yaxis()

fig.tight_layout()

# end of code
