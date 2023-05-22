#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 13:15:18 2023

@author: Ashleigh Womack (UCT)

Relative dispersion and proxy for total deformation of the 2019
Winter Cruise buoys.
First read in Rel_disp_functions.py
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from Rel_disp_functions import separation,rel_disp_0,disp_0,std_disp,rel_disp
from Rel_disp_functions import delta_r_sq,disp_rate,std_disp_rate,rel_disp_var
from matplotlib.dates import DayLocator,DateFormatter

# directory of the processed drifter data
BDIR_P = '../data/'

#%% Read in the buoy drifter data 
theBuoy = 'ISVP1'                                     
drifter_1 = pd.read_csv(BDIR_P+theBuoy+'.csv',index_col='time',parse_dates=True)

theBuoy = 'ISVP2'                                     
drifter_2 = pd.read_csv(BDIR_P+theBuoy+'.csv',index_col='time',parse_dates=True)

theBuoy = 'ISVP3'                                     
drifter_3 = pd.read_csv(BDIR_P+theBuoy+'.csv',index_col='time',parse_dates=True)


#%% Manually resample the data for the cluster dates
#   for every hour or any regular time interval. 
#   - Note this will change the number of indicies in a day. 
time = pd.date_range(start="2019-07-28 05:00:00", 
                         end="2019-08-25 03:00:00", freq="1H")

# find the nearest longitude from buoy to corresponding sea ice contour
def find_ind(array, value):
    idx = np.nanargmin((np.abs(array - value)))
    return idx

dates_1 = np.zeros(len(time))
dates_2 = np.zeros(len(time))
dates_3 = np.zeros(len(time))
i=0
for i in range(len(time)):
    dates_1[i] = find_ind(drifter_1.index,time[i])
    dates_2[i] = find_ind(drifter_2.index,time[i])
    dates_3[i] = find_ind(drifter_3.index,time[i])
    i=i+1

drifter_1 = drifter_1.iloc[dates_1] 
drifter_2 = drifter_2.iloc[dates_2]  
drifter_3 = drifter_3.iloc[dates_3]

#################################################################
#%% Calculation of the deformation of the cluster 
#   Using Rampal and others (2009) method - (small and big lengths)
#   We initally compute this from L0
#################################################################

# First: calculate the separation from L0 between each buoy pair
d1d2 = separation(drifter_1, drifter_2) 
d1d3 = separation(drifter_1, drifter_3)
d2d3 = separation(drifter_2, drifter_3) 

dr2_list_0 = [d1d2['dr2_0'], d1d3['dr2_0'], d1d3['dr2_0']]
freq = 24                                                    # how many indices in a day 
nob  = 3                                                     # number of buoy pairs in each length group

# Second: calculate the mean squared change in separation
RD_0 = rel_disp_0(dr2_list_0, time, nob=nob, freq=freq)      # 24 indices * 1 day

## plot the rel dispersion from L0 vs time
#mpl.rcParams['font.size'] = 35
#PGrid = (3,2)  
#fig = plt.figure(figsize=(52,40))
#ax_RD_0 = plt.subplot2grid(PGrid,(0,0))
#ax_RD_0.plot(RD_0.index[1:], RD_0['dr2_mean'][1:],
#              color='green', linewidth=3)
#ax_RD_0.set_xlabel('Date')
#ax_RD_0.set_ylabel('<Δr$^{2}$> (km$^2$)')
#ax_RD_0.xaxis.set_major_formatter(DateFormatter('%d-%m'))
#ax_RD_0.xaxis.set_major_locator(DayLocator(interval=10))
#ax_RD_0.tick_params(which='major', length=10)
#ax_RD_0.tick_params(which='minor', length=7)
#ax_RD_0.set_yscale("log")
#ax_RD_0.legend(fontsize=27, loc=4)

# Third: calculate the deformation from L0 for each day
D_d1d2 = disp_0(d1d2['L_tau'], time, freq=freq)  
D_d1d3 = disp_0(d1d3['L_tau'], time, freq=freq) 
D_d2d3 = disp_0(d2d3['L_tau'], time, freq=freq) 

# create a list of delta_r2 for all pairs
D_list = [D_d1d2['Deform'], D_d1d3['Deform'], D_d2d3['Deform']]
#standard deformation from L0
stdD  = std_disp(D_list, time, nob=nob, freq=freq)       # 24 indices * 1 day

## plot the std_D vs time
#mpl.rcParams['font.size'] = 35
#PGrid = (3,2)  
#fig = plt.figure(figsize=(52,40))
#ax_std_D = plt.subplot2grid(PGrid,(0,0))
#ax_std_D.plot(stdD.index[1:], stdD['std_D'][1:],
#              color='green', linewidth=3)
#ax_std_D.set_xlabel('Date')
#ax_std_D.set_ylabel('σ$_{Δr/L}$')
#ax_std_D.xaxis.set_major_formatter(DateFormatter('%d-%m'))
#ax_std_D.xaxis.set_major_locator(DayLocator(interval=10))
##ax_std_D.set_yscale("log")

##############################################################
#%% Calculate the two-particle dispersion and deformation rate
#   This based on the change between each time step
##############################################################

dr2_d1d2 = delta_r_sq(drifter_1, drifter_2, time, freq=freq)
dr2_d1d3 = delta_r_sq(drifter_1, drifter_3, time, freq=freq) 
dr2_d2d3 = delta_r_sq(drifter_2, drifter_3, time, freq=freq) 

# create a list of delta_r2 for all pairs
dr2_list = [dr2_d1d2['dr2'], dr2_d1d3['dr2'], dr2_d2d3['dr2']]
RD    = rel_disp(dr2_list, time, nob=nob, freq=freq)

# plot ensemble mean vs time
mpl.rcParams['font.size'] = 35
PGrid = (3,2)  
fig = plt.figure(figsize=(52,40))
ax_dr2 = plt.subplot2grid(PGrid,(0,0))
ax_dr2.plot(RD.index[freq:], RD['dr2_mean'][freq:], color='green', 
            linewidth=3)
ax_dr2.set_xlabel('Date')
ax_dr2.set_ylabel('<Δr$^{2}$> (km$^2$)')
ax_dr2.set_yscale("log")
ax_dr2.tick_params(which='major', length=10)
ax_dr2.tick_params(which='minor', length=7)
ax_dr2.xaxis.set_major_formatter(DateFormatter('%d-%m'))
ax_dr2.xaxis.set_major_locator(DayLocator(interval=10))

# The deformation/dispersion rate
Dr_d1d2 = disp_rate(d1d2['L_tau'], time, freq=freq)  
Dr_d1d3 = disp_rate(d1d3['L_tau'], time, freq=freq) 
Dr_d2d3 = disp_rate(d2d3['L_tau'], time, freq=freq) 

# create a list of D rate for all pairs
Dr_list   = [Dr_d1d2['D_rate'], Dr_d1d3['D_rate'], Dr_d2d3['D_rate']]
stdDr  = std_disp_rate(Dr_list, time, nob=nob, freq=freq)       # 24 indices * 1 day

# plot the std_D vs time
mpl.rcParams['font.size'] = 35
PGrid = (1,1)
fig = plt.figure(figsize=(20,12))
ax_std_Dr = plt.subplot2grid(PGrid,(0,0))
ax_std_Dr.plot(stdDr.index[freq:], stdDr['std_D'][freq:],
              color='green', linewidth=3)
ax_std_Dr.set_xlabel('Date')
ax_std_Dr.set_yscale("log")
ax_std_Dr.set_ylabel("σ$_\dot{D}$"' ($day^{-1}$)')
ax_std_Dr.xaxis.set_major_formatter(DateFormatter('%d-%m'))
ax_std_Dr.xaxis.set_major_locator(DayLocator(interval=10))
ax_std_Dr.tick_params(which='major', length=15)
ax_std_Dr.tick_params(which='minor', length=10)


######################################
#%% Relative dispersion statistics 
# following Lukovich and others (2017)
# this is a varaiance of mean squared 
# change in separation - R^2
######################################

# create a list of all zonal and meridional separations
Lx_list = [d1d2['Lx_tau'], d1d3['Lx_tau'], d2d3['Lx_tau']]
Ly_list = [d1d2['Ly_tau'], d1d3['Ly_tau'], d2d3['Ly_tau']]

nob = 3   # number of buoy pairs
R2  = rel_disp_var(Lx_list, Ly_list, time, nob=nob)

# plot R2 vs time
mpl.rcParams['font.size'] = 35
PGrid = (1,1)
fig = plt.figure(figsize=(20,12))
ax_R2 = plt.subplot2grid(PGrid,(0,0))
ax_R2.plot(R2.index, R2['R2_z'], color ='darkviolet', linewidth=3, label="Zonal") 
ax_R2.plot(R2.index, R2['R2_m'], color ='magenta', linewidth=3, label="Meridional")
ax_R2.plot(R2.index, R2['R2_t'], color='darkturquoise', linewidth=3,label="Total") 
ax_R2.set_xlabel('Date')
ax_R2.set_ylabel('RD$^2$ (km$^2$)')
ax_R2.set_yscale("log")
ax_R2.tick_params(which='major', length=15)
ax_R2.tick_params(which='minor', length=10)
ax_R2.xaxis.set_major_formatter(DateFormatter('%d-%m'))
ax_R2.xaxis.set_major_locator(DayLocator(interval=10))
ax_R2.legend(fontsize=30, loc=4)

#%% Deformation PSD
from scipy import signal
from scipy.signal import find_peaks
import matplotlib.ticker as ticker

d_rate = stdDr                                   

#%% Apply a Fast Fourier Transform to change the domain of the signal from the original
# (time or space domain) to a represenation in the frequency domain. This determines how 
# much variance (power) is contained in the freuqncy bands, which can be associated with 
# external forcing mechanims (e.g. wind forcing).
t = d_rate.index
x = d_rate['std_D'][1:]                                # cluster std_D
dt = 1                                                 # sampling time interval of drifter in hours

NFFT = np.size(x)                                      # FFT length (radix 2 number!)
m    = NFFT/2                                          # number of distinct frequency bins
fc   = 1/(2*dt)                                        # the critical frequency
f    = fc*np.arange(0,int(m)+1)/m                      # the frequency bins (NFFT/2 +1 bins)
X    = np.fft.fft(x,NFFT)                              # FFT of the length NFFT

Pxx  = np.zeros(np.size(x))
for j in range(len(X)):
    Pxx[j] = abs(X[j]*np.conj(X[j]))/NFFT              # periodogram = |X|^2       
Pxx  = Pxx[0:int(m)+1]                                 # from 0:m+1 
Pxx  = np.append(Pxx[0],Pxx[1:int(m)+1]*2)             # from 1:m+1 
print((x**2).sum())                                    # Parseval's theorem -total power computed in time domain
print(Pxx.sum())                                       # must be equal to total power computed in frequency domain
  
#%% Plot of Frequency vs Power 
mpl.rcParams['font.size'] = 25
fig, ax = plt.subplots(figsize=(20,10))
   
def invert(f):
    return (1/f) 

ax.loglog(f[1:int(m)],Pxx[1:int(m)], linewidth=2.5, 
            color='royalblue')   
ax.set_xlabel('Frequency (h$^{-1}$)')
ax.set_ylabel('Power (m$^{2}$ s$^{-2}$)')
# set the period (1/f) to be the secondary xaxis
ax2 = ax.twiny() 
ax2.loglog(1/f[1:int(m)],Pxx[1:int(m)],
               color='royalblue', linewidth=0.1)
ax2.invert_xaxis()
ax2.set_xscale(value='log')
ax2.set_xlabel('Period (hours)')
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax.grid(True, which='both')
ax.legend()

#%% Determine the PSD using Welch's method
# Here the Hann window is used in the Welch spectral object as the data are 
# in the form of a discretised time series. The window dampens the beginning 
# and end of the time series to 0 and better defines the spectral peaks. 

Fs = 1/dt                                              # sampling frequency (1/dt)
# cluster std_D
freqs_x, psd_x = signal.welch(x,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0)     # extracting frequencies and psd of signal

#%% Plot Frequency vs PSD
mpl.rcParams['font.size'] = 35
fig, ax = plt.subplots(figsize=(25,12))
ax.loglog(freqs_x,psd_x, linewidth=3, color='green')
ax.set_xlabel('Frequency (h⁻¹)')
ax.set_ylabel('PSD (m² s⁻² h⁻¹)')
ax.tick_params(which='major', length=10)
ax.tick_params(which='minor', length=4)
ax.legend(loc=3, fontsize=30)

# set the period (1/f) to be the secondary xaxis
ax2 = ax.twiny() 
ax2.loglog(1/freqs_x,psd_x, color='royalblue', linewidth=0.0001)
ax2.invert_xaxis()
ax2.set_xscale(value='log')
ax2.set_xlabel('Period (hours)')
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))


peaksx,freqs_peakx = find_peaks(psd_x, height=None,threshold=0.000005)
np.diff(peaksx)
periodx = 1/freqs_x[peaksx]                            # psd freq and peaks
xc = 1/periodx[12]                                     # change no. for line at inertial period
#ax.axvline(x=xc, color='red', linewidth=3)
ax.set_ylim(1e-7, 1e0)
ax2.grid(True, which='both')
ax.yaxis.grid(which='both')
ax.set_title('PSD σ$_\dot{D}$')

# end of code

