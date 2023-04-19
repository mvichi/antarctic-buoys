#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 12:09:18 2021

@author: Ashleigh Womack

Power Spectral Analysis of the Deformation rate of the cluster
    - Power spectrum of a time series decribes the distribution of power with respect
      to the frequency contained in the analysed signal.
    - Adapted from the book: Modelling Methods for Marine Sciences 
      Glover, D.M., Jenkins, W.J. and Doney, S.C., 2011. Modeling methods for marine science. 
      Cambridge University Press.

First - Run compute the relative_dispersion.py of the cluster to use the deformation proxy.

"""
from scipy import signal
from scipy.signal import find_peaks
from scipy.stats import chi2
import matplotlib.ticker as ticker

#%% 
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
