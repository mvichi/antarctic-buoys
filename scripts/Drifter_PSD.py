#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 12:09:18 2021

@author: Ashleigh Womack

Power Spectral Analysis of the Drifter Velocity
    - Power spectrum of a time series decribes the distribution of power with respect
      to the frequency contained in the analysed signal.
    - Adapted from the book: Modelling Methods for Marine Sciences 
      Glover, D.M., Jenkins, W.J. and Doney, S.C., 2011. Modeling methods for marine science. 
      Cambridge University Press.

"""
#%% 
from scipy import signal
from scipy.signal import find_peaks
from scipy.stats import chi2
import matplotlib.ticker as ticker

#%% First read in the velocity variables from Paper_1_diagnostics
# If the time interval of the data is irregular - resample the data to the best 
# time interval so that a constant sampling frequncy can be applied for the FFT
drifter = drifter.resample('4H').mean()

#%% Apply a Fast Fourier Transform to change the domain of the signal from the original
# (time or space domain) to a represenation in the frequency domain. This determines how 
# much variance (power) is contained in the freuqncy bands, which can be associated with 
# external forcing mechanims (e.g. wind forcing).
t = drifter.index
x = drifter.u[1:]                                       # zonal velocity
y = drifter.v[1:]                                       # meridional velocity
dt = 1   # change this                                  # sampling time interval of drifter in hours

# FFT for zonal component
NFFT = np.size(x)                                       # FFT length (radix 2 number!)
m = NFFT/2                                              # number of distinct frequency bins
fc = 1/(2*dt)                                           # the critical frequency
f = fc*np.arange(0,int(m)+1)/m                          # the frequency bins (NFFT/2 +1 bins)
X = np.fft.fft(x,NFFT)                                  # FFT of the length NFFT

Pxx = np.zeros(np.size(x))
for j in range(len(X)):
    Pxx[j] = abs(X[j]*np.conj(X[j]))/NFFT               # periodogram = |X|^2       
Pxx = Pxx[0:int(m)+1]                                   # from 0:m+1 
Pxx=np.append(Pxx[0],Pxx[1:int(m)+1]*2)                 # from 1:m+1 
print((x**2).sum())                                     # Parseval's theorem -total power computed in time domain
print(Pxx.sum())                                        # must be equal to total power computed in frequency domain
 
# FFT for meridional component
NFFT = np.size(y)                                       # FFT length 
m = NFFT/2                                              # number of discrete frequency bins
fc = 1/(2*dt)                                           # the critical frequency
f = fc*np.arange(0,int(m)+1)/m                          # the frequency bins (NFFT/2 +1 bins)
Y = np.fft.fft(y,NFFT)                                  # FFT of the length NFFT

Pyy = np.zeros(np.size(y))
for j in range(len(Y)):
    Pyy[j] = abs(Y[j]*np.conj(Y[j]))/NFFT               # periodogram = |X|^2   
Pyy = Pyy[0:int(m)+1]                                   # from 0:m+1 
Pyy=np.append(Pyy[0],Pyy[1:int(m)+1]*2)                 # from 1:m+1 
print((y**2).sum())                                     # Parseval's theorem -total power computed in time domain
print(Pyy.sum())                                        # must be equal to total power computed in frequency domain

f_linear = 10**f    
def invert(f):
    return (1/f)  

#%% Plot of Frequency vs Power 
mpl.rcParams['font.size'] = 25
fig, ax = plt.subplots(figsize=(20,10))
   
# Zonal component
freqs = np.fft.fftfreq(len(x))                          # find the frequency peaks 
peaks,freqs = find_peaks(Pxx, height=None,threshold=0.005)
np.diff(peaks)
period = 1/f[peaks]    

# Meridional component
freqs_y = np.fft.fftfreq(len(y))                        # find the frequency peaks 
peaks_y,freqs_y = find_peaks(Pyy, height=None,threshold=0.005)
np.diff(peaks_y)
period_y = 1/f[peaks_y]

ax.loglog(f[1:int(m)],Pxx[1:int(m)], linewidth=2.5, color='royalblue', label='zonal')   
ax.loglog(f[1:int(m)],Pyy[1:int(m)], linewidth=2.5, color='orange', label='meridional') 
ax.set_xlabel('Frequency (h$^{-1}$)')
ax.set_ylabel('Power (m$^{2}$ s$^{-2}$)') 
secax = ax.secondary_xaxis('top',functions=(invert, invert))
secax.set_xscale('linear')
secax.set_xlabel('Period (hours)')
ax.grid(True, which='both')
ax.legend()

#%% Determine the PSD using Welch's method
# Here the Hann window is used in the Welch spectral object as the data are 
# in the form of a discretised time series. The window dampens the beginning 
# and end of the time series to 0 and better defines the spectral peaks. 

Fs = 1/dt                                                    # sampling frequency (1/dt)
# Zonal component
freqs, psd = signal.welch(x,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0)           # extracting frequencies and psd of signal

# Meridional component
freqs_y, psd_y = signal.welch(y,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0)           # extracting frequencies and psd of signal

#%% Plot Frequency vs PSD
mpl.rcParams['font.size'] = 25
fig, ax = plt.subplots(figsize=(20,10))
ax.loglog(freqs,psd, linewidth=2.5, color='royalblue', label='zonal')
ax.loglog(freqs_y,psd_y, linewidth=2.5, color='orange', label='meridional')
ax.set_xlabel('Frequency (h⁻¹)')
ax.set_ylabel('PSD (m² s⁻² h⁻¹)')

secax = ax.secondary_xaxis('top',functions=(invert, invert))
secax.set_xscale(value='linear')
secax.set_xlabel('Period (hours)')
secax.grid()
ax.legend(loc=3)

peaks,freq = find_peaks(psd, height=None,threshold=0.0005)
np.diff(peaks)
period = 1/freqs[peaks]
xc = 1/period[12]                                             # change no. for line at inertial period
plt.axvline(x=xc, color='red', linewidth=2.5)
ax.grid(True, which='both')

#%% Application of Butterworth High-Pass Filter 
# This is a 12th order filter used to remove all frequencies below the cutoff
# frequency of 0.04 Hz which is approximately the daily frequency and allows 
# the sub-daily frequencies to become more significant. 

#dt = 1                                                     # sampling interval in hours
#Fs = 1/dt                                                  # sampling frequency (1/dt)
b, a = signal.butter(12, 0.04, 'high', Fs)                  # 12th order with 0.04 h-1 cutoff frequency            
w,h = signal.freqs(b,a)                                     # frequency response of a digital filter
x_filt = signal.filtfilt(b,a,drifter.u[1:], padtype=None)   # filter data series
y_filt = signal.filtfilt(b,a,drifter.v[1:], padtype=None)   # filter data series

#%% Apply a Fast Fourier Transform to filtered data series
# FFT for x_filt (zonal component)
NFFT = np.size(x_filt)                                      # FFT length (radix 2 number!)
m = NFFT/2                                                  # number of discrete frequency bins
fc = 1/(2*dt)                                               # the critical frequency
f = fc*np.arange(0,int(m)+1)/m                              # the frequency bins (NFFT/2 +1 bins)
X_filt = np.fft.fft(x_filt,NFFT)                            # FFT of the length NFFT

Pxx_filt = np.zeros(np.size(x_filt))
for j in range(len(X_filt)):
    Pxx_filt[j] = abs(X_filt[j]*np.conj(X_filt[j]))/NFFT    # periodogram = |X|^2   
Pxx_filt = Pxx_filt[0:int(m)+1]                             # from 0:m+1 
Pxx_filt = np.append(Pxx_filt[0],Pxx_filt[1:int(m)+1]*2)    # from 1:m+1
print((x_filt**2).sum())                                    # Parseval's theorem -total power computed in time domain
print(Pxx_filt.sum())                                       # must be equal to total power computed in frequency domain
freqs_filt = np.fft.fftfreq(len(x_filt))

# FFT for y_filt (meridional component)
NFFT = np.size(y_filt)                                      # FFT length (radix 2 number!)
m = NFFT/2                                                  # number of discrete frequency bins
fc = 1/(2*dt)                                               # the critical frequency
f = fc*np.arange(0,int(m)+1)/m                              # the frequency bins (NFFT/2 +1 bins)
Y_filt = np.fft.fft(y_filt,NFFT)                            # FFT of the length NFFT

Pyy_filt = np.zeros(np.size(y_filt))
for j in range(len(Y_filt)):
    Pyy_filt[j] = abs(Y_filt[j]*np.conj(Y_filt[j]))/NFFT    # periodogram = |X|^2   
Pyy_filt = Pyy_filt[0:int(m)+1]                             # from 0:m+1 
Pyy_filt = np.append(Pyy_filt[0],Pyy_filt[1:int(m)+1]*2)    # from 1:m+1 
print((y_filt**2).sum())                                    # Parseval's theorem -total power computed in time domain
print(Pyy_filt.sum())                                       # must be equal to total power computed in frequency domain
freqs_filt = np.fft.fftfreq(len(y_filt))

#%% Determine the PSD of filtered data series
# Zonal component
freqs_filt, psd_filt = signal.welch(x_filt,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0) 
# Meridional component
freqs_filt_y, psd_filt_y = signal.welch(y_filt,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0) 

#%% Linear PSD vs with High-Pass Filter
mpl.rcParams['font.size'] = 25
fig, ax = plt.subplots(figsize=(20,10))
ax.plot(freqs,psd, linewidth=3, label='normal')             # this just for the U component
ax.plot(freqs_filt,psd_filt, linewidth=3, label='filtered')
plt.plot(w,abs(h), color='magenta', label='frequency filter response')
plt.xlim(0,0.25)
ax.set_xlabel('Frequency (h⁻¹)')
ax.axis([0, 0.13, 0, 6])
ax.set_ylabel('PSD (m² s⁻² h⁻¹)')
plt.legend()

# end of code
