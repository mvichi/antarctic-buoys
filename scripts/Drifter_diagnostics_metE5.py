#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 13:13:49 2023

Reads preprocessed files and ERA5 data for computing the 
drift response to atmospheric forcing:
- Wind factor, turning angle 
- Correlation coefficients (Pearson and vector)
- residual current
- Power Spectral Analysis of the Drifter Velocity
  Adapted from the book: Modelling Methods for Marine Sciences 
  Glover, D.M., Jenkins, W.J. and Doney, S.C., 2011. Modeling methods for marine science. 
  Cambridge University Press.

This script requires ERA5 hourly data (on single levels).
It must be modified for other met data.
    The variables should include 
        - wind velocity (u,v)
        - 2-m air temperature
        - mslp
    for the same period and region as the drifter of interest.

@author: Ashleigh Womack, Marcello Vichi (UCT)

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy
import math 
import pandas as pd
import xarray as xr
from scipy import signal
from scipy.signal import find_peaks
from scipy.stats import chi2
import matplotlib.ticker as ticker
import matplotlib as mpl

#%%
# Function to correct wind speed to 10 m reference values
def W10(U,z,z0=0.0002,z0u=0.1):
    """
    This function adjusts the wind speed U measured at height z
    to the corrected wind speed Uc which would be indicated
    locally at 10 m above terrain with roughness z0.
    All heights are in m and speeds in m/s.
    z0u is the effective roughness length of the terrain
        upstream of the measurement station
    z0 is roughness length in the application (e.g. over the grid cell)

    World Meteorological Organization 2014
    Guide to Meteorological Instruments and Methods of Observation.
    Geneva, Switzerland,
    World Meteorological Organization, 1128p. (WMO-No.8, 2014)

    """
    CF = 1.  # flow distortion correction (1 if on a mast)
    CT = 1.  # correction factor due to topographic effects (1 on the ocean)
    corr = (np.log(60/z0u)*np.log(10/z0))/(np.log(10/z0u)*np.log(60/z0))
    Uc = U*CF*CT*(np.log(10/z0u)/np.log(z/z0u))*corr
    return Uc

# Function to compute magnitude and direction from u,v vectors
def uv2magdir(u,v):
    """
    Computes magnitude and direction from vector components
    """
    wm = np.sqrt(u**2+v**2)
    wd = 270 - (180/np.pi)*np.arctan2(v,u) # trig angle to wind direction
    wdcorr = np.where(wd<360,wd,wd-360)
    return wm,wdcorr

# Extract the e5 variables 
def e5_variables(drifter, e5_file):
    '''
    Extract the ERA5 variables at the location of the drifter
    
    Parameters
    ----------
    drifter : pandas DataFrame
        DESCRIPTION.
    e5_file : NetCDF4 file
        The file containing ERA5 mslp, 2m-Temp, 10m-winds (u, v)
   Returns
    -------
    drifter : DataFrame
        The updated DataFrame containing ERA5 reanalysis variables
    '''
    # create the variables
    U10     = np.zeros(len(drifter))
    U10_dir = np.zeros(len(drifter))
    temp    = np.zeros(len(drifter))
    mslp    = np.zeros(len(drifter))
    u10     = np.zeros(len(drifter))
    v10     = np.zeros(len(drifter))
    n = 0
    for i,j,t in zip(drifter['longitude (deg)'],drifter['latitude (deg)'],drifter.index):
        u       = e5_file["u10"].sel(longitude=i,latitude=j,time=t,method='nearest').values
        v       = e5_file["v10"].sel(longitude=i,latitude=j,time=t,method='nearest').values
        U10[n], U10_dir[n] = uv2magdir(u,v) 
        temp[n] = e5_file["t2m"].sel(longitude=i,latitude=j,time=t,method='nearest').values
        mslp[n] = e5_file["msl"].sel(longitude=i,latitude=j,time=t,method='nearest').values
        u10[n]  = e5_file["u10"].sel(longitude=i,latitude=j,time=t,method='nearest').values
        v10[n]  = e5_file["v10"].sel(longitude=i,latitude=j,time=t,method='nearest').values
        n=n+1
    # add E5 variables to drifter DataFrame
    drifter['U10_E5 (m/s)']      = U10         # modular speed of wind
    drifter['U10_dir_E5 (degs)'] = U10_dir     # wind direction 
    drifter['t2m_E5 (degC)']     = temp-273.15 # Kelvin to Celcius
    drifter['mslp_E5 (hPa)']     = mslp/100.   # Pa to hPa
    drifter['u10_E5 (m/s)']      = u10         # u-component of wind
    drifter['v10_E5 (m/s)']      = v10         # v-component of wind
    return drifter   

#%% Create quiver plot for wind and drifter speeds and direction
    
def stick_plot(time, u, v, **kw):
    width = kw.pop('width', 0.002)
    headwidth = kw.pop('headwidth', 0)
    headlength = kw.pop('headlength', 0)
    headaxislength = kw.pop('headaxislength', 0)
    angles = kw.pop('angles', 'uv')
    ax = kw.pop('ax', None)
    
    if angles != 'uv':
        raise AssertionError("Stickplot angles must be 'uv' so that"
                             "if *U*==*V* the angle of the arrow on"
                             "the plot is 45 degrees CCW from the *x*-axis.")

    time, u, v = map(np.asanyarray, (time, u, v))
    if not ax:
        fig, ax = plt.subplots()
    
    q = ax.quiver(mdates.date2num(time), [[0]*len(time)], u, v,
                  angles='uv', width=width, headwidth=headwidth,
                  headlength=headlength, headaxislength=headaxislength,
                  **kw)

    ax.axes.get_yaxis().set_visible(False)
    ax.xaxis_date()
    return q


#%% Compute the wind factor, turning angle and vector coefficient of determination   
def drift_response(drifter):
    '''
    Relationship between the ice motions and ERA5 winds - 
    adapted from Kimura and Wakatsuchi (2000), Kimura (2004), 
    Fukamachi et al. (2011)
    Computes the wind factor, turning angle and vector coefficient 
    of determination between the wind speed and ice (drifter) speed.
    We also compute the Pearson coefficient of determination.
    Parameters
    ----------
    drifter : pandas DataFrame
        Drifter DataFrame
    Returns
    -------
    theta_degs, WF, R2v, R2p , p-value: numpy.float64
    '''
    # variables (component - mean)
    u   = drifter['u (m/s)'] - np.mean(drifter['u (m/s)'])
    v   = drifter['v (m/s)'] - np.mean(drifter['v (m/s)'])
    u10 = drifter['u10_E5 (m/s)'] - np.mean(drifter['u10_E5 (m/s)'])
    v10 = drifter['v10_E5 (m/s)'] - np.mean(drifter['v10_E5 (m/s)'])
    
    # sum of multiplication of variables
    u10v  = np.nansum(u10*v)
    v10u  = np.nansum(v10*u)
    u10u  = np.nansum(u10*u)
    v10v  = np.nansum(v10*v)
    u10_2 = np.nansum(np.square(u10)) # u10^2
    v10_2 = np.nansum(np.square(v10))
    u_2   = np.nansum(np.square(u)) 
    v_2   = np.nansum(np.square(v))
    
    # turning angle (angle from the wind to ice drift)
    # it is positive on a cartesian coordinate system
    # but needs to be negative on geographical coordinate system
    #theta = -math.atan((u10v - v10u)/(u10u + v10v)) 
    theta = -math.atan2(u10v - v10u, u10u + v10v) 
    theta_degs = math.degrees(theta) # in degrees (neg due to left deflection)
    
    # wind factor     
    c1 = (np.cos(theta))*(u10u)
    c2 = (np.sin(theta))*(v10u)
    c3 = (np.sin(theta))*(u10v)
    c4 = (np.cos(theta))*(v10v)

    WF = ((c1 + c2 - c3 + c4)/(u10_2 + v10_2))  
    
    # vector correlation
    R2v = np.square((c1 + c2 - c3 + c4)/(np.sqrt(u10_2 + v10_2)*np.sqrt(u_2 + v_2)))

    # Pearson correlation
    r, p = scipy.stats.pearsonr(drifter['U10_E5 (m/s)'][1:], 
                                                 drifter['U (m/s)'][1:])
    R2p = np.square(r)
    
    return theta, theta_degs, WF, R2v, R2p, p

#%% Use the drift response measurements to estimate the underlying currents
def current_estimation(drifter,theta, WF):
    '''
    Estimate the underlying current velocity components -
    This method is taken from Kimura and Wakatsuchi (2000).
    It relates sea-ice velocity (u,v), ERA5 wind velocity 
    (u10,v10) for the same time, the wind factor, turning
    angle and mean ocean current (Cu,Cv).
    Parameters
    ----------
    drifter : pandas DataFrame
        Drifter DataFrame
    theta_degs : numpy.float64
        Mean turning angle in degrees
    WF : numpy.float64
        Mean wind factor
    Returns
    -------
    Cu, Cv : numpy.float64
        Mean ocean current components (m/s)
    C : numpy.float64
        Mean modular current speed (m/s)
    '''
    # mean of the velocity components
    u10 = np.mean(drifter['u10_E5 (m/s)']) # mean u of wind
    v10 = np.mean(drifter['v10_E5 (m/s)']) # mean v of wind
    u   = np.mean(drifter['u (m/s)'])      # mean u of ice
    v   = np.mean(drifter['v (m/s)'])      # mean v of ice
    
    # Estimate the mean ocean current which is derived by subtracting 
    # the wind effect from the sea-ice motion
    # set up the three matrices
    M1 = np.array([[u],[v]])
    # apply the counter-clockwise rotation matrix
    M2 = np.array([[np.cos(theta), np.sin(theta)],
           [-np.sin(theta), np.cos(theta)]])
    M3 = np.array([[u10],[v10]])
    
    # this is the dot product multiplication 
    Cuv = (M1 - WF*(M2.dot(M3)))
    C     = np.sqrt(Cuv[0]**2 + Cuv[1]**2)
    
    return Cuv, C




#%% USER INPUT
BDIR = '../data/' # directory of the drifter data
EDIR = './'# directory of ERA5 netcdf file

# Read in the buoy drifter data 
theBuoy = 'ISVP1'                                     
drifter = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
# dates ticks to be used in the final figure
DAY = 10
E5 = '../ERA5_2019.nc'

#### Other examples (uncomment)
# Winter 2017  
# theBuoy = 'Trident_U1'                                     
# drifter = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
# DAY = 15
# drifter = drifter.loc['2017-08-05':]     # locate any special date range           
 
# theBuoy = '2014S9'                                     
# drifter = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
# DAY = 10
# drifter = drifter.loc['2015-07-06':]      

# theBuoy = '2015P6'                                     
# drifter = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
# DAY = 15
# df=drifter.loc['2015-04-01':'2015-08-15'].copy()
# # quality-check the data using the max criterion (< 2 m/s)
# drifter=drifter_check(df,criterion='max',maxspeed=2.) 

#%% extract coordinates (not used yet)
lats = drifter['latitude (deg)']
lons = drifter['longitude (deg)']
x = drifter['x (km)']
y = drifter['y (km)']
dfo = drifter['distance (km)'] # distance from origin
start = drifter.index.min()
end = drifter.index.max()

#%% Compute the drift response estimates
# Read in the ERA5 data- atmospheric reanalysis  
e5_file = xr.open_dataset(EDIR+E5)

#   Extract the ERA5 atmospheric data at drifter location
drifter = e5_variables(drifter, e5_file)

# compute the kinematic parameters
theta, theta_degs, WF, R2v, R2p = drift_response(drifter)    
# using the WF and theta, compute the residual current speed (m/s)
Cuv, C = current_estimation(drifter,theta, WF)

#%% Define workspace velocity variables 
#   If the time interval of the data is irregular - resample the data to the best 
#   time interval so that a constant sampling frequncy can be applied for the FFT
      
t = drifter.index
x = drifter['u (m/s)'][1:]                                       # zonal velocity
y = drifter['v (m/s)'][1:]                                       # meridional velocity
# NB: HARDCODED VALUE
dt = 1   # change this                                           # sampling time interval of drifter in hours

# Era5 wind u and v components
xE5 = drifter['u10_E5 (m/s)'][1:]                                # Era5 zonal velocity
yE5 = drifter['v10_E5 (m/s)'][1:]                                # Era5 meridional velocity
#%% Apply a Fast Fourier Transform to change the domain of the signal from the original
#   (time or space domain) to a represenation in the frequency domain. This determines how 
#   much variance (power) is contained in the freuqncy bands, which can be associated with 
#   external forcing mechanims (e.g. wind forcing).

# FFT for zonal component
NFFT = np.size(x)                                                # FFT length (radix 2 number!)
m    = NFFT/2                                                    # number of distinct frequency bins
fc   = 1/(2*dt)                                                  # the critical frequency
f    = fc*np.arange(0,int(m)+1)/m                                # the frequency bins (NFFT/2 +1 bins)
X    = np.fft.fft(x,NFFT)                                        # FFT of the length NFFT

Pxx  = np.zeros(np.size(x))
for j in range(len(X)):
    Pxx[j] = abs(X[j]*np.conj(X[j]))/NFFT                        # periodogram = |X|^2       
Pxx  = Pxx[0:int(m)+1]                                           # from 0:m+1 
Pxx  = np.append(Pxx[0],Pxx[1:int(m)+1]*2)                       # from 1:m+1 
print((x**2).sum())                                              # Parseval's theorem -total power computed in time domain
print(Pxx.sum())                                                 # must be equal to total power computed in frequency domain
 
# FFT for meridional component
NFFT = np.size(y)                                                # FFT length 
m    = NFFT/2                                                    # number of discrete frequency bins
fc   = 1/(2*dt)                                                  # the critical frequency
f    = fc*np.arange(0,int(m)+1)/m                                # the frequency bins (NFFT/2 +1 bins)
Y    = np.fft.fft(y,NFFT)                                        # FFT of the length NFFT

Pyy  = np.zeros(np.size(y))
for j in range(len(Y)):
    Pyy[j] = abs(Y[j]*np.conj(Y[j]))/NFFT                        # periodogram = |X|^2   
Pyy  = Pyy[0:int(m)+1]                                           # from 0:m+1 
Pyy  = np.append(Pyy[0],Pyy[1:int(m)+1]*2)                       # from 1:m+1 
print((y**2).sum())                                              # Parseval's theorem -total power computed in time domain
print(Pyy.sum())                                                 # must be equal to total power computed in frequency domain
 
f_linear = 10**f      

#%% Plot of Frequency vs Power 
mpl.rcParams['font.size'] = 25
fig, ax = plt.subplots(figsize=(20,10))
   
def invert(f):
    return (1/f) 

ax.loglog(f[1:int(m)],Pxx[1:int(m)], linewidth=2.5, color='royalblue', label='zonal')   
ax.loglog(f[1:int(m)],Pyy[1:int(m)], linewidth=2.5, color='orange', label='meridional') 
ax.set_xlabel('Frequency (h$^{-1}$)')
ax.set_ylabel('Power (m$^{2}$ s$^{-2}$)')
# set the period (1/f) to be the secondary xaxis
ax2 = ax.twiny() 
ax2.loglog(1/f[1:int(m)],Pxx[1:int(m)], color='royalblue', linewidth=0.1)
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

Fs = 1/dt                                                    # sampling frequency (1/dt)
# Zonal component
freqs_x, psd_x = signal.welch(x,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0)           # extracting frequencies and psd of signal
freqs_xE5, psd_xE5 = signal.welch(xE5,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0)           

# Meridional component
freqs_y, psd_y = signal.welch(y,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0)           
freqs_yE5, psd_yE5 = signal.welch(yE5,fs=Fs, window='hann',scaling='density', 
                          nperseg=256, noverlap=0)           


#%% Plot Frequency vs PSD
mpl.rcParams['font.size'] = 35
fig, ax = plt.subplots(figsize=(25,12))
ax.loglog(freqs_x,psd_x, linewidth=2.5, color='royalblue', label='zonal drift')
ax.loglog(freqs_y,psd_y, linewidth=2.5, color='orange', label='meridional drift')
ax.loglog(freqs_xE5,psd_xE5, linewidth=2.5, color='blue', label='zonal wind')
ax.loglog(freqs_yE5,psd_yE5, linewidth=2.5, color='orangered', label='meridional wind')
ax.set_xlabel('Frequency (h⁻¹)')
ax.set_ylabel('PSD (m² s⁻² h⁻¹)')
ax.tick_params(which='major', length=10)
ax.tick_params(which='minor', length=4)
ax.legend(loc=3, fontsize=30)

# Calculate and the plot the 95% confidence intervals of the PSD
probability = 0.95 
alfa = 1 - probability
# P is the number of estimates in welch function and 
# also the degree of freedom (256).
P = 256
v = 2 * P 
ci = chi2.ppf([1 - alfa/2, alfa/2], v)
ci = v / ci
# drifter confidence
Pxxc_lower_x = psd_x * ci[0]
Pxxc_upper_x = psd_x * ci[1]
Pxxc_lower_y = psd_y * ci[0]
Pxxc_upper_y = psd_y * ci[1]
plt.fill_between(freqs_x, Pxxc_lower_x, Pxxc_upper_x, color = 'lightblue')
plt.fill_between(freqs_y, Pxxc_lower_y, Pxxc_upper_y, color = 'navajowhite')
# era5 confidence
Pxxc_lower_xE5 = psd_xE5 * ci[0]
Pxxc_upper_xE5 = psd_xE5 * ci[1]
Pxxc_lower_yE5 = psd_yE5 * ci[0]
Pxxc_upper_yE5 = psd_yE5 * ci[1]
plt.fill_between(freqs_xE5, Pxxc_lower_xE5, Pxxc_upper_xE5, color = 'cornflowerblue')
plt.fill_between(freqs_yE5, Pxxc_lower_yE5, Pxxc_upper_yE5, color = 'lightsalmon')

# set the period (1/f) to be the secondary xaxis
ax2 = ax.twiny() 
ax2.loglog(1/freqs_x,psd_x, color='royalblue', linewidth=0.0001)
ax2.invert_xaxis()
ax2.set_xscale(value='log')
ax2.set_xlabel('Period (hours)')
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))

peaks,freqs_peak = find_peaks(psd_x, height=None,threshold=0.0005)
np.diff(peaks)
period = 1/freqs_x[peaks]                                    # psd freq and peaks
xc = 1/period[12]                                            # change no. for line at inertial period
ax.axvline(x=xc, color='black', linewidth=3)
ax2.grid(True, which='both')
ax.yaxis.grid(which='both')

# end of code
