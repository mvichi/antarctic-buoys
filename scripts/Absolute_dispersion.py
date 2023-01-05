#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:10:03 2023

@author: Ashleigh Womack

Reads preprocessed files
Calculation of the Cluster Absolute (Single-Particle) Dispersion
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import matplotlib.dates as mdates
import seaborn as sns
from geopy.distance import distance 
from matplotlib.dates import DayLocator,HourLocator,DateFormatter

BDIR_P = '/home/ashleigh/Documents/PhD/data/buoys/Processed_buoys/' # directory of the processed drifter data

#%% Read in the buoy drifter data 
theBuoy = 'ISVP4'                                     
drifter_1 = pd.read_csv(BDIR_P+theBuoy+'.csv',index_col='time',parse_dates=True)
# 24/10/2019 12:00 - 07/12/2019  
# locate the correct timeframe 
drifter_1 = drifter_1.loc[:'2019-12-07']
#%% Read in the buoy drifter data 
theBuoy = 'ISVP5'                                     
drifter_2 = pd.read_csv(BDIR_P+theBuoy+'.csv',index_col='time',parse_dates=True)
# 28/10/2019 10:00 - 09/12/2019  
# locate the correct timeframe 
drifter_2 = drifter_2.loc[:'2019-12-09'] 
#%% Read in the buoy drifter data 
theBuoy = 'ISVP6'                                     
drifter_3 = pd.read_csv(BDIR_P+theBuoy+'.csv',index_col='time',parse_dates=True)
# 29/10/2019 12:00 - 07/12/2019
# locate the correct timeframe 
drifter_3 = drifter_3.loc[:'2019-12-07']
#%% Read in the buoy drifter data 
theBuoy = 'Trident_U4'                                     
drifter_4 = pd.read_csv(BDIR_P+theBuoy+'.csv',index_col='time',parse_dates=True)
# 30/10/2019 12:00 - 02/12/2019 00:00
# locate the correct timeframe 
drifter_4 = drifter_4.loc['2019-10-30 12:00:00':]

#%% Manually resample the data for the cluster dates
# for every hour
time = pd.date_range(start="2019-10-30 12:00:00", 
                         end="2019-12-02 00:00:00", freq="4H")

# find the nearest longitude from buoy to corresponding sea ice contour
def find_ind(array, value):
    idx = np.nanargmin((np.abs(array - value)))
    return idx

dates_1 = np.zeros(len(time))
dates_2 = np.zeros(len(time))
dates_3 = np.zeros(len(time))
dates_4 = np.zeros(len(time))
i=0
for i in range(len(time)):
    dates_1[i] = find_ind(drifter_1.index,time[i])
    dates_2[i] = find_ind(drifter_2.index,time[i])
    dates_3[i] = find_ind(drifter_3.index,time[i])
    dates_4[i] = find_ind(drifter_4.index,time[i])
    i=i+1

drifter_1 = drifter_1.iloc[dates_1] 
drifter_2 = drifter_2.iloc[dates_2] 
drifter_3 = drifter_3.iloc[dates_3]
drifter_4 = drifter_4.iloc[dates_4]
      
#%% compute the new distance from new origin as their start time 
#   is different (cluster time)
drifter_1 = coordinates(drifter_1)
drifter_2 = coordinates(drifter_2)
drifter_3 = coordinates(drifter_3)
drifter_4 = coordinates(drifter_4)

# create a list of all distances from origin (x and y)
dx_list = [drifter_1['x (km)'], drifter_2['x (km)'], drifter_3['x (km)'],
           drifter_4['x (km)']]
dy_list = [drifter_1['y (km)'], drifter_2['y (km)'], drifter_3['y (km)'],
           drifter_4['y (km)']]

#%% absolute dispersion
nob = 4   # number of buoys
def abs_dispersion(dx_list, dy_list, time, nob=nob):
    '''
    Single particle dispersion - 
    adapted from Taylor (1921) and Lukovich et al. (2017, 2021)
    Computes the variance of displacement around the cluster centroid, 
    an indicator of dispersion.
    The method should be computed on an ensemble of floats
    
    A^2 = < (xi(t) - xi(0) - <xi(t) - xi(0)>)**2 >  
    
    where xi is the [x or y] position of drifter i at time t, and the angular 
    brackets denote the ensemble mean over all particles. 
    Parameters
    ----------
    dx_list: list
        Contains the x-direction distances from origin
    dy_list: list
        Contains the y-direction distances from origin
    time: DatetimeIndex
        array of the hourly drifter dates
    nob:  int
        number of drifter pairs 
    Returns
    -------
    A2: Dataframe containing the absolute dispersion for the buoy cluster
    '''
    # time interval
    length = math.ceil(time.size)
    # tau -> change in time from t(0)
    tau    = [abs(time[d]-time[0])/timedelta(days=1) for d in range(length)]
    A2_z   = []
    A2_m   = []
    A2_t   = []
    n=0
    for n in range(length):
        # zonal(x) direction dispersion
        A2_z.append((1/nob)*np.nansum(np.square([l[n] for l in dx_list] 
            - ((1/nob)*(np.nansum([l[n] for l in dx_list], axis=0)))), axis=0))
        # meridional (y) direction disperion
        A2_m.append((1/nob)*np.nansum(np.square([l[n] for l in dy_list] 
            - ((1/nob)*(np.nansum([l[n] for l in dy_list], axis=0)))), axis=0))
        # total absolute dispersion (equal to the components combined)
        A2_t.append(A2_z[n] + A2_m[n])
        n=n+1
    A2 = pd.DataFrame({'time': time, 'delta_t': tau, 'A2_z': A2_z, 
                       'A2_m': A2_m, 'A2_t': A2_t}, 
                            columns=['time', 'delta_t', 'A2_z', 
                                     'A2_m', 'A2_t']) 
    A2.set_index('time',inplace=True)
    return A2

abs_disp  = abs_dispersion(dx_list, dy_list, time, nob=nob)
# remove and redo for isvps as the trident buoy as it did not move as north 
#abs_disp2 = abs_dispersion(dx_list2, dy_list2, time, nob=3)

#%% plot the absolute dispersion for the cluster
mpl.rcParams['font.size'] = 40
PGrid = (1,1)
fig = plt.figure(figsize=(27,14))
ax_A2 = plt.subplot2grid(PGrid, (0,0))
ax_A2.plot(abs_disp.index, abs_disp['A2_z'], color ='darkviolet', linewidth=3, label="Zonal") 
ax_A2.plot(abs_disp.index, abs_disp['A2_m'], color ='magenta', linewidth=3, label="Meridional")
ax_A2.plot(abs_disp.index, abs_disp['A2_t'], color='darkturquoise', linewidth=3,label="Total") 
ax_A2.set_xlabel('Date')
ax_A2.set_yscale("log")
ax_A2.set_ylabel('A$^2$ (km$^2$)')
ax_A2.tick_params(which='major', length=15)
ax_A2.tick_params(which='minor', length=10) 
ax_A2.xaxis.set_major_formatter(DateFormatter('%d-%m'))
ax_A2.xaxis.set_major_locator(DayLocator(interval=5))
ax_A2.legend(fontsize=35, ncol=1, loc=4) 
ax_A2.set_ylim(1e-3, 2e4) 