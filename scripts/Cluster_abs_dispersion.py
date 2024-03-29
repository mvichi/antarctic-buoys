#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 10:55:49 2022
Reads preprocessed files
Calculation of the Cluster Absolute (Single-Particle) Dispersion
      
@author: Ashleigh Womack (UCT)
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import math
from Drifter import coordinates
from datetime import timedelta
from matplotlib.dates import DayLocator,HourLocator,DateFormatter

BDIR = '../data/' # directory of the drifter data

#%% Cluster absolute dispersion  
def cluster_dispersion(dx_list, dy_list, time):
    '''
    Particle dispersion adapted from Taylor (1921) and Lukovich et al. (2017, 2021)
    Computes the variance of displacement around the cluster centroid, 
    as an indicator of dispersion.
    
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
    nob = len(dx_list) # number of buoys in the cluster
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


#%% Read in the buoy drifter data 
theBuoy = 'ISVP1'                                     
drifter_1 = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True) 

theBuoy = 'ISVP2'                                     
drifter_2 = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)

theBuoy = 'ISVP3'                                     
drifter_3 = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
   
# Manually resample the data for the cluster dates for every hour
#   or to the same time interval.
#   This should be the longest time frame where all buoys drifted together.
time = pd.date_range(start="2019-07-28 05:00:00", 
                         end="2019-08-25 03:00:00", freq="1H")

# find the index of the nearest value in an array
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

#%% Compute the new distance from new origin as their start time 
#   is different (cluster time frame)
drifter_1 = coordinates(drifter_1)
drifter_2 = coordinates(drifter_2)
drifter_3 = coordinates(drifter_3)

# create a list of all distances from origin (x and y)
dx_list = [drifter_1['x (km)'], drifter_2['x (km)'], drifter_3['x (km)']]
dy_list = [drifter_1['y (km)'], drifter_2['y (km)'], drifter_3['y (km)']]


#%% Compute and plot the cluster absolute dispersion
abs_disp = cluster_dispersion(dx_list, dy_list, time)

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
ax_A2.tick_params(which='major', length=10)
ax_A2.tick_params(which='minor', length=7)
ax_A2.xaxis.set_major_formatter(DateFormatter('%d-%m'))
ax_A2.xaxis.set_major_locator(DayLocator(interval=10))
ax_A2.legend(fontsize=35, ncol=1, loc=4) 

# end of code
