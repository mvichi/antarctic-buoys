#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 16:48:38 2023

@author: Ashleigh Womack

This script takes the 0 % SIC computed from matlab to
calculate the daily distance of drifter (buoy) from 
the ice edge.
"""

#%% Read in all packages 
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from matplotlib.dates import DayLocator,HourLocator,DateFormatter, MonthLocator
from geopy.distance import geodesic, distance
from datetime import timedelta, datetime, date
from math import radians, cos, sin, asin, sqrt, atan2, pi
from scipy.io import loadmat

######################################################
#%% Distance of drifter from sea ice edge - 0% and 15%
######################################################
def find_ind(array, value):
    ''' 
    Find the matching index of SIC longitude value to
    the drifter longitude value.
    '''
    idx = np.nanargmin((np.abs(array - value)))
    return idx

def find_nearest(array, value):
        '''
        Find the matching SIC longitude value to the
        drifter longitude. 
        '''
        array = np.asarray(array)
        idx = np.nanargmin((np.abs(array - value)))
        return array[idx]

#for the AMSR2 SIC data
def AMSR2_dist(drifter, amsr2_lats, amsr2_lons):
    '''
    Find the AMSR2 ice edge (0% and 15%) latitude for the same (nearest)
    longitude as the drifter. From this, the latitudinal distance 
    between the drifter and ice edge can be calculated.
    Parameters
    ----------
    drifter : pandas DataFrame
        Drifter DataFrame
    lats : float, list or npdarray
        SIC Latitudinal location
    lons : float, list or npdarray
        SIC Longitudinal location
    Returns 
    -------
    drifter: : DataFrame
        The updated DataFrame containing the amrs2 SIC coordinates
    '''
    # drifter lats and lons
    lats = drifter['latitude (deg)']    # value each day
    lons = drifter['longitude (deg)']   # value each day
    N = len(lats)
    
    # open arrays for amsr2 lats and lons
    IL      = np.zeros(N)     # index of 0 or 15 % longitudes
    IE_lats = np.zeros(N)
    IE_lons = np.zeros(N)
    i = 0
    for i in range(N):
        IE_lons[i] = find_nearest(amsr2_lons[i],lons[i])        # find lons value 
        IL[i]      = find_ind(amsr2_lons[i],lons[i])            # find lons index
        IE_lats[i] = amsr2_lats[i][IL[i].astype(int)]           # find lats value 
        i = i+1   
    drifter['amsr2 latitude (deg)']     = IE_lats
    drifter['amsr2 longitude (deg)']    = IE_lons
    
    # Latitudinal distance ~111 km for each degree
    IE_dist = (lats - IE_lats)*111
    drifter['amsr2 edge dist (km)']  = IE_dist
    
    return drifter

#%% Read in the buoy drifter data 
BDIR_P = '/home/..' 

theBuoy = 'ISVP1'                   
drifter_1 = pd.read_csv(BDIR_P+theBuoy+'.csv', index_col='time', parse_dates=True)
drifter_1 = drifter_1.resample('1D').mean(numeric_only)
#%% Latitudinal distance of drifters from the AMSR2 ice edge

# Co-ordinates (contours) of the AMSR2 0% SIC ice edge 
amsr2_dir = '/home/..'
amsr2_file = loadmat(amsr2_dir+'AMSR2_ice_edge.mat', squeeze_me=True)  # change file name
amsr2_coords = amsr2_file['coord']
amsr2_lats = amsr2_coords[:,1,:] 
amsr2_lats[amsr2_lats==0]=np.nan
amsr2_lons = amsr2_coords[:,0,:]
amsr2_lons[amsr2_lons==0]=np.nan

# make sure the indices of the amsr2 lats and lons match the file numbers (dates)
# you want to iterate through. 
drifter_1 = AMRS2_dist(drifter_1, amsr2_lats[:], amsr2_lons[:], 
                       amsr2_lats15[:], amsr2_lons15[:]) 


# end of code