#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 11:52:19 2022

@author: Ashleigh Womack (UCT)

Functions to 
1. compute the relative dispersion of sea ice,
using methods reported by Rampal and others (2008, 2009)
and Girard and others (2009).

2. compute the relative dispersion of sea ice, using methods 
reported by Lukovich and others (2017, 2021). This is a variance
method.

"""
import numpy as np
import pandas as pd
import math
from geopy.distance import distance
from datetime import timedelta
#############################
#%% calculate the deformation
#############################

# 1. Calulate separation between buoy pairs
# 2. Create a list of all dr2
# 3. Calculate the mean squared change in separation from origin
# 4. Calculate the deformation for each buoy pair
# 5. Create a list of the deformations from origin
# 6. Calculate the standard deviation of deformation

def separation(drifter1, drifter2):
    '''
    Calculate the separation L between two drifters, where L(tau) is the 
    separation at time=tau and L0 is the separation at time=0 between
    buoys numbered 1 and 2 with absolute [x or y] positions
    x1 and x2, and with separtaions Y=x1-x2 in the.
    
    Then calculate squared change in separation (delta_r_squared) where 
    delta_r =  delta_r = ||Y(tau)|| - ||Y0|| = L(tau) - L0.
    
    Parameters
    ----------
    drifter1 : pandas DataFrame
        The drifter DataFrame
    drifter2 : pandas DataFrame
        The drifter DataFrame (shorter drifter)
    time: DateTimeIndex
    Returns
    -------
    drifter_pair : A DataFrame containing the L_tau and dr2_0
    '''
    # set the time to the length of the shortest buoy
    t    = drifter2.index
    lat1 = drifter1['latitude (deg)']
    lon1 = drifter1['longitude (deg)']
    lat2 = drifter2['latitude (deg)']
    lon2 = drifter2['longitude (deg)']
    # separation at each time step between the drifter pair -> L(tau)
    L_tau  = [distance([lat1[k], lon1[k]],
                      [lat2[k], lon2[k]]).km for k in range(len(lat2))]
    # zonal and meridional separation
    Lx_tau = [distance([lat1[k], lon1[k]],
                       [lat1[k], lon2[k]]).km for k in range(len(lat2))]
    Ly_tau = [distance([lat1[k], lon1[k]],
                       [lat2[k], lon1[k]]).km for k in range(len(lat2))]        
    # calculate the delta_r_sq from L0
    N = len(t)
    if N==0:
        dr2_0 = np.nan
    else: 
        dr2_0 = np.zeros(N)
        for k in range(N):
            dr2_0[k] = np.square(L_tau[k] - L_tau[0])     
    drifter_pair = pd.DataFrame({'time': t, 'Lx_tau': Lx_tau, 'Ly_tau': Ly_tau,
                                 'L_tau': L_tau, 'dr2_0': dr2_0}, 
                            columns=['time', 'Lx_tau', 'Ly_tau', 'L_tau', 'dr2_0']) 
    drifter_pair.set_index('time',inplace=True)
    return drifter_pair

freq = 24
nob  = 3 
def rel_disp_0(dr2_list_0, time, nob=nob, freq=freq):
    '''
    Compute the ensemble mean of the dr2_0 - squared change in separation
    from L0.
    Parameters
    ----------
    dr2_list_0: list
        list that holds all the dr2_0 for each drifter pair.
    time: DatetimeIndex
        array of the hourly drifter dates
    nob:  int
        number of drifter pairs
    freq: int
        the frequency of window over which calulation is run
    Returns
    -------
    RD_0:
        Dataframe with the relative disperionfrom L0 and corresponding time
    '''
    # time interval
    length = math.ceil(time.size)
    # tau -> change in time from t(0)
    tau    = [abs(time[d]-time[0])/timedelta(days=1) for d in range(length)]
    dr2_mean = []
    n=0
    for n in range((len(time))):
        dr2_mean.append((1/nob)*(np.nansum([l[n] for l in dr2_list_0], axis=0)))
    n=n+1
    RD_0 = pd.DataFrame({'time': time, 'delta_t': tau, 'dr2_mean': dr2_mean}, 
                            columns=['time', 'delta_t','dr2_mean']) 
    RD_0.set_index('time',inplace=True)
    return RD_0  

def disp_0(L_tau, time, freq=freq):
    '''
    Two-particle dispersion - 
    adapted from Rampal and others (2008); Girard and others (2009)
    Computes the dispersion rate between two drifters. From 
    a solid mechanics perspective, it is pertinent to consider the 
    rate D which is analogous deformation rate:
    
    D = delta_r/L0
    
    where delta_r =  delta_r = ||Y(tau)|| - ||Y0|| = L(tau) - L0
    
    Parameters
    ----------
    L_tau: numpy array
        array of distances between two drifters
    tau: numpy array
        array of number of days 
    Returns
    -------
    D: the dispersion rate 
    '''
    # time interval
    length = math.ceil(time.size)
    # tau -> change in time from t(0)
    tau    = [abs(time[d]-time[0])/timedelta(days=1) for d in range(length)]
    N = len(tau)
    if N==0:
        D = np.nan
    else: 
        D = np.zeros(N)
        for k in range(N):
            D[k] = (L_tau[k] - L_tau[0])/(L_tau[0])
    D = pd.DataFrame({'time': time,'Deform': D, 'delta_t': tau}, 
                            columns=['time','Deform', 'delta_t']) 
    D.set_index('time',inplace=True)
    return D

def std_disp(D_list, time, nob=nob, freq=freq):
    '''
    Standard deviation of Deformation D- 
    adapted Rampal and others (2008); Girard and others (2009)
    Computes the standard deviation of dispersion normalised by
    length and time:
        
    std_D = (< (D - < D >)^2 >)^0.5
    
    where D is the dispersion and the angular brackets denote the 
    ensemble mean caluclated over N pairs of buoys initally separated 
    by L0 +- L_tau and over a time interval tau.
    Parameters
    ----------
    D_list: list
        list that holds all the dr2 for each drifter pair.
    time: DatetimeIndex
        array of the hourly drifter dates
    nob:  int
        number of drifter pairs
    freq: int
        the frequency of window over which calulation is run
    Returns
    -------
    std_D: Dataframe containing the std of dispersion
    '''
    length = math.ceil(time.size)
    tau    = [abs(time[d]-time[0])/timedelta(days=1) for d in range(length)]
    sD     = []
    m=0
    for m in range(length):
        sD.append(np.sqrt((1/nob)*np.nansum(np.square([l[m] for l in D_list] 
            - ((1/nob)*(np.nansum([l[m] for l in D_list], axis=0)))), axis=0)))
        m=m+1
    std_D = pd.DataFrame({'time': time, 'delta_t': tau, 'std_D': sD}, 
                            columns=['time', 'delta_t', 'std_D']) 
    std_D.set_index('time',inplace=True)
    return std_D

##############################################################
#%% Calculate the two-particle dispersion and deformation rate
##############################################################
    
# 1. Use separation between buoy pairs from above
# 2. Calculate the delta_r_squared for each buoy pair (moving time)
# 3. Create a list of the new dr2
# 4. Calculate the mean squared change in separation
# 5. Calculate the dispersion/deformation rate
# 6. Create a list of the deformation rates
# 7. Calculate the standard deviation of the deformation rate

def delta_r_sq(drifter1, drifter2, time, freq=freq):
    '''
    Calculation of the change in separation between a pair of 
    buoys (delta_r) as a function of tau - Methods from Rampal and others (2009).
    
    delta_r^2 = np.square(||Y(t+tau)|| - ||Y(t)||) = np.square(L(t+tau) - L(t))
        
    Parameters
    ----------
    drifter1 : pandas DataFrame
        The drifter DataFrame
    drifter2 : pandas DataFrame
        The drifter DataFrame (shorter drifter)
    time : datetime
    freq: int
        How many indices in tau (change in time)
    Returns
    -------
    drifter_pair : A DataFrame containing the L_tau, delta_r and delta_r2
    '''
    # separation at each time step between the drifter pair -> L(tau)
    lat1  = drifter1['latitude (deg)']
    lon1  = drifter1['longitude (deg)']
    lat2  = drifter2['latitude (deg)']
    lon2  = drifter2['longitude (deg)']
    L_tau = [distance([lat1[k], lon1[k]],
                      [lat2[k], lon2[k]]).km for k in range(len(lat2))]
    # set the time and tau (delta_t)
    tau   = [abs((time[l+freq]-time[l]))/timedelta(days=1) for l in range(len(time)-freq)]
    # this time delta_r is a function of time
    # change in separation l(t+tau) - L(t) = delta_r(tau)
    N = len(tau)
    delta_r = np.zeros(N)
    for l in range(N):
        delta_r[l] = (L_tau[l+freq]-L_tau[l])
    # pad with freq len of zeros
    tau     = np.append(np.zeros(freq), tau)
    delta_r = np.append(np.zeros(freq), delta_r)
    # delta_r squared 
    dr2 = np.square(delta_r)
    drifter_pair = pd.DataFrame({'time': time,'delta_t': tau, 
                                 'L_tau': L_tau, 'delta_r': delta_r,
                                 'dr2': dr2}, 
                            columns=['time', 'delta_t', 'L_tau', 'delta_r','dr2']) 
    drifter_pair.set_index('time',inplace=True)
    return drifter_pair

def rel_disp(dr2_list, time, nob=nob, freq=freq):
    '''
    Compute the ensemble mean of the dr2 - mean squared change of separation.
    Parameters
    ----------
    dr2_list: list
        list that holds all the dr2 for each drifter pair.
    time: DatetimeIndex
        array of the hourly drifter dates
    nob:  int
        number of drifter pairs
    freq: int
        the frequency of window over which calulation is run
    Returns
    -------
    RD:
        Dataframe with the relative disperion and corresponding time
    '''
    day = [abs(time[d]-time[0])/timedelta(days=1) for d in range(time.size)]
    tau = [abs((time[l+freq]-time[l]))/timedelta(days=1) for l in range(len(time)-freq)]
    tau = np.append(np.zeros(freq), tau)
    dr2_mean = []
    n=0
    for n in range((len(time))):
        dr2_mean.append((1/nob)*(np.nansum([l[n] for l in dr2_list], axis=0)))
    n=n+1
    RD = pd.DataFrame({'time': time, 'delta_t': tau, 'day no.': day, 
                       'dr2_mean': dr2_mean}, 
                            columns=['time', 'delta_t','day no.', 'dr2_mean']) 
    RD.set_index('time',inplace=True)
    return RD  

def disp_rate(L_tau, time, freq=freq):
    '''
    Two-particle dispersion rate - for a moving delta_r
    Parameters
    ----------
    L_tau: numpy array
        array of distances between two drifters
    tau: numpy array
        array of number of days 
    Returns
    -------
    D: the dispersion rate 
    '''
    # tau (delta_t)
    tau = [abs(time[l+freq]-time[l])/timedelta(days=1) for l in range(len(time)-freq)]
    N = len(tau)
    if N==0:
        D = np.nan
    else: 
        D = np.zeros(N)
        for k in range(N):
            D[k] = (L_tau[k+freq] - L_tau[k])/(L_tau[0]*tau[k])
    # pad with freq len of zeros
    tau = np.append(np.zeros(freq), tau)
    D   = np.append(np.zeros(freq), D)
    disp   = pd.DataFrame({'time': time,'D_rate': D, 'delta_t': tau}, 
                            columns=['time','D_rate', 'delta_t']) 
    disp.set_index('time',inplace=True)
    return disp

def std_disp_rate(D_list, time, nob=nob, freq=freq):
    '''
    Standard deviation of Deformation D- 
    adapted Girard and others (2009)
    Computes the standard deviation of dispersion rate normalised by
    length and time:
        
    std_D = (< (D_rate - < D_rate >)^2 >)^0.5
    
    where D is the dispersion rate and the angular brackets denote the 
    ensemble mean caluclated over N pairs of buoys initally separated 
    by L0 +- L_tau and over a time interval tau.
    Parameters
    ----------
    D_list: list
        list that holds all the dr2 for each drifter pair.
    time: DatetimeIndex
        array of the hourly drifter dates
    nob:  int
        number of drifter pairs
    freq: int
        the frequency of window over which calulation is run
    Returns
    -------
    std_Dr: Dataframe containing the std of dispersion rate
    '''
    day = [abs(time[d]-time[0])/timedelta(days=1) for d in range(time.size)]
    tau = [abs(time[l+freq]-time[l])/timedelta(days=1) for l in range(len(time)-freq)]
    # pad with freq len of zeros
    tau = np.append(np.zeros(freq), tau)
    sD     = []
    m=0
    for m in range(len(tau)):
        sD.append(np.sqrt((1/nob)*np.nansum(np.square([l[m] for l in D_list] 
            - ((1/nob)*(np.nansum([l[m] for l in D_list], axis=0)))), axis=0)))
        m=m+1
    std_Dr = pd.DataFrame({'time': time, 'delta_t': tau, 'day no.':day,'std_D': sD}, 
                            columns=['time', 'delta_t','day no.', 'std_D']) 
    std_Dr.set_index('time',inplace=True)
    return std_Dr

################################
#%% Relative dispersion variance
# methods by Lukovich and others
# (2017, 2021)
################################
    
# 1. read in the separation between buoys from above
# 2. create a list of the zonal and meridional separations
# 3. calculate the ensemble mean of the zonal, total and meridional R2
    
def rel_disp_var(Lx_list, Ly_list, time, nob=nob):
    '''
    Two-particle dispersion - 
    adapted from Lukovich et al. (2017, 2021)
    Computes the variance of two-particle dispersion.
    The method should be computed on an ensemble of floats
    
    R^2 = < (xi(t) - xi+1(t) - <xi(t) - xi+1(t)>)**2 >  
    
    where xi is the [x or y] position of drifter i at time t, and the angular 
    brackets denote the ensemble mean over all particles. 
    For a single particle, the outer ensemble mean is dropped, and the 
    centroid is the mean location of the particle computed over time t-t0 
    Parameters
    ----------
    Lx_list: list
        Contains the zonal separation
    Ly_list: list
        Contains the meridional separation
    time: DatetimeIndex
        array of the hourly drifter dates
    nob:  int
        number of buoy (drifter) pairs 
    Returns
    -------
    R2: Dataframe containing the absolute dispersion for the buoy cluster
    '''
    # time interval
    length = math.ceil(time.size)
    # tau -> change in time from t(0)
    tau    = [abs(time[d]-time[0])/timedelta(days=1) for d in range(length)]
    R2_z   = []
    R2_m   = []
    R2_t   = []
    n=0
    for n in range(length):
        # zonal(x) direction dispersion
        R2_z.append((1/nob)*np.nansum(np.square([l[n] for l in Lx_list] 
            - ((1/nob)*(np.nansum([l[n] for l in Lx_list], axis=0)))), axis=0))
        # meridional (y) direction disperion
        R2_m.append((1/nob)*np.nansum(np.square([l[n] for l in Ly_list] 
            - ((1/nob)*(np.nansum([l[n] for l in Ly_list], axis=0)))), axis=0))
        # total two-partcile dispersion (equal to the components combined)
        R2_t.append(R2_z[n] + R2_m[n])
        n=n+1
    R2 = pd.DataFrame({'time': time, 'delta_t': tau, 'R2_z': R2_z, 
                       'R2_m': R2_m, 'R2_t': R2_t}, 
                            columns=['time', 'delta_t', 'R2_z', 
                                     'R2_m', 'R2_t']) 
    R2.set_index('time',inplace=True)
    return R2

# end of code
