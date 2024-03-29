# -*- coding: utf-8 -*-
"""
Drifter library

@author: Marcello Vichi, Ashleigh Womack (UCT)
"""
import numpy as np
import pandas as pd
from geopy.distance import distance

##################
#%% Drift analysis
##################

def distance_from(lat0,lon0,lats,lons):
    '''
    Compute the geodetic distance from a reference point in meters
    
    Parameters
    ----------
    lat0 : float
        Latitude of the origin
    lon0 : float
        Longitude of the origin
    lats : float, list or npdarray
        Latitudinal location
    lons : float, list or npdarray
        Longitudinal location

    Returns 
    -------
    delta_x : the list of geodetic distances from a point in m

    '''
    delta_x = [distance([lat0, lon0], 
                        [lats[k], lons[k]]).m for k in range(len(lats))]
    return np.array(delta_x)

def coordinates(drifter):
    '''
    Compute coordinates in km from the deployment
    Parameters
    ----------
    drifter : DataFrame
        DESCRIPTION.

    Returns
    -------
    drifter : DataFrame
        The updated DataFrame containing the coordinates
    '''
    lats = drifter['latitude (deg)']
    lons = drifter['longitude (deg)']
    N = len(lats)    
    dist_x = np.array([distance([lats[0], lons[0]], 
                    [lats[0], lons[k]]).km for k in range(N)])
    dist_y = np.array([distance([lats[0], lons[0]], 
                    [lats[k], lons[0]]).km for k in range(N)])
    drifter['x (km)'] = dist_x
    drifter['y (km)'] = dist_y
    drifter['distance (km)'] = np.sqrt(dist_x**2+dist_y**2)
    return drifter

def drift_speed(drifter):
    '''
    Compute drift and speed
    Parameters
    ----------
    drifter : TYPE
        DESCRIPTION.

    Returns
    -------
    drifter : TYPE
        DESCRIPTION.

    '''
    lats = drifter['latitude (deg)'].values
    lons = drifter['longitude (deg)'].values    
    DX = np.zeros(lats.size)
    DY = np.zeros(lats.size)
    for i in range(len(DX)-1):
        DX[i+1]=np.sign(lons[i+1]-lons[i])*distance((lats[i],lons[i]),
                                                    (lats[i],lons[i+1])).meters
        DY[i+1]=np.sign(lats[i+1]-lats[i])*distance((lats[i],lons[i]),
                                                    (lats[i+1],lons[i])).meters
    drifter['drift_x (m)'] = DX
    drifter['drift_y (m)'] = DY
    drifter['drift (m)'] = np.sqrt(DX**2+DY**2)  
    time = drifter.index
    DT = []
    for i in range(len(time)-1):
        delta_t = abs(time[i]-time[i+1]) # use timedelta
        DT.append(delta_t.seconds)
    DT = np.append(0., DT)
    drifter['delta_t (s)'] = DT    
    # Modular speed of the drifter
    drifter['U (m/s)'] = drifter['drift (m)']/DT
    # Zonal and meridional components
    drifter['u (m/s)'] = DX/DT
    drifter['v (m/s)'] = DY/DT
    return drifter

def drifter_check(drifter,criterion='max',quant=0.99, maxspeed=1.5):
    '''
    Quality check:
        1. speed: retain data <= 99th percentile
        
    Parameters
    ----------
    drifter : TYPE
        DESCRIPTION.
    criterion : string
        the type of criterion to apply
        'max'
        'quantile'

    Returns
    -------
    drifter

    '''
    drifter.replace(np.inf,np.nan) # remove all inf values and replace with nan
    speed = drifter['U (m/s)']
    if (criterion=='max'):
        out = (speed > maxspeed)
    else:
        out = (speed > speed.quantile(quant))
    index_out = speed[out].index
    drifter.drop(index_out,inplace=True)
    return drifter


def meander_coeff(drifter,cumulative=False):
    ''' 
    Meander coefficient adapted from Vihma et al. (1996) 
    mc(t) = I(t)/delta_x(t)
    where I is the total trajectory distance travelled by the drifter at time (t) 
    and delta_x is geodesic distance of drifter at time (t) from origin
    '''
    # distance along trajectory
    drift = drifter['drift (m)']
    try:
        I = drift.cumsum() #-drift[0] # in case the trajectory does not start from 0,0
    except IndexError:
        mc_ret = np.nan
    else:
        lats = drifter['latitude (deg)']
        lons = drifter['longitude (deg)']
        delta_x = distance_from(lats[0], lons[0], lats, lons)
        mc = I/delta_x
        if (cumulative):
            mc_ret = mc
        else:
            mc_ret = mc[-1]
    return mc_ret

def meander_coeff_discrete(drifter,freq='D',label='left'):
    '''
    Parameters
    ----------
    drifter : pandas DataFrame
        The drifter DataFrame
    freq : string, optional
        The time period over which to compute the coefficient. 
        The default is 'D'.
    label : the date label 'left' or 'right' for the time interval

    Returns
    -------
    mc: pandas Series contaning the meander coeffient at frequency freq
        The index is located on the last day of the time interval

    '''
    mc = drifter.groupby(pd.Grouper(freq=freq,
                                    label=label)).apply(meander_coeff)
    return mc


#%% 
def abs_dispersion(drifter,cumulative=False,which='x'):
    '''
    Single particle dispersion - 
    adapted from Taylor (1921) and Lukovich et al. (2017, 2021)
    Computes the variance of displacement around the cluster centroid, 
    an indicator of dispersion.
    The method should be computed on an ensemble of floats
    
    A^2 = < (xi(t) - <xi(t)-xi(0)>)**2 >  
    
    where xi is the [x or y] position of drifter i at time t, and the angular 
    brackets denote the ensemble mean over all particles. 
    For a single particle, the outer ensemble mean is dropped, and the 
    centroid is the mean location of the particle computed over time t-t0 

    Parameters
    ----------
    drifter : pandas DataFrame
        Drifter DataFrame
    cumulative : boolean
        Compute cumulative (as a function of time) or value over the period
    which : string
        'x' or 'y' for zonal or meridional variance

    Returns
    -------
    var_x: the variance along direction x [or y] 

    '''
    x = drifter[which+' (km)'].values
    N = len(x)
    if N==0:
        var_x = np.nan
    else:
        var_x = np.zeros(N)
        for l in range(N):
            var_x[l] = np.square(x[l] - np.mean(x[:l]-x[0]))
        if (not(cumulative)):
            var_x = var_x[-1]
    return var_x

def abs_dispersion_discrete(drifter,freq='D',which='x'):
    '''
    Parameters
    ----------
    drifter : pandas DataFrame
        The drifter DataFrame
    freq : string, optional
        The time period over which to compute dispersion. 
        The default is 'D'.

    Returns
    -------
    A2: pandas Series containing the variance at frequency freq

    '''
    A2 = drifter.groupby(pd.Grouper(freq=freq,
                                    label='right')).apply(abs_dispersion,which=which)
    return A2

# end of code
