# -*- coding: utf-8 -*-
"""
Drifter library

@author: Marcello Vichi, Ashleigh Womack
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
        I = drift.cumsum()-drift[0] # in case the trajectory does not start from 0,0
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

def meander_coeff_discrete(drifter,freq='D'):
    '''
    Parameters
    ----------
    drifter : pandas DataFrame
        The drifter DataFrame
    freq : string, optional
        The time period over which to compute the coefficient. 
        The default is 'D'.

    Returns
    -------
    mc: pandas Series contaning the meander coeffient at frequency freq
        The index is located on the last day of the time interval

    '''
    mc = drifter.groupby(pd.Grouper(freq=freq,
                                    label='right')).apply(meander_coeff)
    return mc

def meander_coeff_mean(drifter, n=24):
    ''' 
    Meander coefficient adapted from Vihma et al. (1996) 
    mc(t) = I(t)/delta_x(t)
    where I is the total trajectory distance travelled by the drifter at time (t) 
    and delta_x is geodesic distance of drifter at time (t) from origin.
    
    This method calculates the mc for each day for every index and then averages
    over 24 hours. 
    
    Parameters
    ----------
    drifter : pandas DataFrame
        The drifter DataFrame
    n: number of indices in 1 day (default is 24 - 24 hours per day)
    Returns
    -------
    mc_mean: pandas Dataframe contaning the average meander coeffient for each day
                The index is located on the first day (left) of the time interval.
    '''
    lats = drifter['latitude (deg)']
    lons = drifter['longitude (deg)']
    # displacement between reference point and 24 hours later
    delta_x = [distance([lats[k], lons[k]], 
                        [lats[k+n], lons[k+n]]).km for k in range(len(lats)-n)]
    # distance along trajectory
    DX = np.zeros(lats.size)
    DY = np.zeros(lats.size)
    for i in range(len(DX)-1):
        DX[i+1]=distance((lats[i],lons[i]),(lats[i],lons[i+1])).kilometers
        DY[i+1]=distance((lats[i],lons[i]),(lats[i+1],lons[i])).kilometers
    drifter['dist_x (km)'] = DX
    drifter['dist_y (km)'] = DY
    drifter['dist (km)'] = np.sqrt(DX**2+DY**2)
    dist = drifter['dist (km)']
    # trajectory distance
    # this cumulatively sums the drift distance within 24 hours (last column)
    I = pd.DataFrame([dist[i:i+n].cumsum().tolist() for i in range(len(dist[:-n]))])[n-1] 
    mc = pd.DataFrame(I/delta_x)
    # add the drifters time to the dataframe
    mc['time'] = drifter.index[0+n:]
    mc.set_index('time',inplace=True)
    # resample to daily time interval
    mc_1D_mean = mc.resample('D').mean()
    return mc_1D_mean

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

########################
#%% Atmospheric response
########################

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
    # add to variables to drifter DataFrame
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
    
    q = ax.quiver(date2num(time), [[0]*len(time)], u, v,
                  angles='uv', width=width, headwidth=headwidth,
                  headlength=headlength, headaxislength=headaxislength,
                  **kw)

    ax.axes.get_yaxis().set_visible(False)
    ax.xaxis_date()
    return q


#%% Compute the wind factor, turning angle and vector coefficient of determination
import math    
def drift_response(drifter):
    '''
    Relationship between the ice motions and ERA5 winds - 
    adapted from Kimura and Wakatsuchi (2000), Kimura (2004), 
    Fukamachi et al. (2011)
    Computes the wind factor, turning angle and vector coefficient 
    of determination between the wind speed and ice (drifter) speed.
    We also compute the Pearon coefficient of determination.
    Parameters
    ----------
    drifter : pandas DataFrame
        Drifter DataFrame
    Returns
    -------
    theta_degs, WF, R2v, R2p : numpy.float64
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
    theta_degs = math.degrees(theta) # in degrees (neg due to left delfection)
    
    # wind factor     
    c1 = (np.cos(theta))*(u10u)
    c2 = (np.sin(theta))*(v10u)
    c3 = (np.sin(theta))*(u10v)
    c4 = (np.cos(theta))*(v10v)

    WF = ((c1 + c2 - c3 + c4)/(u10_2 + v10_2))  
    
    # vector coefficient of determination
    R2v = np.square((c1 + c2 - c3 + c4)/(np.sqrt(u10_2 + v10_2)*np.sqrt(u_2 + v_2)))

    # Pearson coefficient of determination
    R2p = np.square(scipy.stats.pearsonr(drifter['U10_E5 (m/s)'][1:], 
                                                 drifter['U (m/s)'][1:]))
    
    return theta, theta_degs, WF, R2v, R2p

# Use the drift repsonse measurements to estimate the underlying currents
def current_estimation(drifter,theta, WF):
    '''
    Estimate the underlying current velcotiy components -
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
        Mean ocean current components
    C : numpy.float64
        Mean modular current speed
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
    # multiply by 100 to get into cm/s
    Cuv = (M1 - WF*(M2.dot(M3)))*100
    C     = np.sqrt(Cuv[0]**2 + Cuv[1]**2)
    
    return Cuv, C

# end of code
