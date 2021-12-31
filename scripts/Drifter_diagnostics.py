#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 11:23:04 2021

Reads preprocessed files

Overview of the main diagnostics of sea-ice drift
      - Drifter Velocity
      - Meander Coefficient (MC)
          1. Cumulative Method
          2. Discrete Method
      - Absolute (Single-Particle) Dispersion
      
@author: Ashleigh Womack, Marcello Vichi
"""
#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from matplotlib.dates import DayLocator,HourLocator,DateFormatter
from geopy.distance import geodesic, distance
from datetime import timedelta, datetime, date
from operator import attrgetter
import math


#%%
def distance_from(lat0,lon0,lats,lons):
    '''
    Parameters
    ----------
    lat0 : TYPE
        DESCRIPTION.
    lon0 : TYPE
        DESCRIPTION.
    lats : TYPE
        DESCRIPTION.
    lons : TYPE
        DESCRIPTION.

    Returns 
    -------
    delta_x : the list of geodetic distances from a point in m

    '''
    from geopy.distance import distance
    delta_x = [distance([lat0, lon0], 
                        [lats[k], lons[k]]).m for k in range(len(lats))]
    return delta_x
#%% Meander coeficients 

def meander_coeff(drifter,cumulative=False):
    ''' 
    Meander coefficient adapted from Vihma et al. (1996) 
    mc(t) = I(t)/delta_x(t)
    where I is the total trajectory distance travelled by the drifter at time (t) 
    and delta_x is geodesic distance of drifter at time (t) from origin
    '''
    # compute trajectory (substitute nan with 0)
    I = drifter['drift (m)'].fillna(0).cumsum()
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

    '''
    mc = drifter.groupby(pd.Grouper(freq=freq)).apply(meander_coeff)
    return mc

#%% Absolute dispersion - adapted from Lukovich et al. (2017, 2021)
# A_squared = <|( (xi(t)-xi(0)) - <(xi(t)-xi(0))> )**2)|>  
# where xi is the lat/lon position of the drifter at some time (t),
# and the angular brackets denote the ensemble mean. 
# The ensemble mean is calculated as: 1/N*sum(...), and is N the number of samples. 
# To determine the dynamic regime of the drifter, calculate the scaling component (alpha)
# where: 
# alpha > 1 - superdiffusive regime 
# alpha = 1 - diffusive regime   
# alpha < 1 - subdiffusive ("trapping") regime    
# Within the superdiffusive category, αlpha = 2 corresponds to a ballistic regime indicative of
# advection, αlpha = 5/3 to an elliptic regime, and αlpha = 5/4 to a hyperbolic regime

def abs_dispersion(drifter):
    lats = drifter['latitude (deg)']
    lons = drifter['longitude (deg)']
    N = len(drifter)
    abs_disp_z = np.zeros(N)  # zonal absolute dispersion
    abs_disp_m = np.zeros(N)  # meridional absolute dispersion
    abs_disp_t = np.zeros(N)  # total absolute dispersion 
    l=0
    for l in range(len(drifter)):
        abs_disp_z[l] = (1/N)*(np.nansum((abs(distance([lats[0],
                     lons[l]],[lats[0],lons[0]]).km - (1/N)*(np.nansum(distance([lats[0],
                         lons[l]],[lats[0],lons[0]]).km))))**2))
        
        abs_disp_m[l] = (1/len(drifter))*(np.nansum((abs(distance([lats[l],
                     lons[0]],[lats[0],lons[0]]).km - (1/len(drifter))*(np.nansum(distance([lats[l],
                         lons[0]],[lats[0],lons[0]]).km))))**2))
        
        abs_disp_t[l] = (1/len(drifter))*(np.nansum((abs(distance([lats[l],
                     lons[l]],[lats[0],lons[0]]).km - (1/len(drifter))*(np.nansum(distance([lats[l],
                         lons[l]],[lats[0],lons[0]]).km))))**2))
        l=l+1
    
        
    # Create sequence of the size of drifter    
    seq = np.arange(drifter.index.size)*4                   # multiply by sampling frequency (change this)          
    alpha, intercept = np.polyfit(abs_disp_t, seq, 1)       # the slope is the alpha
    return alpha



#%% Read in the buoy drifter data 
BDIR = '../data/' # directory of the drifter
theBuoy = 'Trident_U1'                                     
drifter = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
#drifter = drifter.loc[drifter.index[:]]   # locate any special date range             

#%% Compute the meander coefficients
mc_cum = meander_coeff(drifter,cumulative=True)
mc_1D = meander_coeff_discrete(drifter)
mc_5D = meander_coeff_discrete(drifter,freq='5D')


#%% start plotting
im=plt.scatter(lons,lats,c=drifter.index)
cbar=plt.colorbar(im)
cbar.ax.set_yticklabels(pd.to_datetime(cbar.get_ticks()).strftime(date_format='%Y-%m-%d'))

#%% Plot the Velocity, Cumulative MC and Discrete MC
mpl.rcParams['font.size'] = 25
PGrid = (3,2)
f = plt.figure(figsize=(40,30))

# Panel A: Drifter Velocity
ax_Vel = plt.subplot2grid(PGrid, (0,0))    
ax_Vel.plot(drifter.index, drifter['u'], color="royalblue")
ax_Vel.set_ylim(-0.75,0.75)                             # set the velocity limits
ax_Vel.set_ylabel("uᵢ (m s$^{-1}$)", color="royalblue")
ax_Vel.tick_params(axis='y', colors="royalblue")
ax_Vel.text(-0.15, 0.95,'(a)', transform=ax_Vel.transAxes,
      fontsize=25, va='top')
ax_Vel2=plt.twinx(ax_Vel)
ax_Vel2.plot(drifter.index, drifter['v'], color="orange")
ax_Vel2.set_xlabel('Date')
ax_Vel2.xaxis.set_major_formatter(DateFormatter('%m-%d'))
ax_Vel2.set_ylim(-0.75,0.75)                            # set the velocity limits
ax_Vel2.set_ylabel("vᵢ (m s$^{-1}$)", color="orange")
ax_Vel2.tick_params(axis='y', colors="orange")

# Panel B: Meander coefficient  
ax_MC = plt.subplot2grid(PGrid, (1,0))     
ax_MC.plot(drifter.index, MC, linewidth=2, color='red', label='Cumulative MC')
ax_MC.plot(drifter_D.index, MC_D, linewidth=2, color='green', label='Discrete MC')
ax_MC.xaxis.set_major_formatter(DateFormatter('%m-%d'))
ax_MC.set_xlabel('Date')
ax_MC.set_ylabel('Meander Coefficient (M)')
ax_MC.text(-0.15, -0.3,'(b)', transform=ax_Vel.transAxes,
      fontsize=25, va='top')
plt.legend()
  
# Panel C: Absolute (single-particle) dispersion - semilog plot 
ax_AD = plt.subplot2grid(PGrid, (2,0))  
ax_AD.plot(drifter.index, abs_disp_z, color="darkviolet", linewidth=2.5, label='Zonal')
ax_AD.plot(drifter.index, abs_disp_m, color="magenta", linewidth=2.5, label='Meridional')
ax_AD.plot(drifter.index, abs_disp_t, color="darkturquoise", linewidth=2.5, label='Total')
ax_AD.set_xlabel('Date')
ax_AD.xaxis.set_major_formatter(DateFormatter('%m-%d'))
ax_AD.set_yscale("log")
ax_AD.set_ylabel('$A^{2}_{zonal}$, ''$A^{2}_{meridional}$, ''$A^{2}_{total}$')  
ax_MC.text(-0.15, -1.5,'(c)', transform=ax_Vel.transAxes,
      fontsize=25, va='top')
plt.legend(loc=4,fontsize=25)

# end of code