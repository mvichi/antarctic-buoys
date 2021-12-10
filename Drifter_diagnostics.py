#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 11:23:04 2021

@author: Ashleigh Womack

Overview of the main diagnostics of sea-ice drift
      - Drifter Velocity
      - Meander Coefficient (MC)
          1. Cumulative Method
          2. Discrete Method
      - Absolute (Single-Particle) Dispersion
"""
#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from matplotlib.dates import DayLocator,HourLocator,DateFormatter
from geopy.distance import geodesic, vincenty, distance
from datetime import timedelta, datetime, date
from operator import attrgetter
from numpy import meshgrid
from mpl_toolkits import basemap 
import math

#%% Read in the buoy drifter data 
BDIR = '/home/.../'                                     # directory of the drifter
drifter = =pd.read_excel(ISVP_dir+'drifter.xlsx',       # change name of file
                   index_col=4,parse_dates=True,dayfirst=True) # index column with the date + time

# If the date and time are separated in two different columns: 
# Load the drifter date + time coloumns and make them the index 
#BDIR='/home/..../'
#drifter = pd.read_excel(BDIR+'drifter.xlsx')
#drifter['DateTime'] = drifter['Reading_Date'] + pd.to_timedelta(drifter['GPS_Time'].astype(str))
#drifter.set_index('DateTime',inplace=True)

drifter = drifter.loc[drifter.index[:]]                 # locate any special date range             

#%% Compute the speed of the drifter
lats = drifter['latitude'].values                       # change the name of latitude column name
lons = drifter['longitude'].values                      # change the name of longitude column name

DX = np.zeros(lats.size)*np.nan
DY = np.zeros(lats.size)*np.nan
for i in range(len(DX)-1):
    DX[i+1]=np.sign(lons[i+1]-lons[i])*geodesic((lats[i],lons[i]),(lats[i],lons[i+1])).meters
    DY[i+1]=np.sign(lats[i+1]-lats[i])*geodesic((lats[i],lons[i]),(lats[i+1],lons[i])).meters
drifter['DX']=DX
drifter['DY']=DY

# If the format of date has +00 at the end change format to: fmt = '%Y-%m-%d %H:%M:%S%z'   
fmt = '%Y-%m-%d %H:%M:%S' #%z'                             
time = drifter.index[:].astype(str)
DT = []
for i in range(len(time)-1):
    delta_t = abs(datetime.strptime(time[i], fmt)-datetime.strptime(time[i+1], fmt))
    DT.append(delta_t.total_seconds())
DT = np.append(0, DT)

# Modular speed of the drifter 
drifter['U'] = np.sqrt(DX**2+DY**2)/(DT)
# Zonal and meridional components 
drifter['u']=DX/(DT)
drifter['v']=DY/(DT)

#%% Meander coeficient - intergrated response method adapted from Vihma et al. (1996). 
# MC = I/delta_x
# where I is the total trajectory distance travelled by the drifter at time (t)
# and delta_x is the total displacement of the drifter at time (t)
from geopy.distance import geodesic 
from numpy import cumsum

I = []
for j in range(len(drifter)-1):
    dd = distance([lats[j], lons[j]], [lats[j+1], lons[j+1]]).km
    I.append(dd)
drifter['I'] = np.append(0,cumsum(I))

delta_x = []
for k in range(len(drifter)):
    dx = distance([lats[0], lons[0]], [lats[k], lons[k]]).km
    delta_x.append(dx)
drifter['delta_x'] = delta_x

MC = drifter['I']/drifter['delta_x']
drifter['Cumulative MC'] = MC

#%% Meander coefficient - discrete method
# resample to a time interval (e.g. daily, 5-daily, monthly)
drifter_D = drifter.resample('D').mean()                    # resampled to daily
# If there is a of missing data when resampling fill the data
drifter_D = drifter_D.fillna(method="ffill")

lats_D = drifter_D['Latitude.1'].values                    # change the name of latitude column name
lons_D = drifter_D['Longitude.1'].values                   # change the name of longitude column name

# Total trajectory distance travelled by the drifter each day
I_D = []
for l in range(len(drifter_D)-1):
    dd_D = drifter_D.I[l+1]-drifter_D.I[l]
    I_D.append(dd_D)
drifter_D['I_daily'] = np.append(0,I_D)

# total displacement travelled by the drifter each day
delta_x_D = []
for m in range(len(drifter_D)-1):
    dx_D = distance([lats_D[m], lons_D[m]], [lats_D[m+1], lons_D[m+1]]).km
    delta_x_D.append(dx_D)
drifter_D['delta_x_daily'] = np.append(0,delta_x_D)

MC_D = drifter_D['I_daily']/drifter_D['delta_x_daily'] 
drifter_D['Discrete MC'] = MC_D


#%% Absolute dispersion - adapted from Lukovich et al. (2017, 2021)
# A_squared = <|( (xi(t)-xi(0)) - <(xi(t)-xi(0))> )**2)|>  
# where xi is the lat/lon position of the drifter at some time (t),
# and the angular brackets denote the ensemble mean. 
# The ensemble mean is calculated as: 1/N*sum(...), and is N the number of samples. 
abs_disp_z = np.zeros(len(drifter))                      # zonal absolute dispersion
abs_disp_m = np.zeros(len(drifter))                      # meridional absolute dispersion
abs_disp_t = np.zeros(len(drifter))                      # total absolute dispersion 
l=0
for l in range(len(drifter)):
    abs_disp_z[l] = (1/len(drifter))*(np.nansum((abs(distance([lats[0],
                 lons[l]],[lats[0],lons[0]]).km - (1/len(drifter))*(np.nansum(distance([lats[0],
                     lons[l]],[lats[0],lons[0]]).km))))**2))
    
    abs_disp_m[l] = (1/len(drifter))*(np.nansum((abs(distance([lats[l],
                 lons[0]],[lats[0],lons[0]]).km - (1/len(drifter))*(np.nansum(distance([lats[l],
                     lons[0]],[lats[0],lons[0]]).km))))**2))
    
    abs_disp_t[l] = (1/len(drifter))*(np.nansum((abs(distance([lats[l],
                 lons[l]],[lats[0],lons[0]]).km - (1/len(drifter))*(np.nansum(distance([lats[l],
                     lons[l]],[lats[0],lons[0]]).km))))**2))
    l=l+1

# To determine the dynamic regime of the drifter, calculate the scaling component (alpha)
# where: 
# alpha > 1 - superdiffusive regime 
# alpha = 1 - diffusive regime   
# alpha < 1 - subdiffusive ("trapping") regime    
# Within the superdiffusive category, αlpha = 2 corresponds to a ballistic regime indicative of
# advection, αlpha = 5/3 to an elliptic regime, and αlpha = 5/4 to a hyperbolic regime
    
# Create sequence of the size of drifter    
seq = np.arange(drifter.index.size)*4                   # multiply by sampling frequency (change this)          
slope, intercept = np.polyfit(abs_disp_t, seq, 1)       # the slope is the alpha

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