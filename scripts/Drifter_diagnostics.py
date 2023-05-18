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
          3. Daily-averaged Discrete Method
      - Absolute (Single-Particle) Dispersion
      
@author: Ashleigh Womack, Marcello Vichi
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from Drifter import drifter_check
from Drifter import meander_coeff,meander_coeff_discrete
from Drifter import abs_dispersion,abs_dispersion_discrete
import matplotlib.dates as mdates
import seaborn as sns

BDIR = '../data/' # directory of the drifter data


#%% USER INPUT
# Read in the buoy drifter data 

theBuoy = 'ISVP1'                                     
drifter = pd.read_csv(BDIR+theBuoy+'.csv',index_col='time',parse_dates=True)
# dates ticks to be used in the final figure
DAY = 10

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

#%% extract coordinates and compute diagnostics for plotting
lats = drifter['latitude (deg)']
lons = drifter['longitude (deg)']
x = drifter['x (km)']
y = drifter['y (km)']
dfo = drifter['distance (km)'] # distance from origin
start = drifter.index.min()
end = drifter.index.max()

#%% Compute the meander coefficients
mc_cum = meander_coeff(drifter, cumulative=True)
mc_1D = meander_coeff_discrete(drifter)
mc_5D = meander_coeff_discrete(drifter,freq='5D',label='right')

f,ax = plt.subplots()
mc_cum.plot(ax=ax,label='cumulative')
mc_5D.plot(marker='o',ax=ax,label='5D')
mc_1D.plot(marker='o',ax=ax,label='1D')
ax.set_ylabel('Meander coefficient')
plt.legend()

#%% Compute cumulative dispersion - this is for 1 single buoy
abs_disp_x = abs_dispersion(drifter,cumulative=True,which='x')
abs_disp_y = abs_dispersion(drifter,cumulative=True,which='y')
drifter['disp'] = abs_disp_x+abs_disp_y
drifter['disp_x'] = abs_disp_x
drifter['disp_y'] = abs_disp_y


# Compute discrete dispersion
A2_1D_x = abs_dispersion_discrete(drifter,freq='1D',which='x')
A2_1D_y = abs_dispersion_discrete(drifter,freq='1D',which='y')
A2_5D_x = abs_dispersion_discrete(drifter,freq='5D',which='x')
A2_5D_y = abs_dispersion_discrete(drifter,freq='5D',which='y')

# Compute numerical acceleration (first-order forward)
speed = drifter['U (m/s)']
acc = speed.diff()/speed.index.to_series().diff().dt.total_seconds()
acc_1D = acc.rolling('1D').mean()*1.e5 # daily running mean in m/s2*10^5

#
# To determine the dynamic regime of the drifter, 
#     disp ~ T^alpha
# calculate the scaling component (alpha) with nonlinear curve fitting
# where: 
# alpha > 1 - superdiffusive regime 
# alpha = 1 - diffusive regime   
# alpha < 1 - subdiffusive ("trapping") regime    
# Within the superdiffusive category, αlpha = 2 corresponds to a ballistic regime indicative of
# advection, α = 5/3 to an elliptic regime, and αlpha = 5/4 to a hyperbolic regime  
ys = drifter['disp'].dropna()
xs = mdates.date2num(ys.index)
xs = xs - xs[0]
from scipy.optimize import curve_fit
def f(x, a):
    return np.power(x,a)
alpha, pcov = curve_fit(f, xs, ys.values)
# plot the result to check
plt.figure()
plt.plot(xs,ys.values)
plt.plot(xs,f(xs,alpha))


#%% Plot the Velocity, Cumulative MC and Discrete MC

def set_daterange(ax,byday=15):
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=byday))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%d'))
    ax.set_xlim(start.date(),end.date())
    #ax_Vel.tick_params(axis='x', labelrotation = 45)


mpl.rcParams['font.size'] = 11
PGrid = (3,2)
f = plt.figure(figsize=(15,9))

# Panel A: Trajectory
ax_traj = plt.subplot2grid(PGrid, (0,0),colspan=2)  
plt.plot(x,y,'-k',linewidth=0.5)
im = plt.scatter(x,y,c=mdates.date2num(drifter.index),cmap='tab20b')
ax_traj.set_xlabel('km')
ax_traj.set_ylabel('km')
cbar = plt.colorbar(im)
ax_twin = plt.twiny(ax_traj)
ax_twin.plot(lons,y,linewidth=0)
ax_twin.set_xlabel('deg longitude')
loc = mdates.AutoDateLocator()
cbar.ax.yaxis.set_major_locator(loc)
cbar.ax.yaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
ax_traj.text(-0.07, 1,'(a)', transform=ax_traj.transAxes,
      fontsize=18, va='top')

# Panel B: Drifter Velocity
ax_Vel = plt.subplot2grid(PGrid, (1,0))    
ax_Vel.plot(drifter.index, drifter['u (m/s)'], color="royalblue", label='Zonal')
ax_Vel.set_ylim(-1.,1.)                             # set the velocity limits
ax_Vel.set_ylabel("Velocity (m s$^{-1}$)")
#ax_Vel.tick_params(axis='y', colors="royalblue")
# ax_Vel2=plt.twinx(ax_Vel)
# ax_Vel2.plot(drifter.index, drifter['v (m/s)'], color="orange")
# ax_Vel2.set_ylim(-0.75,0.75)                            # set the velocity limits
# ax_Vel2.set_ylabel("v (m s$^{-1}$)", color="orange")
# ax_Vel2.tick_params(axis='y', colors="orange")
ax_Vel.plot(drifter.index, drifter['v (m/s)'], color="orange",label='Meridional')
ax_Vel.grid(axis='y')
plt.legend()
set_daterange(ax_Vel,byday=DAY)
ax_Vel.text(-0.13, 1,'(b)', transform=ax_Vel.transAxes,
      fontsize=18, va='top')

# Panel C: acceleration (alternative)
# ax_AC = plt.subplot2grid(PGrid, (1,1))
# ax_AC.plot(drifter.index,acc_1D, linestyle='-', linewidth=2)
# ax_AC.grid(axis='y')
# ax_AC.set_ylabel('Acceleration (10$^5$ m s$^{-2}$)')  
# set_daterange(ax_AC)
# ax_AC.text(-0.13, 1,'(c)', transform=ax_AC.transAxes,
#       fontsize=18, va='top')
ax_dist = plt.subplot2grid(PGrid, (1,1))
sns.histplot(drifter['U (m/s)'],ax=ax_dist)
ax_dist.text(-0.13, 1,'(c)', transform=ax_dist.transAxes,
      fontsize=18, va='top')

# Panel C: Meander coefficient  
ax_MC = plt.subplot2grid(PGrid, (2,0))     
ax_MC.plot(mc_cum.index, mc_cum, linewidth=2, color='red', label='Cumulative')
ax_MC.plot(mc_5D.index, mc_5D, linewidth=2, color='green', 
           marker='|', label='Discrete (5D)')
ax_MC.plot(mc_1D.index, mc_1D, linewidth=2, color='lightblue', 
           marker='|', label='Discrete (1D)')
ax_MC.set_ylabel('Meander Coefficient')
set_daterange(ax_MC,byday=DAY)
ax_MC.text(-0.13, 1,'(d)', transform=ax_MC.transAxes,
      fontsize=18, va='top')
plt.legend()
  
# Panel D: single-particle dispersion - semilog plot 
ax_AD = plt.subplot2grid(PGrid, (2,1))  
ax_AD.plot(drifter.index, drifter['disp_x'], color="darkviolet", linewidth=2, label='Zonal')
ax_AD.plot(drifter.index, drifter['disp_y'], color="magenta", linewidth=2, label='Meridional')
ax_AD.plot(drifter.index, drifter['disp'], color="darkturquoise", linewidth=2, label='Total')
ax_AD.set_yscale("log")
ax_AD.set_ylabel('Variance (km$^2$)')  
set_daterange(ax_AD,byday=DAY)
ax_AD.text(-0.13, 1,'(e)', transform=ax_AD.transAxes,
      fontsize=18, va='top')
plt.legend(loc=4)

plt.tight_layout()

#%%
# Plot trajectory in km from deployment
# Alternative method
# ticks = pd.date_range(drifter.index.min(),
#                      drifter.index.max(),freq='5D',normalize=True)
# ticks = [date.value for date in ticks]
# f,ax = plt.subplots()
# im = plt.scatter(x,y,c=drifter.index)
# ax.set_xlabel('km')
# ax.set_ylabel('km')
# cbar = plt.colorbar(im,ticks=ticks)
# cbar.ax.set_yticklabels(pd.to_datetime(cbar.get_ticks()).strftime(date_format='%Y-%m-%d'))


# end of code


