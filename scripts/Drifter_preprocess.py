#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Drifter_preprocess

This script prepares the drifter data files from raw data
It computes the geodetic drift from GPS location and standardizes variable 
names

Created on Tue Dec 28 15:03:27 2021

@author: vichi
"""

import numpy as np
import pandas as pd
from geopy.distance import geodesic

# date format = '%Y-%m-%d %H:%M:%S'
fmt = '%Y-%m-%d %H:%M:%S' # exclude %z
ODIR = '../data/'
BDIR='../data/raw/'

#%% Compute drift and speed
def drift_speed(drifter):
    lats = drifter['latitude (deg)'].values
    lons = drifter['longitude (deg)'].values
    
    DX = np.zeros(lats.size)*np.nan
    DY = np.zeros(lats.size)*np.nan
    for i in range(len(DX)-1):
        DX[i+1]=np.sign(lons[i+1]-lons[i])*geodesic((lats[i],lons[i]),
                                                    (lats[i],lons[i+1])).meters
        DY[i+1]=np.sign(lats[i+1]-lats[i])*geodesic((lats[i],lons[i]),
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
    
#%% Trident 1-3
buoy = ['Trident_U1','Trident_U2','Trident_U3']
file = ['Unit1.xlsx','Unit2.xlsx','Unit3.xlsx']
latname='Latitude'
lonname='Longitude'
timename='GPS_Time'
datename='Reading_Date'
for b,f in zip(buoy,file):
    drifter = pd.read_excel(BDIR+f)
    drifter['time'] = drifter[datename] + pd.to_timedelta(drifter[timename].astype(str))
    drifter.set_index('time',inplace=True)
    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
    drifter.drop(timename, axis=1, inplace=True)
    drifter.drop(datename, axis=1, inplace=True)
    drifter['Unit_ID'] = b
    drifter['Temp'] -= 40. # add the offset
    drifter.rename(columns={'Temp':'temperature (degC)'},inplace=True)
    drifter = drift_speed(drifter)
    drifter.to_csv(ODIR+b+'.csv')

#%% Trident 4
BFILE='Unit4.xlsx'
latname='Latitude'
lonname='Longitude'
timename='DateTime'
ID = 'Trident_U4'
drifter = pd.read_excel(BDIR+BFILE)
drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
drifter.set_index('time',inplace=True)
drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
drifter.drop(timename, axis=1, inplace=True)
drifter['Unit_ID'] = ID
drifter.rename(columns={'Temperature':'temperature (degC)'},inplace=True)
drifter = drift_speed(drifter)
drifter.to_csv(ODIR+ID+'.csv')

#%% AWI 
buoy = ['2014S9','2015P21','2015P6','2016P18','2016P28','2019P93']
file = ['2014S9_300234060376490_proc.csv',
        '2015P21_300234062886910_proc.csv',
        '2015P6_300234061771590_proc.csv',
        '2016P18_300234062881940_proc.csv',
        '2016P28_300234062787480_proc.csv',
        '2019P93_300234067701680_proc.csv'
        ]
for b,f in zip(buoy,file):   
    drifter = pd.read_csv(BDIR+f,parse_dates=True)
    drifter['time'] = pd.to_datetime(drifter['time'],format=fmt)
    drifter.set_index('time',inplace=True)
    drifter['Unit_ID'] = b
    drifter = drift_speed(drifter)
    drifter.to_csv(ODIR+b+'.csv')

#%% WIIOS
latname='gps/lat'
lonname='gps/lon'
timename='date/iso'

buoy = ['WIIOS_NYU1','WIIOS_NYU2']
file = ['10501_all_processed.csv',
        '10502_all_processed.csv']

for b,f in zip(buoy,file):
    drifter = pd.read_csv(BDIR+f,parse_dates=True)
    drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
    drifter.set_index('time',inplace=True)
    drifter.drop(timename, axis=1, inplace=True)
    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
    drifter.rename(columns={'temperature/air':'temperature (degC)'},inplace=True)
    drifter['Unit_ID'] = b
    drifter = drift_speed(drifter)
    drifter.to_csv(ODIR+b+'.csv')
#%% SAWS ISVP
latname=' LATITUDE'
lonname=' LONGITUDE'
timename='Data Date (UTC)'

buoy = ['ISVP1','ISVP2','ISVP3']
file = ['300234067003010-300234067003010-20191015T064320UTC.csv',
        '300234066992870-300234066992870-20191015T064314UTC.csv',
        '300234067002060-300234067002060-20191015T064316UTC.csv']

for b,f in zip(buoy,file):
    drifter = pd.read_csv(BDIR+f,parse_dates=True)
    drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
    drifter.set_index('time',inplace=True)
    drifter.drop(timename, axis=1, inplace=True)
    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
    drifter.rename(columns={' SST':'temperature (degC)'},inplace=True)
    drifter.rename(columns={' BP':'barometric_pressure (hPa)'},inplace=True)
    drifter['Unit_ID'] = b
    drifter = drift_speed(drifter.sort_index()) # NOTE index is reversed in raw
    drifter.to_csv(ODIR+b+'.csv')

#%% SAWS ISVP
latname='latitude'
lonname='longitude'
timename='time'

buoy = ['ISVP4','ISVP5','ISVP6']
file = ['300234066433050.xlsx',
        '300234066433051.xlsx',
        '300234066433052.xlsx']

for b,f in zip(buoy,file):
    drifter = pd.read_excel(BDIR+f,parse_dates=True)
    drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
    drifter.set_index('time',inplace=True)
    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
    drifter.rename(columns={'sst':'temperature (degC)'},inplace=True)
    drifter.rename(columns={'slp':'barometric_pressure (hPa)'},inplace=True)        
    drifter['Unit_ID'] = b
    drifter = drift_speed(drifter)
    drifter.to_csv(ODIR+b+'.csv')