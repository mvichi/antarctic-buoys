#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Drifter_preprocess

This script prepares the drifter data files from raw data
It computes the geodetic drift from GPS location and standardizes variable 
names

Created on Tue Dec 28 15:03:27 2021

@author: Marcello Vichi (UCT)
"""

import pandas as pd
from Drifter import drift_speed, coordinates
import scipy.io as sio
from scipy import signal, stats

# date format = '%Y-%m-%d %H:%M:%S'
fmt = '%Y-%m-%d %H:%M:%S' # exclude %z
ODIR = '../data/'
BDIR='../data/raw/'


    
#%% Trident 1-3
buoy = ['Trident_U1','Trident_U2','Trident_U3']
file = ['Unit1.xlsx','Unit2.xlsx','Unit3.xlsx']
latname ='Latitude'
lonname ='Longitude'
timename ='GPS_Time'
datename ='Reading_Date'
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
    drifter = coordinates(drifter)
    drifter.to_csv(ODIR+b+'.csv')

#%% Trident 4 - 2019 spring cruise
#   NOTE that this buoy was switched on ealier to check the its transmission
#   of data but it was only deployed on the 30th October 2019 
#  (Ryan-Keogh and Vichi, 2020). 
BFILE ='Unit4.xlsx'
latname ='Latitude.1'
lonname ='Longitude.1'
timename ='UTC Date Time'
ID = 'Trident_U4'
drifter = pd.read_excel(BDIR+BFILE)
drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]') # change the fmt of datetimeindex
drifter.set_index('time',inplace=True)
drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
drifter.drop(timename, axis=1, inplace=True)
drifter['Unit_ID'] = ID
drifter.rename(columns={'Temperature':'temperature (degC)'},inplace=True)
drifter = drift_speed(drifter)
drifter = coordinates(drifter)
drifter = drifter.resample('4H').mean() # contains missing timestamps
drifter = drifter = drifter.fillna(method='ffill') # interpolate missing data
drifter = drifter.loc['2019-10-30 12:00':]
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
    drifter = coordinates(drifter)
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
    drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]')  # change the fmt of datetimeindex
    drifter.set_index('time',inplace=True)
    drifter.drop(timename, axis=1, inplace=True)
    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
    drifter.rename(columns={'temperature/air':'air temperature (degC)'},inplace=True)
    drifter['Unit_ID'] = b
    drifter = drift_speed(drifter)
    drifter = coordinates(drifter)
    drifter.to_csv(ODIR+b+'.csv')
    
    
#%% SAWS ISVP (winter)
# These data contain missing timestamps.
# N.B. 
# The buoys have been renamed to match their sequence of deployment. 
# In the cruise report (Ryan-Keogh and Vichi, 2020) ISVP2 and ISVP3 are swapped
# ISVP2 was ID 300234066992870 and it is now 300234067002060
latname=' LATITUDE'
lonname=' LONGITUDE'
timename='Data Date (UTC)'

buoy = ['ISVP3','ISVP2']
file = ['300234066992870-300234066992870-20191015T064314UTC.csv',
        '300234067002060-300234067002060-20191015T064316UTC.csv']

for b,f in zip(buoy,file):
    drifter = pd.read_csv(BDIR+f,parse_dates=True)
    drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
    drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]') # change the fmt of datetimeindex
    drifter.set_index('time',inplace=True)
    drifter.drop(timename, axis=1, inplace=True)
    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
    drifter.rename(columns={' SST':'surf temperature (degC)'},inplace=True)
    drifter.rename(columns={' BP':'barometric_pressure (hPa)'},inplace=True)
    drifter.rename(columns={' BPT':'air temperature (degC)'},inplace=True)
    drifter['Unit_ID'] = b
    drifter = drifter.sort_index() # NOTE index is reversed in raw
    drifter = drift_speed(drifter)
    drifter = coordinates(drifter)
    drifter = drifter.resample('1H').mean()
    drifter = drifter.fillna(method='ffill') # interpolate missing data
    drifter.to_csv(ODIR+b+'.csv')
    
# process ISVP1 at 30 minutes frequency
b ='ISVP1_30min'
f ='300234067003010-300234067003010-20191015T064320UTC.csv'
drifter = pd.read_csv(BDIR+f,parse_dates=True)
drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]') # change the fmt of datetimeindex
drifter.set_index('time',inplace=True)
drifter.drop(timename, axis=1, inplace=True)
drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
drifter.rename(columns={' SST':'temperature (degC)'},inplace=True)
drifter.rename(columns={' BP':'barometric_pressure (hPa)'},inplace=True)
drifter.rename(columns={' BPT':'air temperature (degC)'},inplace=True)
drifter['Unit_ID'] = b
drifter = drifter.sort_index() # NOTE index is reversed in raw
drifter = drift_speed(drifter)
drifter = coordinates(drifter)
drifter = drifter.resample('30min').mean()
drifter = drifter.fillna(method='ffill') # interpolate missing data
drifter.to_csv(ODIR+b+'.csv')

# reprocess ISVP1 at 1 hour frequency
b='ISVP1'
f='300234067003010-300234067003010-20191015T064320UTC.csv'
drifter = pd.read_csv(BDIR+f,parse_dates=True)
drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]') # change the fmt of datetimeindex
drifter.set_index('time',inplace=True)
drifter.drop(timename, axis=1, inplace=True)
drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
drifter.rename(columns={' SST':'temperature (degC)'},inplace=True)
drifter.rename(columns={' BP':'barometric_pressure (hPa)'},inplace=True)
drifter.rename(columns={' BPT':'air temperature (degC)'},inplace=True)
drifter['Unit_ID'] = b
drifter = drifter.sort_index() # NOTE index is reversed in raw
drifter = drift_speed(drifter)
drifter = coordinates(drifter)
drifter = drifter.resample('1H').mean() # from 30 min
drifter = drifter.fillna(method='ffill') # interpolate missing data
drifter.to_csv(ODIR+b+'.csv')

#%% SAWS ISVP - 2019 spring cruise
# These data contain missing timestamps.
latname ='latitude'
lonname ='longitude'
timename ='time'

buoy = ['ISVP4','ISVP5','ISVP6']
file = ['300234066433050.xlsx',
        '300234066433051.xlsx',
        '300234066433052.xlsx']

for b,f in zip(buoy,file):
    drifter = pd.read_excel(BDIR+f,parse_dates=True)
    drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
    drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]') # change the fmt of datetimeindex
    drifter.set_index('time',inplace=True)
    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
    drifter.rename(columns={'sst':'temperature (degC)'},inplace=True)
    drifter.rename(columns={'slp':'barometric_pressure (hPa)'},inplace=True)        
    drifter['Unit_ID'] = b
    drifter = drift_speed(drifter)
    drifter = coordinates(drifter)
    drifter = drifter.resample('1H').mean()  
    drifter = drifter.fillna(method='ffill') # interpolate missing data
    drifter.to_csv(ODIR+b+'.csv')

#%% SWIFT buoys - 2019 winter and spring cruises
#   these buoys were deployed and retrieved in winter
#   and then re-deployed in spring and collected.
#   (Thomson et al. 2023). 
#   They have been already processed and downloaded as 
#   .mat files. This allows the buoys to be read into python. 

# convert the matlab time to readable python datetime
def datenum_to_datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60
    return datetime.fromordinal(int(datenum)) \
           + timedelta(days=int(days)) \
           + timedelta(hours=int(hours)) \
           + timedelta(minutes=int(minutes)) \
           + timedelta(seconds=round(seconds)) \
           - timedelta(days=366)
         
fmt = '%Y-%m-%d %H:%M:%S' 
buoy = ['S20_winter','S20_spring',
        'S21_winter','S21_spring']
file = ['SWIFT20_winter.mat','SWIFT20_spring.mat',
        'SWIFT21_winter.mat','SWIFT21_spring.mat']

for b,f in zip(buoy,file):
    # load the .mat files  
    fname = pjoin(BDIR+f)         
    contents = sio.loadmat(fname)
    teststruct = contents['SWIFT']
    teststruct.dtype
    # read in the variables
    Hs = np.transpose(teststruct['sigwaveheight'].astype(float))
    lat = np.transpose(teststruct['lat'].astype(float))
    lon = np.transpose(teststruct['lon'].astype(float))
    matlab_time = np.transpose(teststruct['time'])
    time = np.zeros(len(matlab_time)).astype(datetime) 
    # matlab time to datetime
    i=0
    for i in range(len(matlab_time)):
        time[i] = datenum_to_datetime(float(matlab_time[i])).strftime(fmt)
        i=i+1
    # create dataframe for the swift buoys
    drifter = pd.DataFrame(data=[time, lat, lon, Hs]).T
    drifter.columns =['time', 'latitude (deg)', 'longitude (deg)', 'sig wave height (m)']
    # remove the brackets around the values in each column
    drifter['latitude (deg)'] = drifter['latitude (deg)'].str.get(0)
    drifter['longitude (deg)'] = drifter['longitude (deg)'].str.get(0)
    drifter['sig wave height (m)'] = drifter['sig wave height (m)'].str.get(0)
    drifter['time'] = drifter['time'].astype('datetime64[ns]')
    drifter.set_index('time',inplace=True)
    drifter['Unit_ID'] = b
    drifter = drift_speed(drifter)
    drifter = coordinates(drifter)
    drifter.to_csv(ODIR+b+'.csv')
    
#%% SHARC buoys - 2022 winter cruise
#latname ='GPS Latitude'
#lonname ='GPS Longitude'
#timename ='UTC Time'   

#buoy = ['SHARC_B1','SHARC_B2', 'SHARC_B3', 'SHARC_B4',
#        'SHARC_B5', 'SHARC_B6']
#file = ['SHARC_SB1_ICE21.csv','SHARC_SB2_ICE22.csv', 
#        'SHARC_SB3_ICE23.csv','SHARC_SB4_ICE11.csv', 
#        'SHARC_SB5_ICE12.csv','SHARC_SB6_ICE00.csv']

#fmt = '%d/%b/%Y %H:%M:%S' # different format to the other buoy clusters
#for b,f in zip(buoy,file):
#    drifter = pd.read_csv(BDIR+f,parse_dates=True)
#    drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
#    drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]') # change the fmt of datetimeindex
#    drifter.set_index('time',inplace=True)
#    drifter = drifter[::-1]  
#    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
#    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
#    drifter = drifter.drop(drifter[drifter['latitude (deg)']==0].index)
##    drifter.rename(columns={'sst':'temperature (degC)'},inplace=True)
##    drifter.rename(columns={'slp':'barometric_pressure (hPa)'},inplace=True)        
#    drifter['Unit_ID'] = b
#    drifter = drift_speed(drifter)
#    drifter = coordinates(drifter)
#    drifter.to_csv(ODIR+b+'.csv')

#%% wave buoys - 2022 winter cruise
#latname ='Latitude (deg)'
#lonname ='Longitude (deg)'
#timename ='Calc_time (UTC)'   

#buoy = ['wave_B1','wave_B2', 'wave_B3']
#file = ['wave_lp2.csv','wave_lp4.csv', 
#        'wave_lp8.csv']

#fmt = '%Y.%m.%d %H:%M:%S' # different format to the other buoy clusters
#for b,f in zip(buoy,file):
#    drifter = pd.read_csv(BDIR+f,parse_dates=True)
#    drifter['time'] = pd.to_datetime(drifter[timename],format=fmt)
#    drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]') # change the fmt of datetimeindex
#    drifter.set_index('time',inplace=True)
#    drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
#    drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
##    drifter.rename(columns={'sst':'temperature (degC)'},inplace=True)
##    drifter.rename(columns={'slp':'barometric_pressure (hPa)'},inplace=True)        
#    drifter['Unit_ID'] = b
#    drifter = drift_speed(drifter)
#    drifter = coordinates(drifter)
#    drifter.to_csv(ODIR+b+'.csv')
    
#%% box buoy - 2022 winter cruise

#BFILE ='box_buoy.csv'
#latname ='lat'
#lonname ='lon'
#timename ='time'
#ID = 'box_buoy'

#fmt = "%Y-%m-%dT%H:%M:%S"
#drifter = pd.read_csv(BDIR+BFILE)
#drifter['time'] = pd.to_datetime(drifter[timename])#,format=fmt)
#drifter.time = drifter['time'].dt.strftime('%Y-%m-%d %H:%M:%S').astype('datetime64[ns]') # change the fmt of datetimeindex
#drifter.set_index('time',inplace=True)
#drifter.rename(columns={latname:'latitude (deg)'},inplace=True)
#drifter.rename(columns={lonname:'longitude (deg)'},inplace=True)
#drifter.drop(timename, axis=1, inplace=True)
#drifter['Unit_ID'] = ID
#drifter = drift_speed(drifter)
#drifter = coordinates(drifter)
#drifter.to_csv(ODIR+ID+'.csv')

# end of code
