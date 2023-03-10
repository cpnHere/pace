"""
********************************************
 Chamara Rajapakshe (2023)
 chamara.rajapakshe@nasa.gov/charaj1@umbc.edu
********************************************
Collocate MODIS and ARM ground-based observations
- To compare cloud water path between MODIS and ARM mwrret products
"""
import numpy as np
import time,sys,os,glob
from pyhdf.SD import SD,SDC
import pandas as pd
import xarray as xr
from calendar import monthrange
import datetime as dtm
from cpnCommonlib import vprint, progress_bar, haversine
from cpnMODISlib import readvalue, get_doy
from math import pi
def parallax_correction(_sZA,_sAA,_cth,_lat,H=705e3,R_e = 6371e3):
    """
    Perform correction for parallax effect
    _sZA = sensor zenith angle
    _sAA = sensor azimuth angle
    _cth = cloud top height in m
    _lat = latitude of the considering location
    H    = MODIS terra height in m
    R_e  = Earth radius in m
    """
    tan_a = H*np.tan(_sZA/180.*pi/(H-_cth))
    dx = _cth*tan_a*np.sin(_sAA/180.*pi)
    dy = _cth*tan_a*np.cos(_sAA/180.*pi)
    dLon = dx/R_e/np.cos(_lat/180.*pi) #longitudinal displacement
    dLat = dy/R_e #latitudinal displacement
    return dLon,dLat
def mod06_l2_arm_mwrret(yr,mn,e_max=90,verbose=False):
    """
    Collocate and save MOD06_L2 cloud water path (cwp) and ARM cwp
    yr: year
    mn: month
    e_max: maximum number of data entries per month (assuming max three MODIS granules per day)
    """
    paths ={'mod06':'data/modis/MOD06_L2/region-100_34_-95_39/%d'%yr,\
            'mod03':'data/modis/MOD03/region-100_34_-95_39/%d'%yr,\
            'arm'  :'data/arm/MWRRET/%d/ascii-csv'%yr}
    out_d = {'year'      :np.zeros(e_max,dtype=int),'doy':np.zeros(90,dtype=int),\
             'mod_time'  :np.zeros(e_max,dtype=float),\
             'mod_cwp'   :np.zeros(e_max,dtype=float),\
             'mod_cwp_un':np.zeros(e_max,dtype=float),\
             'arm_time'  :np.zeros(e_max,dtype=float),\
             'arm_cwp'   :np.zeros(e_max,dtype=float),\
             'arm_cwp_un':np.zeros(e_max,dtype=float)} # output data
    arm_lon,arm_lat = -97.458,36.605 # longitude and latitude of the ground-based observation (ARM)
    out_file = "data/outputs/MOD06_L2-arm_MWRRET_-100_34_-95_39_%d%02d_plx.csv"%(yr,mn) # output file name
    t0 = time.time()
    e_i = 0 #data entry row index
    n_d = monthrange(yr,mn)[1]
    N = 0 # for number of processed files test
    for i in np.arange(1,n_d+1,1):
        hdf = {}
        doy = get_doy(yr,mn,i)
        for m6_file in glob.glob(paths['mod06']+'/MOD06_L2.A%d%03d.????.061.?????????????.hdf'%(yr,doy)):
            tm_s = m6_file.split('/',)[-1].split('.')[2] # time string ex. '1745'
            mo06_f = m6_file
            mo03_f = glob.glob(paths['mod03']+'/MOD03.A%d%03d.%s.061.?????????????.hdf'%(yr,doy,tm_s))[0]
            try:
                arm_f  = glob.glob(paths['arm']+'/sgpmwrret1liljclouC1.c2.%d%02d%02d.??????.csv'%(yr,mn,i))[0]
            except Exception as e:
                try:
                    arm_f = glob.glob(paths['arm']+'/sgpmwrret1liljclouC1.c2.%d%02d%02d.??????.custom.csv'%(yr,mn,i))[0]
                except Exception as e:
                    print(e)
                    continue
            vprint("mo06_f: %s"%mo06_f,verbose)
            vprint("mo03_f: %s"%mo03_f,verbose)
            vprint("arm_f:  %s"%arm_f ,verbose)
            hdf['mo06'] = SD(mo06_f,SDC.READ)
            hdf['mo03'] = SD(mo03_f,SDC.READ)
            cwp   = readvalue(hdf['mo06'],'Cloud_Water_Path') # g/m^2
            cwp_u = readvalue(hdf['mo06'],'Cloud_Water_Path_Uncertainty') # percent
            cth   = readvalue(hdf['mo06'],'Cloud_Top_Height') # m (rounded to nearest 50m) 
            m6_lat = hdf['mo06'].select('Latitude').get() # to get cth
            m6_lon = hdf['mo06'].select('Longitude').get() # to get cth
            #Cloud Water Path Relative Uncertainty (Percent)from both best points and points identified as cloud edge at 1km resolution or partly cloudy at 250m based on the Cloud_Water_Path result
            lat = hdf['mo03'].select('Latitude').get()
            lon = hdf['mo03'].select('Longitude').get()
            sZA = readvalue(hdf['mo03'],'SensorZenith',verbose) # sensor zenith angle
            sAA = readvalue(hdf['mo03'],'SensorAzimuth',verbose) # sensor azimuth angle
            #finding closest lat-lon
            vprint("Warning!Altitude ignored when finding in collocation",verbose)
            val = haversine(arm_lon,arm_lat,lon.reshape(-1),lat.reshape(-1))
            val_i = val.argmin() # selected value's index
            _sZA = sZA.reshape(-1)[val_i]
            _sAA = sAA.reshape(-1)[val_i]
            _val = haversine(arm_lon,arm_lat,m6_lon.reshape(-1),m6_lat.reshape(-1))
            _val_i = _val.argmin()
            _cth = cth.reshape(-1)[_val_i]
            dLon,dLat = parallax_correction(_sZA,_sAA,_cth,arm_lat)
            val = haversine(arm_lon+dLon,arm_lat+dLat,lon.reshape(-1),lat.reshape(-1))
            val_i = val.argmin() # selected value's index
            arm = pd.read_csv(arm_f)
            mod_time = 60.*int(tm_s[0:2])+int(tm_s[2:]) #in minutes
            arm_time = np.array([60.*int(t.split(' ')[1].split(':')[0])+int(t.split(' ')[1].split(':')[1]) for t in arm.time]) # in minutes
            arm_i = np.argmin(abs(mod_time-arm_time))
            out_d['year'      ][e_i] = yr
            out_d['doy'       ][e_i] = doy
            out_d['mod_time'  ][e_i] = mod_time                 # minutes from the midnight
            out_d['mod_cwp'   ][e_i] = cwp.reshape(-1)[val_i]   # g/m^2
            out_d['mod_cwp_un'][e_i] = cwp_u.reshape(-1)[val_i] # percent
            out_d['arm_time'  ][e_i] = arm_time[arm_i]          # minutes from the midnight
            out_d['arm_cwp'   ][e_i] = arm.be_lwp[arm_i]        # best estimate of the mass content of cloud liquid water [g/m^2]
            out_d['arm_cwp_un'][e_i] = arm.qc_be_lwp[arm_i]     # 0 means no failed QC test
            e_i += 1
            if e_i >= e_max:
                print("ERROR! e_i exceeded e_max: %d"%e_max)
            N += 1
        progress_bar(i,1,n_d,1,t0)
        vprint('',verbose)
    if N>0:
        print()
        out_d = pd.DataFrame.from_dict(out_d)
        out_d = out_d.drop(index=list(np.arange(e_i,e_max,1)))
        out_d.to_csv(out_file)
    else:
        print("ERROR! No files found!")
if __name__ == '__main__':
    year  = int(sys.argv[1])
    month = int(sys.argv[2])
    mod06_l2_arm_mwrret(year,month)
