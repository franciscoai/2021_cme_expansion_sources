#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compares the properties of various fenomena associated to a CME:
- Magnetic AR morphological and magnetic properties
- Ca filament morphological properties
- White light CME 3D morphological properties from GCS forward modeling

@author: iglesias
"""
from cProfile import label
from genericpath import exists
from turtle import color, title
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import datetime as dt
from scipy.io import readsav


# aux functions
def correct_ang(iang, flip=0):
    # corrects the input angle to be in the range -90 to 90 only
    # param: ang: input angles in deg (numpy array)
    ang = np.array(iang)
    ang[(ang < 20)] += 180
    # if flip == 1:
    #     ang[(ang < 0)] += 180
    # else:
    #     ang[(ang > 90) & (ang <= 270) ] -= 180
    #     ang[(ang > 270)] -= 360 
    return ang


# main
local_path = os.getcwd()
repo_path = local_path # os.path.dirname(os.path.dirname(local_path))
opath = repo_path + '/output_data/compare_all_prop'

# filaments data (fil)
fil_file = repo_path + '/output_data/filaments.csv'

# GCS data (gcs)
gcs_dir = repo_path + '/GCSs'

# Magnetic AR data (ar)
ar_file = repo_path + '/input_data/ar.csv'

#arcades prop
arcades_file = repo_path + '/output_data/arcades/props/all_arcades_props.csv'

# constants
# dates and sources [gcs,fil,ar] to flip tilt angle signs if angle >80 deg
sign_corr_dates = {'20130209': [1, 0, 0], '20130517': [1, 0, 0], '20130502': [0, 0, 1], '20110317': [0, 1, 0]}

##########

# read fil
df_fil = pd.read_csv(fil_file, header=0, delimiter=',')
df_fil['date'] = pd.to_datetime(df_fil['date'], format='%Y-%m-%d %H:%M:%S')

# reads ar
df_ar = pd.read_csv(ar_file, header=0, delimiter=',')
df_ar['ar_time_tilt'] = pd.to_datetime(df_ar['ar_time_tilt'], format='%Y/%m/%d %H:%M:%S')

# reads arcades 
df_arcades = pd.read_csv(arcades_file, header=0, delimiter=',')
df_arcades['date'] = pd.to_datetime(df_arcades['date'], format='%Y-%m-%d %H:%M:%S')

# reads only the gcs .sav files for the dates in fil
gcs = []
all_gcs_dates = os.listdir(gcs_dir)
all_gcs_dates = [i for i in all_gcs_dates if i.startswith('GCS_')]
gcs_dates = []
for d in df_fil['date']:
    for gcs_d in all_gcs_dates:
        if str.split(gcs_d, '_')[1] == d.strftime('%Y%m%d'):
            print('Reading gcs for date...' + gcs_d)
            gcs_savs = os.listdir(os.path.join(gcs_dir, gcs_d))
            gcs_savs = [i for i in gcs_savs if i.endswith('.sav')]
            gcs_savs_data = []
            for gcs_sav in gcs_savs:
                bla = readsav(os.path.join(gcs_dir, gcs_d, gcs_sav))
                gcs_savs_data.append(bla)
            gcs.append(gcs_savs_data)
            gcs_dates.append(d.strftime('%Y%m%d'))
gcs_keys = bla.sgui.dtype.names  # all keys in sav files
# computes secondary GCS parameters
#  awl=2.*(sgui.HAN + asin(sgui.RAT)) *!RAdeg
#  awd=2.*asin(sgui.RAT) *!RAdeg

# plots
os.makedirs(opath, exist_ok=True)
exit
# tilt angles
df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
gcs_vs_ar_val = []
gcs_vs_fil_all = []
gcs_vs_arcades_all = []
diff_all = [] # last gcs tilt - mean arcade tilt
diff_dates = []
diff_gcs_param_all = []
for d in df_fil['date']:
    d_str = d.strftime('%Y%m%d')
    print('Plotting tilt angles for date...' + d_str)
    if d_str in gcs_dates:
        if d_str in sign_corr_dates:
            sign_corr = sign_corr_dates[d_str]
        else:
            sign_corr = [0,0,0]
        gcs_key = 'ROT'           
        gcs_idx = gcs_dates.index(d_str)
        gcs_val = np.rad2deg([float(sav.sgui[gcs_key]) for sav in gcs[gcs_idx]]).astype(float)
        gcs_times = np.array([dt.datetime.strptime(sav.sgui['ERUPTIONDATE'][0].decode('UTF-8'),
                                                   '%Y-%m-%dT%H:%M:%S.%f') for sav in gcs[gcs_idx]])
        ind = np.argsort(gcs_times)
        gcs_times = gcs_times[ind]
        gcs_val = gcs_val[ind]         
        gcs_val = correct_ang(gcs_val, flip=1)
        plt.plot(gcs_times, gcs_val, '*k', label='GCS')  
        x1 = d.to_pydatetime()
        y1 = df_fil.loc[df_fil['date'] == d, 'tilt_angle']         
        if len(y1) > 0:
            y1 = np.array(y1).astype(float)[0]
            y1 = correct_ang(y1, flip=1)
            plt.plot(x1, y1, 'sk', label='FIL')
            gcs_vs_fil_all.append([gcs_val[0], y1])
            gcs_times = np.concatenate((gcs_times, x1), axis=None)
            gcs_val = np.concatenate((gcs_val, y1), axis=None)
        x2 = df_ar.loc[df_ar['ar_time_tilt'].dt.date == d.date(), 'ar_time_tilt']
        y2 = df_ar.loc[df_ar['ar_time_tilt'].dt.date == d.date(), 'ar_tilt']
        if len(y2) > 0:
            y2 = np.array(y2).astype(float)[0]
            y2 = correct_ang(y2, flip=1)
            plt.plot(x2, y2, '^k', label='AR')
            gcs_vs_ar_val.append([gcs_val[0], y2])
            gcs_times = np.concatenate((gcs_times, [i.to_pydatetime() for i in x2]), axis=None)
            gcs_val = np.concatenate((gcs_val, y2), axis=None)       
        #adds arcade mean titls
        x3 = df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'date']
        y3 = df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'tilt mean [deg]']
        # subtracts 180 if y3 is larger than 90
        y3 = np.array(y3).astype(float)
        y3 = correct_ang(y3, flip=0)
        # gcs rotation vs arcade
        cdiff = float(np.abs(gcs_val[-1] - np.mean(y3)))
        if cdiff > 90:
            cdiff= 180 - cdiff
        diff_all.append(cdiff)
        diff_dates.append(d)
        gcs_param = np.array((df_ar.loc[df_ar['datetimes'].dt.date == d.date(), 'gcs_lat_vel']))
        if len(gcs_param) > 0:
            gcs_param = np.nanmean(gcs_param)
        diff_gcs_param_all.append(float(gcs_param))
        plt.plot(x3, y3, 'ok', label='Arcades')
        gcs_vs_arcades_all.append([gcs_val[0], y3[0]])
        ind = np.argsort(gcs_times)
        gcs_times = gcs_times[ind]
        gcs_val = gcs_val[ind]
        #plt.plot(gcs_times, gcs_val, '--k')
        plt.title(str(d))
        plt.ylabel('Tilt angle [deg]')
        plt.xlabel('Date')        
        plt.legend()
        plt.tight_layout()
        plt.minorticks_on()
        #plt.ylim([-180., 180.])
        plt.tick_params('both', which='both')
        plt.grid(which='both')
        plt.savefig(opath+'/'+d_str+'_'+gcs_key+'.png')
        plt.close()

gcs_vs_fil_all = np.array(gcs_vs_fil_all)
gcs_vs_ar_val = np.array(gcs_vs_ar_val)
gcs_vs_arcades_all = np.array(gcs_vs_arcades_all)
plt.scatter(gcs_vs_fil_all[:, 0], gcs_vs_fil_all[:, 1], color='b', label='gcs vs fil')
plt.scatter(gcs_vs_ar_val[:, 0], gcs_vs_ar_val[:, 1], color='r', label='gcs vs ar')
plt.scatter(gcs_vs_arcades_all[:, 0], gcs_vs_arcades_all[:, 1], color='g', label='gcs vs arcades')
plt.ylabel('fil/ar [deg]')
plt.xlabel('gcs [deg]')  
#plt.xlim([-90, 90])
#plt.ylim([-90, 90])
plt.legend()
plt.savefig(opath+'/gcs_tilt_vs_fil-ar-arcade.png')
plt.close()

#plots diff vs date
diff_all = np.array(diff_all)
plt.scatter(diff_gcs_param_all, diff_all)
plt.savefig(opath+'/diff_titl_gcs_vs_gcs_param.png')

# plots df_ar['gcs_axial_vel'] vs df_arcades['width fit vel [km/s]'] for each date
print('Plotting velocities...')
all_x=[]
all_y=[]
for d in df_ar['datetimes']:
    d_str = d.strftime('%Y%m%d')
    x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_axial_vel']))
    if len(x) > 1:
        x = np.nanmean(x)
    x = float(x)
    y = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
    all_x.append(x)
    all_y.append(y)
plt.scatter(all_x, all_y)
plt.xlabel('gcs_axial_vel [km/s]')
plt.ylabel('arcade width fit vel [km/s]')
plt.savefig(opath+'/gcs_axial_vel_vs_width_fit_vel.png')
plt.close()

# plots df_ar['gcs_radial_vel_at6'] vs df_arcades['width fit vel [km/s]'] for each date
df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
all_x=[]
all_y=[]
for d in df_ar['datetimes']:
    d_str = d.strftime('%Y%m%d')
    x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_radial_vel_at6']))
    if len(x) > 1:
        x = np.nanmean(x)
    x = float(x)
    y = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
    all_x.append(x)
    all_y.append(y)
plt.scatter(all_x, all_y)
plt.xlabel('gcs_radial_vel_at6 [km/s]')
plt.ylabel('arcade width fit vel [km/s]')
plt.savefig(opath+'/gcs_radial_vel_at6_vs_width_fit_vel.png')
plt.close()

# plots df_ar['gcs_radial_vel_at6'] vs df_arcades['length fit vel [km/s]'] for each date
df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
all_x=[]
all_y=[]
for d in df_ar['datetimes']:
    d_str = d.strftime('%Y%m%d')
    x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_radial_vel_at6']))
    if len(x) > 1:
        x = np.nanmean(x)
    x = float(x)
    y = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length fit vel [km/s]'])
    all_x.append(x)
    all_y.append(y)
plt.scatter(all_x, all_y)
plt.xlabel('gcs_radial_vel_at6 [km/s]')
plt.ylabel('arcade length fit vel [km/s]')
plt.savefig(opath+'/gcs_radial_vel_at6_vs_length_fit_vel.png')
plt.close()

# plots df_ar['gcs_lat_vel'] vs df_arcades['width fit vel [km/s]'] for each date
df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
all_x=[]
all_y=[]
for d in df_ar['datetimes']:
    d_str = d.strftime('%Y%m%d')
    x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_lat_vel']))
    if len(x) > 1:
        x = np.nanmean(x)
    x = float(x)
    y = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
    all_x.append(x)
    all_y.append(y)
plt.scatter(all_x, all_y)
plt.xlabel('gcs_lat_vel [km/s]')
plt.ylabel('arcade width fit vel [km/s]')
plt.savefig(opath+'/gcs_lat_vel_vs_width_fit_vel.png')
plt.close()


# plots df_ar['gcs_awl_awd_ratio'] vs df_arcades['width fit vel [km/s]'] for each date
df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
all_x=[]
all_y=[]
for d in df_ar['datetimes']:
    d_str = d.strftime('%Y%m%d')
    x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_awl_awd_ratio']))
    if len(x) > 1:
        x = np.nanmean(x)
    x = float(x)
    y = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]']) 
    y /= np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length fit vel [km/s]'])
    all_x.append(x)
    all_y.append(y)
plt.scatter(all_x, all_y)
plt.xlabel('gcs_awl_awd_ratio')
plt.ylabel('arcade width fit vel /length fit vel')
plt.savefig(opath+'/gcs_awl_awd_ratio_vs_width_over_length_fit_vel.png')
plt.close()