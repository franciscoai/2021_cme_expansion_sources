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

ROT_DIFF = 140 # defines the limit to correct the tilt angle in case the diference is larger than this value

# aux functions
def correct_ang(ang_to_correct, reference_ang, do=True):
    ''' 
    Returns a scalar: +180, 180 or 0, used to correct ang_to_correct in case its 
    diference is wrt reference_ang is larger than ROT_DIFF deg
    '''
    if (np.abs(ang_to_correct - reference_ang) > ROT_DIFF) & do:
        if ang_to_correct < 0:
            return 180
        else:
            return -180   
    else:
        return 0      

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
do_tilt_corr = True

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

# plots
os.makedirs(opath, exist_ok=True)

# extracts the following keys=['ERUPTIONDATE', 'LAT','LON', 'ROT','HAN', 'RAT'] from gcs and saves them in a .csv file
last_gcs_value =[]
for cd in gcs_dates:
    print('Extracting GCS data for date...' + cd)
    gcs_key = ['ERUPTIONDATE', 'LAT','LON', 'ROT','HAN', 'RAT','HGT']       
    gcs_idx = gcs_dates.index(cd)
    gcs_val = []
    for sav in gcs[gcs_idx]:
        gcs_val.append([sav.sgui[key] for key in gcs_key])
    gcs_val = np.array(gcs_val)[:,:,0]
    #convert from rad to deg only  LAN, LON, ROT and HAN
    for i in range(1,5):
        gcs_val[:,i] = np.rad2deg(gcs_val[:,i].astype(float))
    # sort by ERUPTIONDATE
    gcs_times = np.array([dt.datetime.strptime(sav.sgui['ERUPTIONDATE'][0].decode('UTF-8'),'%Y-%m-%dT%H:%M:%S.%f') for sav in gcs[gcs_idx]])
    ind = np.argsort(gcs_times)
    gcs_times = gcs_times[ind]
    gcs_val = gcs_val[ind]
    np.savetxt(opath+'/'+cd+'_gcs.csv', gcs_val, delimiter=',', header='DATE, LAT, LON, ROT, HAN, RAT, HGT', fmt='%s')
    last_gcs_value.append(gcs_val[-1,:])
# saves all the last values in a .csv file
np.savetxt(opath+'/last_gcs_values.csv', np.array(last_gcs_value), delimiter=',', header='DATE, LAT, LON, ROT, HAN, RAT, HGT', fmt='%s')

# GCS ratio of R/OH (Axial_Radius/Apex)
for cd in gcs_dates:
    print('Plotting R, OH for date...' + cd)
    gcs_idx = gcs_dates.index(cd)
    k = np.array([float(sav.sgui['RAT']) for sav in gcs[gcs_idx]])
    alpha = np.array([float(sav.sgui['HAN']) for sav in gcs[gcs_idx]])
    oh = np.array([float(sav.sgui['HGT']) for sav in gcs[gcs_idx]])
    h_leg = oh * (1-k) / (1/np.cos(alpha) + np.tan(alpha))
    beta = h_leg / np.cos(alpha)
    rho = h_leg * np.tan(alpha)
    r = k*(beta+rho)/(1-k**2)
    gcs_times = np.array([dt.datetime.strptime(sav.sgui['ERUPTIONDATE'][0].decode('UTF-8'),'%Y-%m-%dT%H:%M:%S.%f') for sav in gcs[gcs_idx]])
    ind = np.argsort(gcs_times)
    gcs_times = gcs_times[ind]
    oh = oh[ind]
    r = r[ind]
    h_leg = h_leg[ind]
    plt.plot(gcs_times, oh, '*k', label='OH')
    plt.plot(gcs_times, r, 'sk', label='R')
    oc = oh -r
    plt.plot(gcs_times, oc, '^k', label='OH-R')
    ocm1 = oh -r -1
    plt.plot(gcs_times, ocm1, 'ok', label='OH-R-1')
    # horizontal tick labes rotated 45 deg
    plt.xticks(rotation=45)    
    plt.title(cd)
    plt.xlabel('Date')
    plt.ylabel('R, OH, OH-R, OH-R-1')
    plt.minorticks_on()
    plt.tick_params('both', which='both')
    plt.grid(which='both')
    plt.legend()
    # adds another y axis to the right
    ax2 = plt.gca().twinx()
    roocm1 = r / h_leg
    ax2.plot(gcs_times, roocm1, '-+k')
    ax2.set_ylabel('R/(OH-R-1)')
    plt.tight_layout()   
    plt.savefig(opath+'/'+cd+'_R_OH.png')
    plt.close()
    # saves value in a .csv file
    np.savetxt(opath+'/'+cd+'_R_OH.csv', np.array([gcs_times, r, oh,oc,ocm1,roocm1]).T, delimiter=',', header='DATE, R, OH, OH-R, OH-R-1, R/(OH-R-1)', fmt='%s')
    
# plots GCS 'RAT' only
for cd in gcs_dates:
    print('Plotting RAT for date...' + cd)
    gcs_key = 'RAT'        
    gcs_idx = gcs_dates.index(cd)
    gcs_val = np.array([float(sav.sgui[gcs_key]) for sav in gcs[gcs_idx]])
    gcs_times = np.array([dt.datetime.strptime(sav.sgui['ERUPTIONDATE'][0].decode('UTF-8'),'%Y-%m-%dT%H:%M:%S.%f') for sav in gcs[gcs_idx]])
    ind = np.argsort(gcs_times)
    gcs_times = gcs_times[ind]
    gcs_val = gcs_val[ind] 
    plt.plot(gcs_times, gcs_val, '*k')
    plt.title(cd)
    plt.ylabel('RAT')
    plt.xlabel('Date')
    plt.tight_layout()
    plt.minorticks_on()
    plt.tick_params('both', which='both')
    plt.grid(which='both')
    plt.savefig(opath+'/'+cd+'_RAT.png')
    plt.close()

# tilt angles
gcs_ind_to_use=0#-1 # which gcs value to use in the scatter plot, -1 is the last one in time
arcade_ind_to_use=0#-1 # which arcade value to use in the scatter plot, -1 is the last one in time

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
        gcs_key = 'ROT'
        gcs_idx = gcs_dates.index(d_str)
        gcs_val = np.rad2deg([float(sav.sgui[gcs_key]) for sav in gcs[gcs_idx]]).astype(float)
        gcs_times = np.array([dt.datetime.strptime(sav.sgui['ERUPTIONDATE'][0].decode('UTF-8'),'%Y-%m-%dT%H:%M:%S.%f') for sav in gcs[gcs_idx]])
        ind = np.argsort(gcs_times)
        gcs_times = gcs_times[ind]
        gcs_val = gcs_val[ind]         
        plt.plot(gcs_times, gcs_val, '*k', label='GCS')  
        x1 = d.to_pydatetime()
        y1 = df_fil.loc[df_fil['date'] == d, 'tilt_angle']         
        if len(y1) > 0:
            y1 = np.array(y1).astype(float)[0]
            y1 += correct_ang(np.median(y1),np.median(gcs_val), do=do_tilt_corr)
            plt.plot(x1, y1, 'sk', label='FIL')
            gcs_vs_fil_all.append([gcs_val[arcade_ind_to_use], y1, d_str])
            #gcs_times = np.concatenate((gcs_times, x1), axis=None)
            #gcs_val = np.concatenate((gcs_val, y1), axis=None)
        x2 = df_ar.loc[df_ar['ar_time_tilt'].dt.date == d.date(), 'ar_time_tilt']
        y2 = df_ar.loc[df_ar['ar_time_tilt'].dt.date == d.date(), 'ar_tilt']
        if len(y2) > 0:
            y2 = np.array(y2).astype(float)[0]
            y2 += correct_ang(np.median(y2),np.median(gcs_val), do=do_tilt_corr)
            # put y2 in -90 to 90 range
            if y2 > 90:
                y2 = y2 - 180
            if y2 < -90:
                y2 = y2 + 180
            plt.plot(x2, y2, '^k', label='AR')
            gcs_vs_ar_val.append([gcs_val[arcade_ind_to_use], y2, d_str])
            #gcs_times = np.concatenate((gcs_times, [i.to_pydatetime() for i in x2]), axis=None)
            #gcs_val = np.concatenate((gcs_val, y2), axis=None)       
        #adds arcade mean titls
        x3 = df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'date']
        y3 = np.array(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'tilt mean [deg]']).astype(float)
        if len(y3) > 0:     
            # corrects if the diference is > 140 deg
            y3 += correct_ang(np.median(y3),np.median(gcs_val), do=do_tilt_corr)
            # gcs rotation vs arcade
            cdiff = float(np.abs(gcs_val[arcade_ind_to_use] - y3[arcade_ind_to_use]))
            if cdiff > ROT_DIFF:
                if do_tilt_corr:
                    print('****************ERROR, you should not be here. Date...' + d_str)
                cdiff= 180 - cdiff
            diff_all.append(cdiff)
            diff_dates.append(d)
            gcs_param = np.array((df_ar.loc[df_ar['datetimes'].dt.date == d.date(), 'gcs_lat_vel']))
            if len(gcs_param) > 0:
                gcs_param = np.nanmean(gcs_param)
                diff_gcs_param_all.append(float(gcs_param))
                plt.plot(x3, y3, 'ok', label='Arcades')
                gcs_vs_arcades_all.append([gcs_val[arcade_ind_to_use], y3[arcade_ind_to_use], d_str])            
            else:
                diff_gcs_param_all.append(None)
                gcs_vs_arcades_all.append([None, None]) 
                print('No gcs_lat_vel for date...' + d_str)
                pass
        #ind = np.argsort(gcs_times)
        #gcs_times = gcs_times[ind]
        #gcs_val = gcs_val[ind]
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
        plt.savefig(opath+'/'+d_str+'_'+gcs_key+'_ang_corr'+str(do_tilt_corr)+'.png')
        plt.close()

#scatter plots
gcs_vs_ar_val_arr = np.array(gcs_vs_ar_val)[:,0:2].astype(float)
gcs_vs_fil_all_arr = np.array(gcs_vs_fil_all)[:,0:2].astype(float)
gcs_vs_arcades_all_arr = np.array(gcs_vs_arcades_all)[:,0:2].astype(float)

plt.scatter(gcs_vs_ar_val_arr[:, 1], gcs_vs_ar_val_arr[:, 0], color='r', label='GCS vs. AR')
plt.scatter(gcs_vs_fil_all_arr[:, 1], gcs_vs_fil_all_arr[:, 0], color='b', label='GCS vs. Filament')
plt.scatter(gcs_vs_arcades_all_arr[:, 1], gcs_vs_arcades_all_arr[:, 0], color='g', label='GCS vs. PEA')

plt.xlabel('Source Region Tilt [deg]')
plt.ylabel('GCS Tilt [deg]') 
plt.tight_layout()
plt.xlim([-110, 110])
plt.ylim([-110, 110])
plt.grid(which='both')
plt.legend()
plt.savefig(opath+'/gcs_tilt_vs_fil-ar-arcade_clean.png') 
# annotate dates in each point
# for i, txt in enumerate(gcs_vs_ar_val):
#     plt.annotate(txt[2], (gcs_vs_ar_val_arr[i, 1], gcs_vs_ar_val_arr[i, 0]))
# for i, txt in enumerate(gcs_vs_fil_all):
#     plt.annotate(txt[2], (gcs_vs_fil_all_arr[i, 1], gcs_vs_fil_all_arr[i, 0]))
for i, txt in enumerate(gcs_vs_arcades_all):
    plt.annotate(txt[2], (gcs_vs_arcades_all_arr[i, 1], gcs_vs_arcades_all_arr[i, 0]))
# grid for minor ticks behind the plot
plt.minorticks_on()
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