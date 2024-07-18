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
import scipy.stats
# aux functions


def correct_ang(iang, flip=0):
    # corrects the input angle to be in the range -0 to 180 only
    # param: ang: input angles in deg (numpy array)
    ang = np.array(iang)
    ang[iang > 270] -= 360
    ang[iang > 180] -= 180
    ang[iang > 90] -= 180
    ang[iang < -270] += 360
    ang[iang < -180] += 180
    ang[iang < -90] += 180
    if flip == 1:
        ang[np.abs(ang) > 75] *= -1
    #ang[ang < 0] += 180
    return np.abs(ang)

def drop_by_index(lista_original,indices_to_drop):
    """
    Esta funcion toma 2 listas, la lista que quiero filtrar. La lista de indices_to_drop contiene los indices de la primera que
    quiero eliminar.
    La salida es entonces la lista original sin lista_original[indices_to_drop]
    """
    if not indices_to_drop:
        return lista_original
    if indices_to_drop:
        return [element for i, element in enumerate(lista_original) if i not in indices_to_drop]

# main
local_path = os.getcwd()
repo_path = local_path # os.path.dirname(os.path.dirname(local_path))
opath = repo_path + '/output_data/compare_all_prop'
hpath = repo_path + '/output_data/compare_all_prop_hebe'

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

# set requiered plots
plot_tilt = False
pea_prop_over_sep_speed_vs_gcs_lat_over_axial_vel = True
pea_prop_sep_speed_ratio_vs_gcs_awl_awd_ratio = True
pea_prop_speed_vs_gcs_awl = True
pea_prop_speed_vs_gcs_lat_speed = True
pea_sep_speed_vs_gcs_awd = True
pea_sep_speed_vs_gcs_axial_speed = True
pea_length_vs_gcs_lat_speed = False
#lista de CMEs que no deben plotearse
#poner el numero de ID
#Si no queremos rechazar ninguna, poner --> None
list_rejected_cmes_id = None 
#por ejemplo si quiero quitar el evento 1
#list_rejected_cmes_id = [0]

#Mas colores en https://i.stack.imgur.com/nCk6u.jpg

if plot_tilt:
    # tilt angles
    gcs_vs_ar_val = []
    gcs_vs_fil_all = []
    gcs_vs_arcades_all = []
    for d in df_fil['date']:
        d_str = d.strftime('%Y%m%d')
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
            gcs_val = correct_ang(gcs_val, flip=sign_corr[0])
            plt.plot(gcs_times, gcs_val, '*k', label='GCS')
            x1 = d.to_pydatetime()
            y1 = df_fil.loc[df_fil['date'] == d, 'tilt_angle']
            if len(y1) > 0:
                y1 = np.array(y1).astype(float)[0]
                y1 = correct_ang(y1, flip=sign_corr[1])
                plt.plot(x1, y1, 'sk', label='FIL')
                gcs_vs_fil_all.append([gcs_val[0], y1])
                gcs_times = np.concatenate((gcs_times, x1), axis=None)
                gcs_val = np.concatenate((gcs_val, y1), axis=None)
            x2 = df_ar.loc[df_ar['ar_time_tilt'].dt.date == d.date(), 'ar_time_tilt']
            y2 = df_ar.loc[df_ar['ar_time_tilt'].dt.date == d.date(), 'ar_tilt']
            if len(y2) > 0:
                y2 = np.array(y2).astype(float)[0]
                y2 = correct_ang(y2, flip=sign_corr[2])
                plt.plot(x2, y2, 'ok', label='AR')
                gcs_vs_ar_val.append([gcs_val[0], y2])
                gcs_times = np.concatenate((gcs_times, [i.to_pydatetime() for i in x2]), axis=None)
                gcs_val = np.concatenate((gcs_val, y2), axis=None)
            #adds arcade mean titls
            x3 = df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'date']
            y3 = df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'tilt mean [deg]']
            # subtracts 180 if y3 is larger than 90
            y3 = np.array(y3).astype(float)
            y3 = correct_ang(y3, flip=0)
            plt.plot(x3, y3, 'ok', label='Arcades')
            gcs_vs_arcades_all.append([gcs_val[0], y3[0]])
            ind = np.argsort(gcs_times)
            gcs_times = gcs_times[ind]
            gcs_val = gcs_val[ind]
            plt.plot(gcs_times, gcs_val, '--k')
            plt.title(str(d))
            plt.ylabel('Tilt angle [deg]')
            plt.xlabel('Date')        
            plt.legend()
            plt.tight_layout()
            plt.minorticks_on()
            plt.ylim([-180., 180.])
            plt.tick_params('both', which='both')
            plt.grid(which='both')
            plt.savefig(opath+'/'+d_str+'_'+gcs_key+'.png')
            plt.close()
    #plot fils vs ar
    gcs_vs_fil_all = np.array(gcs_vs_fil_all)
    gcs_vs_ar_val = np.array(gcs_vs_ar_val)
    gcs_vs_arcades_all = np.array(gcs_vs_arcades_all)
    plt.scatter(gcs_vs_fil_all[:, 0], gcs_vs_fil_all[:, 1], color='b', label='gcs vs fil')
    plt.scatter(gcs_vs_ar_val[:, 0], gcs_vs_ar_val[:, 1], color='r', label='gcs vs ar')
    plt.scatter(gcs_vs_arcades_all[:, 0], gcs_vs_arcades_all[:, 1], color='g', label='gcs vs arcades')
    plt.ylabel('fil/ar [deg]')
    plt.xlabel('gcs [deg]')  
    plt.xlim([-90, 90])
    plt.ylim([-90, 90])
    plt.legend()
    plt.savefig(opath+'/scatter.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots PEA propagation speed / PEA separation speed VS GCS AW_L vs. GCS AW_D for each date
if pea_prop_sep_speed_ratio_vs_gcs_awl_awd_ratio:
    ajuste_lineal3 = True

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_awl']))
        x /= np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_awd']))
        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length fit vel [km/s]'])
        y /= np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS AW$_L$ / AW$_D$', fontsize=14)
    ax.set_ylabel('PEA prop./sep. speed ', fontsize=14)
    ax.set_title(' ', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal3:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="upper right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_prop_sep_speed_ratio_vs_gcs_awl_awd_ratio')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots PEA propagation speed/PEA separation speed VS. GCS lat vel / gcs axial vel  ---> Done Hebe
# y vs x
ajuste_lineal = True
if pea_prop_over_sep_speed_vs_gcs_lat_over_axial_vel:
    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')

        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_lat_vel']))
        x /= np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_axial_vel']))

        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)
        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length fit vel [km/s]'])
        y /= np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]']) 
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]
    #breakpoint()
    
    #se filtran los eventos que queremos descartar
    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value and linear regresion
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS lateral (L) / axial (D) exp. speed', fontsize=14)
    ax.set_ylabel('PEA prop. / sep. speed', fontsize=14)
    ax.set_title(' ', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="upper right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_prop_over_sep_speed_vs_gcs_lat_over_axial_vel.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots PEA propagation speed vs. GCS AW_L for each date  --> Done Hebe
if pea_prop_speed_vs_gcs_awl:
    ajuste_lineal1 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_awl']))

        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS AW$_L$', fontsize=14)
    ax.set_ylabel('PEA propagation speed [km s$^{-1}$]', fontsize=14)
    #ax.set_title('Title', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal1:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="upper right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_prop_speed_vs_gcs_awl.png')
    plt.close()
#------------------------------------------------------------------------------------------------------
# plots PEA propagation speed vs. GCS lateral (L) speed for each date  --> Done Hebe
if pea_prop_speed_vs_gcs_lat_speed:
    ajuste_lineal1 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_lat_vel']))
        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS lateral (L) speed [km s$^{-1}$]', fontsize=14)
    ax.set_ylabel('PEA propagation speed [km s$^{-1}$]', fontsize=14)
    #ax.set_title('Title', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal1:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="upper right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_prop_speed_vs_gcs_lat_speed.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots PEA separation speed vs. GCS AW_D for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_awd']))

        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS AW$_D$', fontsize=14)
    ax.set_ylabel('PEA separation speed [km s$^{-1}$]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="lower right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_sep_speed_vs_gcs_awd.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots mean lenght vs. GCS AW_D for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for id in df_ar['ID']:
        d_str = d.strftime('%Y%m%d')
        d = [i for i in df_ar.loc[df_ar['ID'] == id, 'datetimes']][0]
        x = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awd']))
        x = float(x)
        y = [i for i in df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length2 [arcsec]']]
        if len(y) > 0:
            total_y = 0
            for i in y:
                total_y += np.sum([float(j) for j in str.split(str.split(str.split(i,'[')[1],']')[0])])
            total_y /= len(y)
            all_x.append(x)
            all_y.append(total_y)
        else:
            print('Warning, could not find length1 for case ' + d_str)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS AW$_D$', fontsize=14)
    ax.set_ylabel('PEA mean length [acsec]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="lower right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_mean_length_vs_gcs_awd.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots mean lenght vs. GCS AW_D for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for id in df_ar['ID']:
        d_str = d.strftime('%Y%m%d')
        d = [i for i in df_ar.loc[df_ar['ID'] == id, 'datetimes']][0]
        x = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awd']))
        x = float(x)
        y = [i for i in df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width [arcsec]']]
        if len(y) > 0:
            total_y = 0
            for i in y:
                total_y += np.mean([float(j) for j in str.split(str.split(str.split(i,'[')[1],']')[0])])
            total_y /= len(y)
            all_x.append(x)
            all_y.append(total_y)
        else:
            print('Warning, could not find length1 for case ' + d_str)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS AW$_D$', fontsize=14)
    ax.set_ylabel('PEA mean width [arcsec]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="lower right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_mean_width_vs_gcs_awd.png')
    plt.close()
#------------------------------------------------------------------------------------------------------
# plots mean lenght vs. GCS AW_D for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for id in df_ar['ID']:
        d_str = d.strftime('%Y%m%d')
        d = [i for i in df_ar.loc[df_ar['ID'] == id, 'datetimes']][0]
        x = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_radial_vel_at6']))
        x = float(x)
        y = [i for i in df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width [arcsec]']]
        if len(y) > 0:
            total_y = 0
            for i in y:
                total_y += np.mean([float(j) for j in str.split(str.split(str.split(i,'[')[1],']')[0])])
            total_y /= len(y)
            all_x.append(x)
            all_y.append(total_y)
        else:
            print('Warning, could not find length1 for case ' + d_str)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('gcs_radial_vel_at6 [km s$^{-1}$]', fontsize=14)
    ax.set_ylabel('PEA mean width [arcsec]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="best", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_mean_width_vs_gcs_radial_vel_at6.png')
    plt.close()
    
    #------------------------------------------------------------------------------------------------------
# plots mean lenght vs. GCS AW_D for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for id in df_ar['ID']:
        d_str = d.strftime('%Y%m%d')
        d = [i for i in df_ar.loc[df_ar['ID'] == id, 'datetimes']][0]
        x = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_radial_vel_at2']))
        x = float(x)
        y = [i for i in df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width [arcsec]']]
        if len(y) > 0:
            total_y = 0
            for i in y:
                total_y += np.mean([float(j) for j in str.split(str.split(str.split(i,'[')[1],']')[0])])
            total_y /= len(y)
            all_x.append(x)
            all_y.append(total_y)
        else:
            print('Warning, could not find length1 for case ' + d_str)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('gcs_radial_vel_at2 [km s$^{-1}$]', fontsize=14)
    ax.set_ylabel('PEA mean width [arcsec]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="best", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_mean_width_vs_gcs_radial_vel_at2.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots mean PEA length vs. GCS AW_L/AW_D for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_x_sd = []
    all_y=[]
    all_y_sd = []
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for id in df_ar['ID']:
        d_str = d.strftime('%Y%m%d')
        d = [i for i in df_ar.loc[df_ar['ID'] == id, 'datetimes']][0]
        awl = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awl']))
        #awd = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awd']))
        x = awl #/ awd
        #awd_err = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awd_err']))
        awl_err = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awl_err']))
        # error propagation of the ratio
        x_err = awl_err
        if len(x) > 1:
            x = np.nanmean(x)
            breakpoint()
        x = float(x)
        y1 = [i for i in df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length1 [arcsec]']]
        y2 = [i for i in df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length2 [arcsec]']]
        if len(y1) > 0:
            curr_y1 = []
            curr_y2 = []
            for i in y1:
                curr_y1.append(np.sum([float(j) for j in str.split(str.split(str.split(i,'[')[1],']')[0])]))
            for i in y2:
                curr_y2.append(np.sum([float(j) for j in str.split(str.split(str.split(i,'[')[1],']')[0])]))            
            ym = np.mean([curr_y1,curr_y2])
            ysd = np.std([curr_y1,curr_y2])
            all_x.append(x)
            all_y.append(ym)
            all_y_sd.append(ysd)
            all_x_sd.append(x_err)
        else:
            print('Warning, could not find length1 for case ' + d_str)

    all_x = all_x [0:11]
    all_y = all_y [0:11]
    all_y_sd = all_y_sd [0:11]
    all_x_sd = all_x_sd [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    all_y_sd_filtered = drop_by_index(all_y_sd,list_rejected_cmes_id)
    all_x_sd_filtered = drop_by_index(all_x_sd,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS AW$_L$ [deg]', fontsize=14)
    ax.set_ylabel('PEA mean length [arcsec]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        # scatter plot error bars, and a dot in the mean value. Use thin grey lines styles for  all the error bars
        ax.errorbar(all_x_filtered[contador],all_y_filtered[contador],xerr=all_x_sd_filtered[contador],yerr=all_y_sd_filtered[contador],c='grey',fmt='o',elinewidth=0.5, capsize=2, capthick=0.5)
        # plot the points only on top in z
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],zorder=10)   

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="best", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_mean_length_vs_gcs_awl.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots PEA separation speed vs. GCS radial  (d) speed for each date ---> Done hebe
if pea_sep_speed_vs_gcs_axial_speed:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_radial_vel_at6']))

        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('gcs_radial_vel_at6 [km s$^{-1}$]', fontsize=14)
    ax.set_ylabel('PEA separation speed [km s$^{-1}$]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="lower right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_sep_speed_vs_gcs_radial_vel_at6.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots PEA separation speed vs. GCS radial  (d) speed for each date ---> Done hebe
if pea_sep_speed_vs_gcs_axial_speed:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_radial_vel_at2']))

        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS radial speed at 2 R$_{\odot}$ [km s$^{-1}$]', fontsize=14)
    ax.set_ylabel('PEA separation speed [km s$^{-1}$]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):

        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="lower right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_sep_speed_vs_gcs_radial_vel_at2.png')
    plt.close()
#------------------------------------------------------------------------------------------------------
# plots PEA separation speed vs. GCS radial  (d) speed for each date ---> Done hebe
if pea_sep_speed_vs_gcs_axial_speed:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_mean_acc1to6']))

        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('gcs_mean_acc1to6 [m s$^{-2}$]', fontsize=14)
    ax.set_ylabel('PEA separation speed [km s$^{-1}$]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="lower right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_sep_speed_vs_gcs_mean_acc1to6.png')
    plt.close()


#------------------------------------------------------------------------------------------------------
# plots PEA separation speed vs. GCS axial (d) speed for each date ---> Done hebe
if pea_sep_speed_vs_gcs_axial_speed:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_axial_vel']))

        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'width fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS axial speed [km s$^{-1}$]', fontsize=14)
    ax.set_ylabel('PEA separation speed [km s$^{-1}$]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="lower right", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_sep_speed_vs_gcs_axial_speed.png')
    plt.close()

#     #------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# plots mean PEA length vs. GCS AW_L/AW_D for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_x_sd = []
    all_y=[]
    all_y_sd = []
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for id in df_ar['ID']:
        d_str = d.strftime('%Y%m%d')
        d = [i for i in df_ar.loc[df_ar['ID'] == id, 'datetimes']][0]
        awl = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awl']))
        awd = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awd']))
        x = awl / awd
        awd_err = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awd_err']))
        awl_err = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awl_err']))
        # error propagation of the ratio
        x_err = x * np.sqrt((awd_err/awd)**2 + (awl_err/awl)**2)
        if len(x) > 1:
            x = np.nanmean(x)
            breakpoint()
        x = float(x)
        y1 = [i for i in df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length1 [arcsec]']]
        y2 = [i for i in df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length2 [arcsec]']]
        if len(y1) > 0:
            curr_y1 = []
            curr_y2 = []
            for i in y1:
                curr_y1.append(np.sum([float(j) for j in str.split(str.split(str.split(i,'[')[1],']')[0])]))
            for i in y2:
                curr_y2.append(np.sum([float(j) for j in str.split(str.split(str.split(i,'[')[1],']')[0])]))            
            ym = np.mean([curr_y1,curr_y2])
            ysd = np.std([curr_y1,curr_y2])
            all_x.append(x)
            all_y.append(ym)
            all_y_sd.append(ysd)
            all_x_sd.append(x_err)
        else:
            print('Warning, could not find length1 for case ' + d_str)

    all_x = all_x [0:11]
    all_y = all_y [0:11]
    all_y_sd = all_y_sd [0:11]
    all_x_sd = all_x_sd [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    all_y_sd_filtered = drop_by_index(all_y_sd,list_rejected_cmes_id)
    all_x_sd_filtered = drop_by_index(all_x_sd,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS AW$_L$/AW$_D$', fontsize=14)
    ax.set_ylabel('PEA mean length [arcsec]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        # scatter plot error bars, and a dot in the mean value. Use thin grey lines styles for  all the error bars
        ax.errorbar(all_x_filtered[contador],all_y_filtered[contador],xerr=all_x_sd_filtered[contador],yerr=all_y_sd_filtered[contador],c='grey',fmt='o',elinewidth=0.5, capsize=2, capthick=0.5)
        # plot the points only on top in z
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],zorder=10)   

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="best", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_mean_length_vs_gcs_awl_awd_ratio.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots mean filament length vs. GCS AW_L/AW_D for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for id in df_ar['ID']:
        dd = [i for i in df_ar.loc[df_ar['ID'] == id, 'datetimes']][0]
        x = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awl']))
        x /= np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_awd']))
        
        y = df_fil.loc[df_fil['id'] == id, 'lenght']
        
        if len(y)!=0:
            all_x.append(float(x))
            all_y.append(float(y)/1000)


#    all_x = all_x [0:11]
 #   all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS AW$_L$/AW$_D$', fontsize=14)
    ax.set_ylabel('filament length [Mm]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="best", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_fillength_vs_gcs_awl_awd_ratio.png')
    plt.close()

#------------------------------------------------------------------------------------------------------
# plots magnetic properties vs. GCS acc 1to6 for each date ---> Done hebe
if pea_sep_speed_vs_gcs_awd:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for id in df_ar['ID']:
        d_str = d.strftime('%Y%m%d')
        d = [i for i in df_ar.loc[df_ar['ID'] == id, 'datetimes']][0]
        x = np.array((df_ar.loc[df_ar['ID'] == id, 'gcs_lat_vel']))
        y = [i for i in df_ar.loc[df_ar['ID'] == id, 'ar_unsigned_flux']]
        
        if len(y)!=0 and ~np.isnan(y[0]):
            all_x.append(float(x))
            all_y.append(float(y[0]))
      #      all_y.append(np.fromstring(y[0],dtype='float64'))

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    # linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    # slope = linear_regresion.slope
    # intercept = linear_regresion.intercept
    # pearson = linear_regresion.rvalue
    # r_square = pearson*pearson
    # p_value = linear_regresion.pvalue
    # stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS lateral speed km s$^{-1}$]', fontsize=14)
    ax.set_ylabel('AR unsigned flux [Mx]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="best", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/ar_unsigned_flux_vs_gcs_lat_vel.png')
    plt.close()

 #------------------------------------------------------------------------------------------------------
# plots PEA propagation speed vs. GCS radial (d) speed for each date ---> Done hebe
if pea_sep_speed_vs_gcs_axial_speed:
    ajuste_lineal2 = False

    df_ar['datetimes'] = pd.to_datetime(df_ar['Date'], format='%d/%m/%Y')
    all_x=[]
    all_y=[]
    colors = ['k','saddlebrown','brown','r','sandybrown','darkkhaki','lawngreen','mediumspringgreen','mediumturquoise','royalblue','b','mediumblue']
    labels = ['20101212','20101214','20110317','20110605','20130123','20130129','20130209','20130424','20130502','20130517','20130527','20130608']
    for d in df_ar['datetimes']:
        d_str = d.strftime('%Y%m%d')
        x = np.array((df_ar.loc[df_ar['datetimes'] == d, 'gcs_radial_vel_at2']))

        if len(x) > 1:
            x = np.nanmean(x)
        x = float(x)

        y  = np.mean(df_arcades.loc[df_arcades['date'].dt.date == d.date(), 'length fit vel [km/s]'])
        all_x.append(x)
        all_y.append(y)

    all_x = all_x [0:11]
    all_y = all_y [0:11]

    all_x_filtered = drop_by_index(all_x,list_rejected_cmes_id)
    all_y_filtered = drop_by_index(all_y,list_rejected_cmes_id)
    colors_filtered = drop_by_index(colors,list_rejected_cmes_id)
    labels_filtered = drop_by_index(labels,list_rejected_cmes_id)
    #pearson coef and p-value
    linear_regresion = scipy.stats.linregress(all_x_filtered, all_y_filtered)
    slope = linear_regresion.slope
    intercept = linear_regresion.intercept
    pearson = linear_regresion.rvalue
    r_square = pearson*pearson
    p_value = linear_regresion.pvalue
    stdev = linear_regresion.stderr

    fig,ax = plt.subplots()
    #ejes y titulo
    ax.set_xlabel('GCS radial speed at 2 R$_{\odot}$ [km s$^{-1}$]', fontsize=14)
    ax.set_ylabel('PEA propagation speed [km s$^{-1}$]', fontsize=14)
    ax.set_title('', fontsize=18)
    for contador in range(len(all_x_filtered)):
        ax.scatter(all_x_filtered[contador],all_y_filtered[contador],c=colors_filtered[contador],label=labels_filtered[contador],alpha=0.9)

    if ajuste_lineal2:
        x_axis = np.linspace(np.min(all_x_filtered),np.max(all_x_filtered),10)
        ax.plot(x_axis, intercept + slope*x_axis, 'r', label=f'$r2 = {r_square:.2f}$')
    ax.legend(loc="best", prop={'size': 8})
    ax.grid(True)
    plt.savefig(hpath+'/pea_prop_speed_vs_gcs_radial_vel_at2.png')
    plt.close()   
