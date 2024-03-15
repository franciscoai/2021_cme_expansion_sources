#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computes length [km], vel of expansion and tilt angle [deg] from several arcades footpoints selected manually.
Writes a new table

@author: iglesias
"""
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunpy.coordinates.utils import GreatArc
#import sunpy.map
#from sunpy.data.sample import AIA_171_IMAGE
import pandas as pd
import datetime as dt
#from sunpy.coordinates import frames

local_path = os.getcwd()
input_dir = local_path + '/output_data/arcades'
odir = local_path + '/output_data/arcades/props'
arcsec2km = 725.27 # 1 arcsec = 725.27 km

##########
# get .csv files from input_dir
files = [os.path.join(input_dir,f) for f in os.listdir(input_dir) if f.endswith('.csv')]
print('Files found:', files)
# process each file and save the output
for f in files:
    df = pd.read_csv(f, header=0, delimiter=',')
    df['date'] = pd.to_datetime(df['file'], format='AIA%Y%m%d_%H%M%S_0193.fits')
    # delete rows with '[]' in colum 'lon [arcsec]'
    df = df[df['lon [arcsec]'] != '[]']
    # reset id to start from 0
    df.reset_index(drop=True, inplace=True)
    # convert the following columns from string to np.array of floats
    df['lon [arcsec]'] = df['lon [arcsec]'].apply(lambda x: np.array(x[1:-1].split(',')).astype(float))
    df['lat [arcsec]'] = df['lat [arcsec]'].apply(lambda x: np.array(x[1:-1].split(',')).astype(float))

    # for each row, separeate the footpoints in two halves and saves then in a new column for each footpoint
    df['lon1 [arcsec]'] = df['lon [arcsec]'].apply(lambda x: x[:int(len(x)/2)])
    df['lon2 [arcsec]'] = df['lon [arcsec]'].apply(lambda x: x[int(len(x)/2):])
    df['lat1 [arcsec]'] = df['lat [arcsec]'].apply(lambda x: x[:int(len(x)/2)])
    df['lat2 [arcsec]'] = df['lat [arcsec]'].apply(lambda x: x[int(len(x)/2):])
    
    # get sky coord
    df['fp1'] = df.apply(lambda x: [SkyCoord(x['lon1 [arcsec]'][i], x['lat1 [arcsec]'][i], unit='arcsec',frame="heliographic_carrington", observer="earth",  obstime=x['date']) for i in range(len(x['lon1 [arcsec]']))], axis=1)
    df['fp2'] = df.apply(lambda x: [SkyCoord(x['lon2 [arcsec]'][i], x['lat2 [arcsec]'][i], unit='arcsec',frame="heliographic_carrington", observer="earth",  obstime=x['date']) for i in range(len(x['lon2 [arcsec]']))], axis=1)

    # for each row, calculate the great arc distance distance between each footpoint (arcade width) and saves it in an array in a new column
    df['width [arcsec]'] = df.apply(lambda x: np.array([float(GreatArc(x['fp1'][i], x['fp2'][i]).distances()[1].value) for i in range(len(x['fp1']))]), axis=1)

    # for each row and each footpoint group, computes the total distance by adding the euclidean distance between conscutive footpoints
    df['length1 [arcsec]'] = df.apply(lambda x: np.array([np.sqrt((x['lon1 [arcsec]'][i]-x['lon1 [arcsec]'][i+1])**2 + (x['lat1 [arcsec]'][i]-x['lat1 [arcsec]'][i+1])**2) for i in range(len(x['lon1 [arcsec]'])-1)]), axis=1)
    df['length2 [arcsec]'] = df.apply(lambda x: np.array([np.sqrt((x['lon2 [arcsec]'][i]-x['lon2 [arcsec]'][i+1])**2 + (x['lat2 [arcsec]'][i]-x['lat2 [arcsec]'][i+1])**2) for i in range(len(x['lon2 [arcsec]'])-1)]), axis=1)

    # gets the mean of the two arcade lengths for each row
    df['length [arcsec]'] = df.apply(lambda x: np.mean([np.sum(x['length1 [arcsec]']), np.sum(x['length2 [arcsec]'])]), axis=1)
    df['length std [arcsec]'] = df.apply(lambda x: np.std([np.sum(x['length1 [arcsec]']), np.sum(x['length2 [arcsec]'])]), axis=1)

    # compute the arcade tilt, by fiting a line to the footpoints and getting the angle of the line
    # for each row, fit a line to the footpoints and get the angle of the line
    df['tilt [deg]'] = np.nan
    # for i in range(len(df)):
    #     # get the footpoints
    #     fp1 = SkyCoord(df['lon1 [arcsec]'][i], df['lat1 [arcsec]'][i], unit='arcsec',frame="heliographic_carrington", observer="earth",  obstime=df['date'][i])
    #     fp2 = SkyCoord(df['lon2 [arcsec]'][i], df['lat2 [arcsec]'][i], unit='arcsec',frame="heliographic_carrington", observer="earth",  obstime=df['date'][i])
    #     # get the line that fits the footpoints
    #     line = GreatArc(fp1, fp2)
    #     # get the angle of the line
    #     df['tilt [deg]'][i] = line.angle.deg

    ### PLOTS

    # creates odir/every file name
    odir_ev = os.path.join(odir, f.split('/')[-1].replace('.csv', ''))
    os.makedirs(odir_ev, exist_ok=True)

    # scatter plot of the footpoints, each group with a diff color
    for i in range(len(df)):
        plt.figure(figsize=(12,7))
        plt.scatter(df['lon1 [arcsec]'][i], df['lat1 [arcsec]'][i], label=df['date'][i], color='red')
        plt.scatter(df['lon2 [arcsec]'][i], df['lat2 [arcsec]'][i], label=df['date'][i], color='blue')
        plt.gca().set_aspect('equal')
        plt.xlabel('Longitude [arcsec]')
        plt.ylabel('Latitude [arcsec]')
        plt.title('Arcades footpoints')
        if i == 0:
            all_lon_flatten = np.array([])
            all_lat_flatten = np.array([])
            for j in range(len(df)):
                all_lon_flatten = np.append(all_lon_flatten, df['lon [arcsec]'][j])
                all_lat_flatten = np.append(all_lat_flatten, df['lat [arcsec]'][j])
            xlim = [np.min(all_lon_flatten), np.max(all_lon_flatten)]
            ylim = [np.min(all_lat_flatten), np.max(all_lat_flatten)]
            plt.xlim(xlim)
            plt.ylim(ylim)
        else:
            plt.xlim(xlim)
            plt.ylim(ylim)
        
        oimage = os.path.join(odir_ev, f.split('/')[-1].replace('.csv', '_footpoints'+str(i)+'.png'))
        plt.savefig(oimage)
        plt.close()

    # plot the mean and stddev (error bar) of arcade width vs date
    plt.figure(figsize=(12,7))
    plt.errorbar(df['date'], df['width [arcsec]'].apply(np.mean), yerr=df['width [arcsec]'].apply(np.std), fmt='o', label='measured')
    # fits a line to the data and add the slope 
    # date to sec
    dates_sec = (df['date'] - dt.datetime(1995,1,1)).dt.total_seconds()
    z = np.polyfit(dates_sec, df['width [arcsec]'].apply(np.mean), 1)
    p = np.poly1d(z)
    # Velocity from arcsec/date to km/s.
    vel = z[0]*arcsec2km
    plt.plot(df['date'], p(dates_sec), label=f'fitted: {vel:.2f} km/s')
    plt.xlabel('Date')
    plt.ylabel('Width [arcsec]')
    plt.title('Mean and stddev of arcade width vs date')
    plt.xticks(rotation=45)
    plt.legend()
    oimage = os.path.join(odir_ev, f.split('/')[-1].replace('.csv', '_width_vs_date.png'))
    plt.savefig(oimage)
    plt.close()    

    # plot the mean and stddev (error bar) of arcade length vs date
    plt.figure(figsize=(12,7))
    plt.errorbar(df['date'], df['length [arcsec]'], yerr=df['length std [arcsec]'], fmt='o', label='measured')
    z = np.polyfit(dates_sec, df['length [arcsec]'], 1)
    p = np.poly1d(z)
    # Velocity from arcsec/date to km/s.
    vel = z[0]*arcsec2km
    plt.plot(df['date'], p(dates_sec), label=f'fitted: {vel:.2f} km/s')
    plt.xlabel('Date')
    plt.ylabel('Length [arcsec]')
    plt.title('Mean and stddev of arcade length vs date')
    plt.xticks(rotation=45)
    plt.legend()
    oimage = os.path.join(odir_ev, f.split('/')[-1].replace('.csv', '_length_vs_date.png'))
    plt.savefig(oimage)
    plt.close()

    # plot the tilt angle vs date
    plt.figure(figsize=(12,7))
    plt.scatter(df['date'], df['tilt [deg]'])
    plt.xlabel('Date')
    plt.ylabel('Tilt [deg]')
    plt.title('Tilt angle vs date')
    plt.xticks(rotation=45)
    oimage = os.path.join(odir_ev, f.split('/')[-1].replace('.csv', '_tilt_vs_date.png'))
    plt.savefig(oimage)
    plt.close()

print('Done :-)')