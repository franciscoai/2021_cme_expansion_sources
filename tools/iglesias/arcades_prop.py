#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computes length [km], vel of expansion and tilt angle [deg] from several arcades footpoints selected manually.
Writes a new table

@author: iglesias

Output fields:

LengthX are the lengths of the consecutive segments formed by each group of selected points (the sum is the total len)


"""
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunpy.coordinates.utils import GreatArc
#import sunpy.map
#from sunpy.data.sample import AIA_171_IMAGE
import pandas as pd
import datetime as dt
from scipy.optimize import curve_fit
#from sunpy.coordinates import frames

def linear_func(y, a, b):
    return a*y + b

########## MAIN ##########
# constants
local_path = os.getcwd()
input_dir = local_path + '/output_data/arcades'
odir = local_path + '/output_data/arcades/props'
arcsec2km = 725.27 # 1 arcsec = 725.27 km

# get .csv files from input_dir
files = sorted([os.path.join(input_dir,f) for f in os.listdir(input_dir) if f.endswith('.csv')])
print('Files found:', files)
# process each file and save the output
all_df = pd.DataFrame() # to save all the dataframes
for f in files:
    df = pd.read_csv(f, header=0, delimiter=',')
    try:
        df['date'] = pd.to_datetime(df['file'], format='AIA%Y%m%d_%H%M%S_0193.fits')
    except:
        try:
            df['date'] = pd.to_datetime(df['file'], format='%Y%m%d_%H%M%S_14euA.fts')
        except:
            df['date'] = pd.to_datetime(df['file'], format='%Y%m%d_%H%M%S_14euB.fts')
    # delete rows with '[]' in colum 'lon [arcsec]'
    df = df[df['lon [arcsec]'] != '[]']
    # reset id to start from 0
    df.reset_index(drop=True, inplace=True)
    # convert the following columns from string to np.array of floats
    df['lon [arcsec]'] = df['lon [arcsec]'].apply(lambda x: np.array(x[1:-1].split(',')).astype(float))
    df['lat [arcsec]'] = df['lat [arcsec]'].apply(lambda x: np.array(x[1:-1].split(',')).astype(float))

    # for each row, if the number of points is even separeates the footpoints in two halves and saves each group in a new column
    # if the number of points is odd, separeates the footpoints in three groups and saves each group in a new column
    df['lon1 [arcsec]'] = df.apply(lambda x: x['lon [arcsec]'][:len(x['lon [arcsec]'])//2] if len(x['lon [arcsec]'])%2==0 else x['lon [arcsec]'][:len(x['lon [arcsec]'])//3], axis=1)
    df['lat1 [arcsec]'] = df.apply(lambda x: x['lat [arcsec]'][:len(x['lat [arcsec]'])//2] if len(x['lat [arcsec]'])%2==0 else x['lat [arcsec]'][:len(x['lat [arcsec]'])//3], axis=1)
    df['lon2 [arcsec]'] = df.apply(lambda x: x['lon [arcsec]'][len(x['lon [arcsec]'])//2:] if len(x['lon [arcsec]'])%2==0 else x['lon [arcsec]'][len(x['lon [arcsec]'])//3:2*len(x['lon [arcsec]'])//3], axis=1)
    df['lat2 [arcsec]'] = df.apply(lambda x: x['lat [arcsec]'][len(x['lat [arcsec]'])//2:] if len(x['lat [arcsec]'])%2==0 else x['lat [arcsec]'][len(x['lat [arcsec]'])//3:2*len(x['lat [arcsec]'])//3], axis=1)
    df['lon3 [arcsec]'] = df.apply(lambda x: x['lon [arcsec]'][len(x['lon [arcsec]'])//2:]*np.nan if len(x['lon [arcsec]'])%2==0 else x['lon [arcsec]'][2*len(x['lon [arcsec]'])//3:], axis=1)
    df['lat3 [arcsec]'] = df.apply(lambda x: x['lat [arcsec]'][len(x['lat [arcsec]'])//2:]*np.nan if len(x['lat [arcsec]'])%2==0 else x['lat [arcsec]'][2*len(x['lat [arcsec]'])//3:], axis=1) 

    # get sky coord
    df['fp1'] = df.apply(lambda x: [SkyCoord(x['lon1 [arcsec]'][i], x['lat1 [arcsec]'][i], unit='arcsec',frame="heliographic_stonyhurst", observer="earth",  obstime=x['date']) for i in range(len(x['lon1 [arcsec]']))], axis=1)
    df['fp2'] = df.apply(lambda x: [SkyCoord(x['lon2 [arcsec]'][i], x['lat2 [arcsec]'][i], unit='arcsec',frame="heliographic_stonyhurst", observer="earth",  obstime=x['date']) for i in range(len(x['lon2 [arcsec]']))], axis=1)
    df['fp3'] = df.apply(lambda x: [SkyCoord(x['lon3 [arcsec]'][i], x['lat3 [arcsec]'][i], unit='arcsec',frame="heliographic_stonyhurst", observer="earth",  obstime=x['date']) for i in range(len(x['lon3 [arcsec]']))], axis=1)

    # for each row, calculate the great arc distance distance between each footpoint (arcade width) and saves it in an array in a new column
    df['width [arcsec]'] = df.apply(lambda x: np.array([float(GreatArc(x['fp1'][i], x['fp2'][i]).distances()[-1].value) for i in range(len(x['fp1']))]), axis=1)

    # for each row and each footpoint group, computes the total distance by adding the euclidean distances between conscutive footpoints
    df['length1 [arcsec]'] = df.apply(lambda x: np.array([np.sqrt((x['lon1 [arcsec]'][i]-x['lon1 [arcsec]'][i+1])**2 + (x['lat1 [arcsec]'][i]-x['lat1 [arcsec]'][i+1])**2) for i in range(len(x['lon1 [arcsec]'])-1)]), axis=1)
    df['length2 [arcsec]'] = df.apply(lambda x: np.array([np.sqrt((x['lon2 [arcsec]'][i]-x['lon2 [arcsec]'][i+1])**2 + (x['lat2 [arcsec]'][i]-x['lat2 [arcsec]'][i+1])**2) for i in range(len(x['lon2 [arcsec]'])-1)]), axis=1)
    df['length3 [arcsec]'] = df.apply(lambda x: np.array([np.sqrt((x['lon3 [arcsec]'][i]-x['lon3 [arcsec]'][i+1])**2 + (x['lat3 [arcsec]'][i]-x['lat3 [arcsec]'][i+1])**2) for i in range(len(x['lon3 [arcsec]'])-1)]), axis=1)
    # adds df['length3 [arcsec]'] to df['length2 [arcsec]'] where df['length3 [arcsec]'] is not nan, otherwise adds 0
    df['length2+3 [arcsec]'] = df.apply(lambda x: np.array([x['length2 [arcsec]'][i] + x['length3 [arcsec]'][i] if not np.isnan(x['length3 [arcsec]'][i]) else x['length2 [arcsec]'][i] for i in range(len(x['length2 [arcsec]']))]), axis=1)

    # compute the arcade tilt angle by fitting a line to each group of points
    df['tilt mean [deg]'] = np.nan
    df['tilt sd [deg]'] = np.nan
    df['fp4'] = df.apply(lambda x: [SkyCoord(x['lon2 [arcsec]'][i], x['lat1 [arcsec]'][i], unit='arcsec',frame="heliographic_stonyhurst", observer="earth",  obstime=x['date']) for i in range(len(x['lon2 [arcsec]']))], axis=1)
    all_fits = []
    all_tilts = []
    for i in range(len(df)):
        p0= df['fp1'][i]
        # fits a line to lat and lon and gets the tilt from the line slope
        lat =  np.array([p.lat.arcsec for p in p0])
        lon =  np.array([p.lon.arcsec for p in p0])
        # check if the expected line is too vertical, if so switches lon and lat
        if np.abs(np.max(lon) - np.min(lon)) < np.abs(np.max(lat) - np.min(lat)):
            invert_fit = True
            z1 = curve_fit(linear_func, lat, lon)[0]
            tilt_p0 = np.arctan(1./z1[0])*180/np.pi            
        else:
            invert_fit = False
            z1 = curve_fit(linear_func, lon, lat)[0]
            tilt_p0 = np.arctan(z1[0])*180/np.pi              
        # same for p1
        p1= df['fp2'][i]
        lat =  np.array([p.lat.arcsec for p in p1])
        lon =  np.array([p.lon.arcsec for p in p1])
        if invert_fit:
            z2 = curve_fit(linear_func, lat, lon)[0]
            tilt_p1 = np.arctan(1./z2[0])*180/np.pi
        else:
            z2 = curve_fit(linear_func, lon, lat)[0]
            tilt_p1 = np.arctan(z2[0])*180/np.pi         
        # same for p2
        p2= df['fp4'][i]
        lat =  np.array([p.lat.arcsec for p in p2])
        lon =  np.array([p.lon.arcsec for p in p2])
        if invert_fit:
            z3 = curve_fit(linear_func, lat, lon)[0]
            tilt_p2 = np.arctan(1./z3[0])*180/np.pi
        else:
            z3 = curve_fit(linear_func, lon, lat)[0]
            tilt_p2 = np.arctan(z3[0])*180/np.pi
        #save all fits
        all_fits.append([z1, z2, z3, invert_fit])
        # computes the tilt angle as the mean angle of the two similar arcades
        # if the angles are mroe than 90 deg apart, adds 180 deg to the angles with the largest absolut value
        diff = abs(tilt_p0 - tilt_p1)
        print('**Date:', df['date'][i])
        print('tilt_p0', tilt_p0, 'tilt_p1', tilt_p1)
        if diff > 90:
            if np.abs(tilt_p0) > np.abs(tilt_p1):
                if tilt_p0 < 0:
                    tilt_p0 += 180
                else:
                    tilt_p0 -= 180
            else:
                if tilt_p1 < 0:
                    tilt_p1 += 180
                else:
                    tilt_p1 -= 180
        diff = np.abs(tilt_p0 - tilt_p1)
        mean_ang = np.mean([tilt_p0, tilt_p1])
        df['tilt mean [deg]'][i] = mean_ang
        df['tilt sd [deg]'][i] = diff       
        all_tilts.append([tilt_p0, tilt_p1, tilt_p2])        
    # corrects tilts that are >140 deg appart wrt median
    median_tilt = df['tilt mean [deg]'].median()
    for i in range(len(df)):
        if np.abs(df['tilt mean [deg]'][i] - median_tilt) > 140:
            if df['tilt mean [deg]'][i] < 0:
                df['tilt mean [deg]'][i] += 180
            else:
                df['tilt mean [deg]'][i] -= 180
    ### PLOTS

    # creates odir/every file name with event date at the end
    odir_ev = os.path.join(odir, f.split('/')[-1].replace('.csv', '')+'_'+str(df['date'][0].date()))
    os.makedirs(odir_ev, exist_ok=True)

    # scatter plot of the footpoints, each group with a diff color
    for i in range(len(df)):
        plt.figure(figsize=(12,7))
        for j in range(len(df['fp1'][i])):
            fp1= df['fp1'][i][j]
            fp2= df['fp2'][i][j]
            fp3= df['fp3'][i][j]
            plt.scatter(fp1.lon.arcsec, fp1.lat.arcsec, color='r')
            plt.scatter(fp2.lon.arcsec, fp2.lat.arcsec, color='b')
            plt.scatter(fp3.lon.arcsec, fp3.lat.arcsec, color='g')
            arc = GreatArc(fp1, fp2, points=20)
            arc_lon = arc.coordinates().lon.arcsec
            arc_lat = arc.coordinates().lat.arcsec
            plt.plot(arc_lon, arc_lat, color='k')
        plt.gca().set_aspect('equal')
        plt.xlabel('Longitude [arcsec]')
        plt.ylabel('Latitude [arcsec]')
        # overplots fit line for p0
        z1 = all_fits[i][0]
        inverted = all_fits[i][3]
        if inverted:
            lon = np.array([fp.lat.arcsec for fp in df['fp1'][i]])
            p = np.poly1d(z1)
            plt.plot(p(lon), lon, '--r')
        else:
            lon = np.array([fp.lon.arcsec for fp in df['fp1'][i]])
            p = np.poly1d(z1)
            plt.plot(lon, p(lon), '--r')                         
        # overplots fit line for p1
        z2 = all_fits[i][1]
        if inverted:
            lon = np.array([fp.lat.arcsec for fp in df['fp2'][i]])
            p = np.poly1d(z2)
            plt.plot(p(lon), lon, '--b')
        else:
            lon = np.array([fp.lon.arcsec for fp in df['fp2'][i]])
            p = np.poly1d(z2)
            plt.plot(lon, p(lon), '--b')   
        plt.title('Arcades footpoints for Date: ' + str(df['date'][i]) + 
                  '\n Angles: Red '+format(all_tilts[i][0], '.2f')+' Blue '+format(all_tilts[i][1], '.2f'))

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
    #saves fit results in pandas
    df['width fit par0'] = [z[0] for i in range(len(df))]
    df['width fit par1'] = [z[1] for i in range(len(df))]
    df['width fit vel [km/s]'] = [vel for i in range(len(df))]
    plt.plot(df['date'], p(dates_sec), label=f'fitted: {vel:.2f} km/s', color='k')
    plt.xlabel('Date')
    plt.ylabel('Width [arcsec]')
    plt.title('Mean and stddev of arcade width vs date')
    plt.xticks(rotation=45)
    plt.legend()
    oimage = os.path.join(odir_ev, f.split('/')[-1].replace('.csv', '_width_vs_date.png'))
    plt.savefig(oimage)
    plt.close()    

    # plot the arcade lengths of each group vs date
    plt.figure(figsize=(12,7))
    plt.scatter(df['date'], df['length1 [arcsec]'].apply(np.sum), color='r', label='group 1')
    plt.scatter(df['date'], df['length2 [arcsec]'].apply(np.sum) , color='b', label='group 2')
    plt.scatter(df['date'], df['length3 [arcsec]'].apply(np.sum) , color='g', label='group 3')
    plt.scatter(df['date'], df['length2+3 [arcsec]'].apply(np.sum) , color='k', label='group 2+3')
    # fits a line to df['length2+3 [arcsec]'].
    z = np.polyfit(dates_sec, df['length2+3 [arcsec]'].apply(np.sum), 1)
    p = np.poly1d(z)
    # Velocity from arcsec/date to km/s.
    vel = z[0]*arcsec2km
    df['length fit par0'] = [z[0] for i in range(len(df))]
    df['length fit par1'] = [z[1] for i in range(len(df))]
    df['length fit vel [km/s]'] = [vel for i in range(len(df))]
    plt.plot(df['date'], p(dates_sec), label=f'fitted: {vel:.2f} km/s', color='k')
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
    plt.errorbar(df['date'], df['tilt mean [deg]'], yerr=df['tilt sd [deg]'], fmt='o', label='measured')
    # fits a line and shows the line full equation in the plot label, use sci notaiton
    z = np.polyfit(dates_sec, df['tilt mean [deg]'], 1)
    p = np.poly1d(z)
    df['tilt fit par0'] = [z[0] for i in range(len(df))]
    df['tilt fit par1'] = [z[1] for i in range(len(df))]
    df['tilt fit mean value [deg]'] = [p(dates_sec).mean() for i in range(len(df))]
    plt.plot(df['date'], p(dates_sec), label=f'vel [deg/s]: {z[0]:.2e}', color='k')
    plt.xlabel('Date')
    plt.ylabel('Tilt [deg]')
    plt.title('Tilt angle vs date' + '\n Mean value: {:.3f} deg'.format(df['tilt mean [deg]'].mean()))
    plt.legend()
    plt.xticks(rotation=45)
    oimage = os.path.join(odir_ev, f.split('/')[-1].replace('.csv', '_tilt_vs_date.png'))
    plt.savefig(oimage)
    plt.close()
    #adds event id to df
    df['event'] = float(f.split('/')[-1].replace('.csv', '').split('_')[1])
    # appends to main df
    all_df = all_df.append(df, ignore_index=True)

# saves the main df in csv
all_df.to_csv(os.path.join(odir, 'all_arcades_props.csv'))

print('Done :-)')