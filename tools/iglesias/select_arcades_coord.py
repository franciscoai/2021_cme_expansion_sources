# Script to read the arcades initial and final times from the main .csv file and iteratively read images in that time range and allow the user to select points on the images
# using select_img_points.py

import pandas as pd
import os
from select_img_points import SelectImgPoints
import numpy as np

# Constants
# Event to process
id = 1
# cadence of the differential images in seconds
cadence = 60.*15
# Path to the main .csv file
csv= os.getcwd() + '/input_data/ar.csv'
databse= '/gehme/data'
# odir in /output_data of the current repository
odir = os.getcwd() + '/output_data/arcades'
aia_instrument_path = '/gehme/data/sdo/aia/l1/193'
# file structure within aia is:
# aia_instrument_path + 'YYMMDD/file
euvia_instrument_path = '/gehme/data/stereo/secchi/L0/a/img/euvi/'
euvib_instrument_path = '/gehme/data/stereo/secchi/L0/b/img/euvi/'
#file structure within euvia and euvib are: 
#euvia_instrument_path + 'YYYMMDD/preped/file
#euvib_instrument_path + 'YYYMMDD/preped/file

# Functions
# Read the main .csv file
def read_main_csv(file):
    df = pd.read_csv(file)
    return df

# Selects only the event id specified by the user.
def select_event(df, id):
    event = df[df['ID'] == id]
    return event

# Get the "EUV loops data start time" and "EUV loops data end time" times of the selected event in format 'YYYY/MM/DD HH:MM:SS' and returns them in datetime objects
def get_event_times(event):
    start_time = pd.to_datetime(event['EUV loops data start time'].values[0])
    end_time = pd.to_datetime(event['EUV loops data end time'].values[0])
    return start_time, end_time

# for a given file name and instrument, returns the datetime object extracted from the file name.
# formats are different for each instrument
# AIA: aia.lev1.193A_2010-12-12T01_35_07.84Z.image_lev1.fits
# EUVI-A: YYYYMMDD_HHmmSS_14euA.fts
# EUVI-B: YYYYMMDD_HHmmSS_14euB.fts
def get_time_from_file(file, instrument):
    name = file.split('/')[-1]
    if instrument == 'AIA':
        time = name.split('193A_')[1].split('.')[0]
        time = pd.to_datetime(time, format='%Y-%m-%dT%H_%M_%S')
    elif instrument == 'EUVI-A':
        time = name.split('.')[0].split('_14e')[0]
        time = pd.to_datetime(time, format='%Y%m%d_%H%M%S')
    elif instrument == 'EUVI-B':
        time = name.split('.')[0].split('_14e')[0]
        time = pd.to_datetime(time, format='%Y%m%d_%H%M%S')
    return time

# Get all the .fits files that have a datetime within start_time and end_time from the database of the instrument specified
def get_fits_files_paths(start_time, end_time, instrument, database):
    # get YYYYMMDD from start_time
    current_day = start_time.strftime('%Y%m%d')
    if instrument == 'AIA':
        instrument_path = aia_instrument_path
        instrument_path = os.path.join(instrument_path, current_day)        
    elif instrument == 'EUVI-A':
        instrument_path = euvia_instrument_path
        instrument_path = os.path.join(instrument_path, current_day,'preped')   
    elif instrument == 'EUVI-B':
        instrument_path = euvib_instrument_path
        instrument_path = os.path.join(instrument_path, current_day,'preped')            
    # Get all the .fits files that have a datetime within start_time and end_time from the database of the instrument specified
    files = [os.path.join(instrument_path, file) for file in os.listdir(instrument_path) if file.endswith('.fits') or file.endswith('.fts')]
    files = [file for file in files if start_time <= get_time_from_file(file, instrument) <= end_time]
    return sorted(files)

def filter_fits_files(fits_files, cadence):
    # get the time of all files
    times = [get_time_from_file(file, instrument) for file in fits_files]
    # get the time difference between consecutive files
    time_diff = [j-i for i, j in zip(times[:-1], times[1:])]
    # get the time difference in seconds
    time_diff = [i.total_seconds() for i in time_diff]
    # starting from the first, select only files that are at least cadence seconds apart. If the jump is 1.5*cadence print a warning
    new_files = [fits_files[0]]
    flag=0
    for i, t in enumerate(time_diff):
        if flag==0:
            acc_t = t
        else:
            acc_t += t
        if acc_t >= cadence:
            new_files.append(fits_files[i+1])
            flag=0
        elif acc_t >= 1.5*cadence:
            print(f'Warning: time difference between {fits_files[i]} and {fits_files[i+1]} is {t} seconds')
        else:
            flag=1
    return new_files
# main
# Read the main .csv file
df = read_main_csv(csv)
# Select only the event id specified by the user
event = select_event(df, id)
# Get the "EUV loops data start time" and "EUV loops data end time" times of the selected event
start_time, end_time = get_event_times(event)
# Get the .fits files paths of the selected event from the database of the instrument specified in "Instrument" column
instrument = event['Instrument'].values[0]
fits_files = get_fits_files_paths(start_time, end_time, instrument, databse)
if len(fits_files) == 0:
    print(f'Error: No fits files found for event {id} in the specified time range')
    os._exit(0)
if cadence > 0:
    fits_files = filter_fits_files(fits_files, cadence)
if len(fits_files) == 0:
    print(f'Error: No fits files found for event {id} in the specified time range with the specified cadence')
    os._exit(0)
# create the output file
output_file = os.path.join(odir, f'event_{id}_points.csv')
# get roi from df column ROI box [arcsec], expressed in arcsec. Format is xmin/xmax/ymin/ymax
roi = event['ROI box [arcsec]'].values[0]
if type(roi) == str:
    roi = [int(i) for i in roi.split('/')]
else:
    roi = None
# Create an instance of the SelectImgPoints that selects points on the images and saves them to a .csv file
select_img_points = SelectImgPoints(fits_files, output_file, diff=True, roi=roi)
# Select points
select_img_points.select_points()
# Print the path of the output file
print(f'Points saved to {output_file}')






