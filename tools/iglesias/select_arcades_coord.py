# Script to read the arcades initial and final times from the main .csv file and iteratively read images in that time range and allow the user to select points on the images
# using select_img_points.py


from select_img_points import SelectImgPoints
import pandas as pd
import os
import numpy as np


# Constants
# Event to process
id = 4
overwrite = False # if True, the output file will be overwritten if it already exists
# time difference of the differential images in seconds
img_time_diff = 60.*30.
# minima cadence of the differential images in seconds, use None to keep all the images
cadence = 60.*20.
# Amount of sigmas to control image color scale range
color_scl = 3
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
        time = name.split('AIA')[1].split('_0193')[0]
        time = pd.to_datetime(time, format='%Y%m%d_%H%M%S')
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
        instrument_path = os.path.join(instrument_path, current_day,'preped')    
    elif instrument == 'EUVI-A':
        instrument_path = euvia_instrument_path
        instrument_path = os.path.join(instrument_path, current_day,'preped')   
    elif instrument == 'EUVI-B':
        instrument_path = euvib_instrument_path
        instrument_path = os.path.join(instrument_path, current_day,'preped')   
    # get YYYYMMDD from end_time          
    last_day = end_time.strftime('%Y%m%d')
    # if start_time and end_time are in different days, get all the files from the start day and the end day
    if current_day != last_day:
        # if current_day does not exists, skips it
        if os.path.exists(instrument_path):
            files = [os.path.join(instrument_path, file) for file in os.listdir(instrument_path) if file.endswith('.fits') or file.endswith('.fts')]
            files = [file for file in files if start_time <= get_time_from_file(file, instrument) <= end_time]
        else:
            files = []
            print(f'Warning: No files found for {instrument} in {current_day}')
        # get all the files from the last day if it exists
        instrument_path = instrument_path.replace(current_day, last_day)
        if os.path.exists(instrument_path):
            files += [os.path.join(instrument_path, file) for file in os.listdir(instrument_path) if file.endswith('.fits') or file.endswith('.fts')]
            files = [file for file in files if start_time <= get_time_from_file(file, instrument) <= end_time]
        else:
            print(f'Warning: No files found for {instrument} in {last_day}')
        return sorted(files)
    else:
        files = [os.path.join(instrument_path, file) for file in os.listdir(instrument_path) if file.endswith('.fits') or file.endswith('.fts')]
        files = [file for file in files if start_time <= get_time_from_file(file, instrument) <= end_time]
        return sorted(files)
    
def filter_fits_files(fits_files, img_time_diff, cadance=None):
    # get the time of all files
    times = [get_time_from_file(file, instrument) for file in fits_files]
    # get the time difference between consecutive files
    time_diff = [j-i for i, j in zip(times[:-1], times[1:])]
    # get the time difference in seconds
    time_diff = [i.total_seconds() for i in time_diff]
    # For each element in fits_files returns two elements, the element itself and the next file that has 
    # an accumulated time difference greater than img_time_diff
    new_files = []
    acc_time_diff = 0
    for i in range(len(time_diff)):
        new_files.append(fits_files[i])
        acc_time_diff += time_diff[i]
        for j in range(i+1, len(time_diff)):
            if acc_time_diff >= img_time_diff:
                new_files.append(fits_files[j])
                acc_time_diff = 0
                break
            else:
                acc_time_diff += time_diff[j]          

    # keeps only pairs of consecutive files with cadence greater than cadence
    final_files = new_files.copy()
    if cadance is not None:
        # computes the new time_diff
        times = [get_time_from_file(file, instrument) for file in new_files[0::2]]  
        new_time_diff = [j-i for i, j in zip(times[:-1], times[1:])]
        new_time_diff = [i.total_seconds() for i in new_time_diff]
        # keeps only pairs of consecutive files with cadence greater than cadence
        final_files = [new_files[0],new_files[1]]
        acc_new_time_diff=new_time_diff[0]
        for i in range(1,len(new_time_diff)-1):
            if  acc_new_time_diff >= cadance:
                final_files.append(new_files[2*i])
                final_files.append(new_files[2*i+1])
                acc_new_time_diff=new_time_diff[i]
            else:
                acc_new_time_diff+=new_time_diff[i]
    #print side by side the filenames only of the pairs
    print('Selected ' + str(len(final_files)//2) + ' file pairs')
    for i in range(0,len(final_files),2):
        print(final_files[i].split('/')[-1], final_files[i+1].split('/')[-1])
    return final_files
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
if img_time_diff > 0:
    fits_files = filter_fits_files(fits_files, img_time_diff, cadance=cadence)
if len(fits_files) == 0:
    print(f'Error: No fits files found for event {id} in the specified time range with the specified img_time_diff')
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
select_img_points = SelectImgPoints(fits_files, output_file, diff=True, roi=roi, overwrite=overwrite, color_scl=color_scl)
# Select points
select_img_points.select_points()
# Print the path of the output file
print(f'Points saved to {output_file}')






