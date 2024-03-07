# class that takes an input list of .fits files paths  and for each file, one at a time:
# 1- Read file as Sunpy MAPS object and plots it on screen using Sunpy built in functions (depending on the instrument specified in fits header)
# 2- Allows the user to click on the image to select points
# 3- Converts the selected pixel values to Carrington coordinates using the WCS information in the .fits file and Sunpy built in functions
# 4- Saves the selected points to a .csv file

import sunpy.map
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.colors as colors
import os

class SelectImgPoints:
    def __init__(self, fits_files, output_file,diff=False, coord_type='heliographic_stonyhurst'): # heliographic_carrington or heliographic_stonyhurst
        self.fits_files = fits_files
        self.output_file = output_file
        self.points = []
        self.diff = diff
        self.coord_type = coord_type
        # create odirectory if it does not exist
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

    def select_points(self):
        if self.diff:
        # if diff, then reas two consecutive files and computes the difference before plotting and selecting points
            for i in range(len(self.fits_files)-1):
                # Read the .fits file as a Sunpy MAPS object
                map_seq = sunpy.map.Map([self.fits_files[i], self.fits_files[i+1]], sequence=True)
                # Plot the image
                map = sunpy.map.Map(map_seq[0].quantity - map_seq[1].quantity, map_seq[0].meta)
                map.plot(norm=colors.Normalize(vmin=-200, vmax=200), cmap='Greys_r')
                # Allow the user to click on the image to select points
                points = plt.ginput(n=-1, timeout=0)
                # Convert the selected pixel values to Carrington coordinates using the WCS information in the .fits file and Sunpy built in functions
                wcs = map.wcs
                wcs_points= [wcs.pixel_to_world(i[0], i[1]) for i in points]
                carrington_points = [SkyCoord(i.Tx, i.Ty, frame=self.coord_type, obstime=map.date) for i in wcs_points]
                carrington_lon = [i.lon.arcsec for i in carrington_points]
                carrington_lat = [i.lat.arcsec for i in carrington_points]
                # Save the selected points to a common list including the file basename
                self.points.append([self.fits_files[i].split('/')[-1], carrington_lon, carrington_lat])
                plt.close()
                # prompts the user if they want to continue selecting points or stop
                if i < len(self.fits_files)-2:
                    cont = input('Do you want to continue selecting points? (y/n): ')
                    if cont.lower() == 'n':
                        break
            # Save the selected points to a .csv file   
            df = pd.DataFrame(self.points, columns=['file', 'lon [arcsec]', 'lat [arcsec]'])
            df.to_csv(self.output_file, index=False)
        else:
            for fits_file in self.fits_files:
                # Read the .fits file as a Sunpy MAPS object
                map = sunpy.map.Map(fits_file)
                # Plot the image
                map.plot()
                # Allow the user to click on the image to select points
                points = plt.ginput(n=-1, timeout=0)
                # Convert the selected pixel values to Carrington coordinates using the WCS information in the .fits file and Sunpy built in functions
                wcs = map.wcs
                wcs_points= [wcs.pixel_to_world(i[0], i[1]) for i in points]
                carrington_points = [SkyCoord(i.Tx, i.Ty, frame=self.coord_type, obstime=map.date) for i in wcs_points]
                carrington_lon = [i.lon.arcsec for i in carrington_points]
                carrington_lat = [i.lat.arcsec for i in carrington_points]
                # Save the selected points to a common list including the file basename
                self.points.append([fits_file.split('/')[-1], carrington_lon, carrington_lat])
                plt.close()
                # prompts the user if they want to continue selecting points or stop
                if i < len(self.fits_files)-2:
                    cont = input('Do you want to continue selecting points? (y/n): ')
                    if cont.lower() == 'n':
                        break                
            # Save the selected points to a .csv file   
            df = pd.DataFrame(self.points, columns=['file', 'lon [arcsec]', 'lat [arcsec]'])
            df.to_csv(self.output_file, index=False)


