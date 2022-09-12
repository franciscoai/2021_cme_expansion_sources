#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computes filament length and angle from two manual points measured

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
import sunpy.map
from sunpy.data.sample import AIA_171_IMAGE
import pandas as pd
import datetime as dt
from sunpy.coordinates import frames

local_path = os.getcwd()
repo_path = os.path.dirname(os.path.dirname(local_path))
file_path = repo_path + '/input_data/filament_coords.csv'


##########

df = pd.read_csv(file_path, header=0, delimiter=',')
df['date'] = pd.to_datetime(df['date'],format='%Y%m%d_%H:%M:%S')

p0 = SkyCoord(df['lon1']*u.deg, df['lat1']*u.deg,
              frame="heliographic_carrington",
              observer="earth",
              obstime=df['date'])  # "2010/01/01T00:00:30")
p1 = SkyCoord(df['lon2']*u.deg, df['lat2']*u.deg,
              frame="heliographic_carrington",
              observer="earth",
              obstime=df['date'])  # "2010/01/01T00:00:30")
#get arc lenght
arc=[]
for i in range(len(p0)):
    arc.append(GreatArc(p0[i], p1[i], points=2))


df['lenght'] = np.array([i.distances()[1].value for i in arc]).flatten()
print('******Lenghts:', df['lenght'] )
df['tilt_angle'] = None # Compute  properly
print('******Tilt Angles:', df['tilt_angle'])
df['inner_angle'] = np.rad2deg(np.array([i.inner_angles()[1].value for i in arc]).flatten())
print('******Inner Angles:', df['inner_angle'] )
df['tilt'] = np.rad2deg(np.array([i.inner_angles()[1].value for i in arc]).flatten())
print('******Tilt:', df['tilt'] )

#plots
arc=[]
for i in range(len(p0)):
    arc.append(GreatArc(p0[i], p1[i], points=20))

for i in range(len(arc)):
    plt.scatter(arc[i].coordinates().lon,arc[i].coordinates().lat, color='r', s=2)
plt.scatter([df['lon1'],df['lon2']],[df['lat1'],df['lat2']], color='b')
plt.show()

# m = sunpy.map.Map(AIA_171_IMAGE) 
# ax = plt.subplot(projection=m)  
# m.plot(axes=ax)  
# for i in range(len(arc)):
#     ax.plot_coord(arc[i].coordinates().transform_to(frames.Helioprojective), color='r')
# plt.show()
