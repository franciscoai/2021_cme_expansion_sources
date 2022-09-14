#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compares the properties of various fenomena associated to a CME:
- Magnetic AR morphological and magnetic properties
- Ca filament morphological properties
- White light CME 3D morphological properties from GCS forward modeling

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
df['date'] = pd.to_datetime(df['date'], format='%Y%m%d_%H:%M:%S')

p0 = SkyCoord(df['lon1']*u.deg, df['lat1']*u.deg,
              frame="heliographic_carrington",
              observer="earth",
              obstime=df['date'])
p1 = SkyCoord(df['lon2']*u.deg, df['lat2']*u.deg,
              frame="heliographic_carrington",
              observer="earth",
              obstime=df['date'])
p2 = SkyCoord(df['lon2']*u.deg, df['lat1']*u.deg,  # third point of the spherical triangle
              frame="heliographic_carrington",
              observer="earth",
              obstime=df['date'])

# get lenght of the triangle sides in km
a = []
b = []
c = []  # filaments lenght
for i in range(len(p0)):
    a.append(GreatArc(p1[i], p2[i], points=3))
    b.append(GreatArc(p0[i], p2[i], points=3))
    c.append(GreatArc(p0[i], p1[i], points=3))
df['lenght'] = np.array([float(i.distances()[2].value) for i in c]).flatten()
print('******Filament lenght [km]:', df['lenght'])

# filament (c) tilt angles
tilt = []
for i in range(len(c)):
    #angular sides
    a_ang = float(a[i].inner_angles()[2].value)  
    b_ang = float(b[i].inner_angles()[2].value)  
    c_ang = float(c[i].inner_angles()[2].value) 
    ctilt = np.rad2deg(np.arccos((np.cos(a_ang) - np.cos(b_ang) * np.cos(c_ang)) / (np.sin(b_ang) * np.sin(c_ang))))
    if df['lat2'][i] > df['lat1'][i]:
        tilt.append(ctilt)
    else:
        tilt.append(-ctilt)

df['tilt_angle'] = tilt
print('******Tilt Angles:', df['tilt_angle'])

# plots
a = []
b = []
c = [] 
for i in range(len(p0)):
    a.append(GreatArc(p1[i], p2[i], points=20))
    b.append(GreatArc(p0[i], p2[i], points=20))
    c.append(GreatArc(p0[i], p1[i], points=20))

for i in range(len(c)):
    plt.scatter(a[i].coordinates().lon, a[i].coordinates().lat, color='r', s=2)
    plt.scatter(b[i].coordinates().lon, b[i].coordinates().lat, color='g', s=2)
    plt.scatter(c[i].coordinates().lon, c[i].coordinates().lat, color='b', s=2)
    plt.scatter([df['lon1'][i]], [df['lat1'][i]], color='y')
    plt.scatter([df['lon2'][i]], [df['lat2'][i]], color='k')
    plt.scatter([ df['lon2'][i]], [df['lat1'][i]], color='c')  
    plt.gca().set_aspect('equal')
plt.show()

# m = sunpy.map.Map(AIA_171_IMAGE)
# ax = plt.subplot(projection=m)
# m.plot(axes=ax)
# for i in range(len(c)):
#     ax.plot_coord(c[i].coordinates().transform_to(frames.Helioprojective), color='r')
# plt.show()
