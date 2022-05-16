import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import pandas as pd
import lu_pil as lp
import ipdb
import astropy.units as u
from astropy.coordinates import SkyCoord

sharp = sunpy.map.Map("/gehme/data/sdo/hmi/sharps/20101214/*Br.fits", sequence=True)

shtest = sharp[0]

#Testing a submap
#shtest = shtest.submap(SkyCoord(170 * u.deg, 10 * u.deg, frame=shtest.coordinate_frame),top_right=SkyCoord(185 * u.deg, 20 * u.deg, frame=shtest.coordinate_frame))

coords = lp.magpol_barycenters(shtest)
print(coords)

# Seeing where the computed barycenter is
fig = plt.figure(frameon=False)
ax = plt.axes([0, 0, 1, 1])

ax.set_axis_off()
norm = shtest.plot_settings["norm"]
ax.imshow(shtest.data, norm=norm, cmap=shtest.plot_settings["cmap"], origin="lower")
negative_pol = plt.Circle(coords["Neg"], 6, color="blue", fill=False)
positive_pol = plt.Circle(coords["Pos"], 6, color="red", fill=False)
ax.add_artist(negative_pol)
ax.add_artist(positive_pol)
#plt.show()

pil_test = lp.pil_computation(shtest)

ipdb.set_trace()