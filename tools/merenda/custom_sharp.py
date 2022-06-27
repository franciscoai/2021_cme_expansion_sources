# import os
# from SHARPs import calculate_crop_sharpkeys as cs
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
import astropy.units as apu
from lu_tools import AreaSelector, cornerl_order

# Select files to
sharps_dir = "/gehme/data/sdo/hmi/sharps/20101214/"

br_files = glob(sharps_dir + "*Br.fits")
br_files.sort()
bt_files = glob(sharps_dir + "*Bt.fits")
bt_files.sort()
bp_files = glob(sharps_dir + "*Bp.fits")
bp_files.sort()
br_err_files = glob(sharps_dir + "*Br_err.fits")
br_err_files.sort()
bt_err_files = glob(sharps_dir + "*Bt_err.fits")
bt_err_files.sort()
bp_err_files = glob(sharps_dir + "*Bp_err.fits")
bp_err_files.sort()
confdis_files = glob(sharps_dir + "*conf_disambig.fits")
confdis_files.sort()
bitmap_files = glob(sharps_dir + "*bitmap.fits")
bitmap_files.sort()
los_files = glob(sharps_dir + "*magnetogram.fits")
los_files.sort()

# we have to select which subregion we need and then we can
# adapt bobras script to make sharps csv out of submap
# select area to crop
print("Choose the area to crop from the Br map")
preview = sunpy.map.Map(br_files[0])
preview.plot(title="Use mouse to select and press Q after selection is made")
areatoselect = AreaSelector(plt.gcf(), plt.gca())
plt.show()

# Give crop coordinates in bottomleft, topright order in order
# to use them with sunpy's map submap method.
final_coords = cornerl_order(areatoselect.data_coords,"bltp","np")
bottomleft = preview.pixel_to_world(final_coords[0][0] * apu.pix,
                                   final_coords[0][1] * apu.pix)
topright = preview.pixel_to_world(final_coords[1][0] * apu.pix,
                                  final_coords[1][1] * apu.pix)

# Make the submap # this is only for testing purposes.
crop = preview.submap(bottom_left=bottomleft, top_right=topright)
crop.peek()

print(f"Selected coords are: Bottom left corner:\n  Pix coords:  {final_coords[0]}\n"
      f"WCS coords:  {bottomleft},\nTop Right:\n  Pix coords:  {final_coords[1]}\n "
      f"WCS coords:  {topright}")








