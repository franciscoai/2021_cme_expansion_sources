from tkinter.messagebox import NO
import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.data.sample
import sunpy.map
from astropy.io import fits
from astropy.coordinates import SkyCoord

def save_img(ifile,opath,point=None):
    """
    
    :point: Optional point to plot in the image [x,y] in arcsecs. Default None
    """
    
    hmi_map= sunpy.map.Map(ifile)
    hmi_map.plot_settings['cmap'] = "hmimag"
    hmi_map.plot_settings['norm'] = plt.Normalize(-1500, 1500)
    fig = plt.figure(figsize=(6, 5))
    ax2 = fig.add_subplot(1, 1, 1, projection=hmi_map)
    hmi_map.plot(axes=ax2)
    if point is not None: 
        map_coord = (point * u.arcsec)
        ax2.plot(map_coord[0].to('deg'), map_coord[1].to('deg'), 'o', color='white', transform=ax2.get_transform('world'))     
        print(f'Map coordinate [{map_coord[0]}, {map_coord[1]}]')
    plt.savefig(opath + "/" + os.path.basename(ifile) +".png")
    return



