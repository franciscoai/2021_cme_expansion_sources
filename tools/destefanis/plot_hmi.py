from tkinter.messagebox import NO
import os
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.data.sample
import sunpy.map
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from sunpy.coordinates import frames

def save_img(ifile,opath,point=None, point_px=None):
    """
    Gerardo Jose Destefanis (g.destefanis@alumno.um.edu.ar), 2022.04.05

    From a fits file generates a HMI map. Also, plots a map and saves the image by png file.
    Lets, given a point in arcseconds, knows the correspondent pair of pixels in the image. 
    Can realice the inverse operation, is say, given a point joined in pixels, knows it coordinates in arcseconds.

    INPUTS:

    ifile: path of the hmi fits file. 

    opath: path where the image is saved. 

    OPTIONAL INPUTS: 

    point: Optional point to plot in the image [x,y] in arcsecs. Default None

    point_px: Optional point of the image [x,y] in pixels. Default None 

    OUTPUTS:

    The plot of the image is save in path. 
    
    Coordinates of the points joined (in arcseconds and/or pixels). 

    """
    
    hmi_map= sunpy.map.Map(ifile)
    hmi_map.plot_settings['cmap'] = "hmimag"
    hmi_map.plot_settings['norm'] = plt.Normalize(-1500, 1500)
    fig = plt.figure(figsize=(6, 5))
    ax2 = fig.add_subplot(1, 1, 1, projection=hmi_map)
    hmi_map.plot(axes=ax2)

    #Print pixel coordinates for a given point given in arcsec units
    #Also plots point in the final image
    if point is not None:
        #Disclaimer:  Chequear que sistema de coordenadas estamos usando en el modulo frames hay varios definidos.
        point_coord = hmi_map.world_to_pixel(SkyCoord(point[0],point[1], frame = frames.HeliographicStonyhurst, unit = "arcsec"))
        print("The point coord in pixels units are:", point_coord.x,point_coord.y)

        map_coord = (point * u.arcsec)
        ax2.plot(map_coord[0].to('deg'), map_coord[1].to('deg'), 'o', color='white', transform=ax2.get_transform('world'))
        print(f'Map coordinate [{map_coord[0]}, {map_coord[1]}]')

    #Print pixel coordinates for a given point in pixels units
    if point_px is not None:
        point_px_coord = hmi_map.pixel_to_world(point_px[0]*u.pixel,point_px[1]*u.pixel)
        print("The pixel coord are:", point_px_coord.sistema_coordenadas.lat,point_px_coord.sistema_coordenadas.lon)

    plt.savefig(opath + "/" + os.path.basename(ifile) +".png")
    return



