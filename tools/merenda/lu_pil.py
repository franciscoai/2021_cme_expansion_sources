# Polarity Inversion Line computation
# Author: From Poisson et al 2014 paper coded by Luciano A. Merenda

from scipy.optimize import minimize
from scipy.integrate import dblquad
import numpy as np


def sign(f):

    """Sign function:
            maps positive values to 1,
                               0 to 0,
             and negative values to -1
    """

    # need to add a more computational float compare
    if f > 0:
        return 1
    elif f == 0:
        return 0
    else:
        return -1


def magpol_barycenters(mag_map):

    """
       Calculates the pixel coordinates of the "Magnetic Polarity barycenter"
       for a given magnetogram"
       inputs:
            mag_map: HMI sunpy.map.Map or sunpy.map.Map.submap

       output:
            dictionary with positive baricenter (x,y) pixel coordinates
            and negative baricenter (x,y) pixel coordinatss
    """

    # take the magnetic field data:
    data = mag_map.data
    # header = mag_map.meta

    rows, columns = data.shape
    x_neg_sum, x_pos_sum, y_neg_sum, y_pos_sum, tot_sum_pos, tot_sum_neg = (0., 0., 0., 0., 0., 0.)

    for i in range(0, rows-1):
        for j in range(0, columns-1):

            pix = data[i, j]
            apix = np.abs(pix)

            if sign(pix) == 1:
                x_pos_sum += j * apix
                y_pos_sum += i * apix
                tot_sum_pos += apix

            elif sign(pix) == -1:
                x_neg_sum += j * apix
                y_neg_sum += i * apix
                tot_sum_neg += apix

    x_pos = x_pos_sum / tot_sum_pos
    y_pos = y_pos_sum / tot_sum_pos
    x_neg = x_neg_sum / tot_sum_neg
    y_neg = y_neg_sum / tot_sum_neg

    return {"Pos": (x_pos, y_pos), "Neg": (x_neg, y_neg)}


def integral_to_minimize(s, theta,magmap, rp, rn):


    x_min = 0
    y_min = 0
    y_max, x_max = magmap.shape

    r_sy = rn[0] + s * (rp[0] - rn[0])
    r_sx = rn[1] + s * (rp[1] - rn[1])

    integral = dblquad(lambda x, y: np.square(np.abs(sign((x - r_sx) * np.sin(theta) + (y - r_sy) * np.cos(theta)) * np.abs(magmap[int(x), int(y)]) - magmap[int(x), int(y)])),
                       x_min,
                       x_max - 1,
                       lambda y: y_min,
                       lambda y: y_max - 1)

    return integral[0]


def pil_computation(mag_map):

    """

    :param mag_map:
    :return:
    """

    mag_barycenters = magpol_barycenters(mag_map)
    magnetogram = mag_map.data

    # take the coordinates of the magnetic polarities barycenters
    r_p = mag_barycenters["Pos"]
    r_n = mag_barycenters["Neg"]

    # Initial guess
    s_o = 0.5
    theta_o = np.pi / 4
    init_guess = np.array([s_o, theta_o])

    pil_comp = minimize(integral_to_minimize, x0=init_guess,
                        args=(magnetogram, r_p, r_n,),
                        bounds=[(0, 1), (0, np.pi)])

    return pil_comp
