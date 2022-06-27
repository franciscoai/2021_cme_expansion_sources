# Polarity Inversion Line computation
# Author: From Poisson et al 2014 paper coded by Luciano A. Merenda

from scipy.optimize import minimize
from scipy.integrate import dblquad
import numpy as np


def sign(f):

    # Sign function: maps positive values to 1, 0 to 0 and negative values to -1
    if f > 0:
        return 1
    elif f == 0:
        return 0
    else:
        return -1


def pil_to_minimize(s, t, magmap, rp, rn):

    x_min = 0
    y_min = 0
    y_max, x_max = magmap.shape

    r_sy = rn[0] + s * (rp[0] - rn[0])
    r_sx = rn[1] + s * (rp[1] - rn[1])
    theta = t

    integral = dblquad(lambda y, x: np.square(np.abs(sign((x - r_sx) * np.sin(theta) +
                                                          (y - r_sy) * np.cos(theta)) * np.abs(magmap[int(x), int(y)])
                                                     - magmap[int(x), int(y)])),
                       y_min,
                       y_max - 1,
                       lambda x: x_min,
                       lambda x: x_max - 1)

    return integral[0]


def pil_computation(mag_map):

    mag_barycenters = magpol_barycenters(mag_map)
    magnetogram = mag_map.data

    # take the coordinates of the magnetic polarities barycenters
    r_p = mag_barycenters["Pos"]
    r_n = mag_barycenters["Neg"]

    # Initial guess
    s_o = 1.5
    theta_o = np.pi / 4.0
    init_guess = np.array([s_o, theta_o])

    pil_comp = minimize(pil_to_minimize, x0=init_guess,
                        args=(magnetogram, r_p, r_n),
                        method='BFGS',)

    return pil_comp


def magpol_barycenters(mag_map):

    # Calculate the pixel coordinates of the "Magnetic Polarity barycenter" for a given magnetogram"
    # mag_map: HMI sunpy.map.Map or sunpy.map.Map.submap

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

    # Debug:
    # print(x_pos, y_pos, x_neg, y_neg, x_pos_sum, x_neg_sum, y_pos_sum, y_neg_sum, tot_sum_pos, tot_sum_neg)

    return {"Pos": (x_pos, y_pos), "Neg": (x_neg, y_neg)}
