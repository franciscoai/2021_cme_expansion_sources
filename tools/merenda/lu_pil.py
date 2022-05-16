# Polarity Inversion Line computation
# Author: From Poisson et al 2014 paper coded by Luciano A. Merenda

import ipdb
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


def pil_intersection(x, y, s, theta, r_p, r_n):

    r_sy = r_n[0] + s * (r_p[0] - r_n[0])
    r_sx = r_n[1] + s * (r_p[1] - r_n[1])

    return (x - r_sx) * np.sin(theta) + (y - r_sy) * np.cos(theta)


def pilf_to_integrate(x, y, s=0.0, theta=np.pi, mag=np.zeros((3, 4)), rp=(0, 0), rn=(0, 0)):

    print(x,y,s,theta,type(mag),rp,rn)
    y = int(y)
    x = int(x)
    return sign(np.square(pil_intersection(x, y, s, theta, rp, rn)) * np.abs(mag[y, x]) - mag[y, x])


def pil_to_minimize(s, theta, magmap=np.zeros((3, 4)), rp=(0, 0), rn=(0, 0)):

    bo = magmap
    x_min = 0
    y_min = 0
    y_max, x_max = bo.shape

    moreargs = [s, theta, magmap, rp, rn]

    integral = dblquad(pilf_to_integrate, y_min, y_max - 1, lambda x: x_min, lambda x: x_max - 1, args=moreargs)

    return integral


def pil_computation(mag_map):

    mag_barycenters = magpol_barycenters(mag_map)
    magnetogram = mag_map.data

    # take the coordinates of the magnetic polarities barycenters
    r_p = mag_barycenters["Pos"]
    r_n = mag_barycenters["Neg"]

    # Define the vector of the origin of the pil intersectrion
    s_o = 0.5
    theta_o = np.pi / 4.0
    ipdb.set_trace()
    pil_comp = minimize(pil_to_minimize, x0=np.array([s_o, theta_o]), args=[magnetogram, r_p, r_n])

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

    print(x_pos, y_pos, x_neg, y_neg, x_pos_sum, x_neg_sum, y_pos_sum, y_neg_sum, tot_sum_pos, tot_sum_neg)

    return {"Pos": (x_pos, y_pos), "Neg": (x_neg, y_neg)}
