#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module provides plotting support for fits data cubes

@author: iglesias
"""
import cv2
import numpy as np
import copy
import scipy.ndimage as ndimage
import matplotlib
import matplotlib.pyplot as plt
import logging
import matplotlib as mpl

def save_img(img, path, rng=[-3, 3], minmax=[-1., -1.], resize=-1, intord=0,
                   title='No title defined', xlabel='Pixel', ylabel='Pixel',
                   cblabel='Signal [DN]', fontsz=25, doplot=0, drng=[0, 0, 0, 0], pdf=None, figsize=[0, 0]):
    """
    Francisco  A. Iglesias (franciscoaiglesias@hotmail.com), 2017.04.25

    Generic plotting routine: Saves as an image the plot of the input numpy array 
    keeping the original size, adding a color bar and pixel axes. Note that:
        -No interpolation is done by default
        -No resizeing is done by default
        -The output image is larger to fit the color bar and the axes'
        -All the pixels in the input are plotted, no more no less (if resize=-1)
        -NaN values are ignored to compute image statistics. This is useful to draw 
         lines or areas on top of the image

    INPUTS
    img(Numpy arr): Array to plot, 2D [x,y] or 3D [x,y,z]. The diff z are saved as separated 
                    images with _z at the end.
    path(string): Full path to the output file. The file name defines the output image
                  format. Options are: '.png' (tested), '.eps', '.pdf', etc
                  if pdf is set, the the image is not saved on disc but a pdf.savefig() is run                  

    OPTIONAL INPUTS
    rng([int,int]): Range of the color bar in fractions of the spatial
                    std deviaiton [min, max].
    minmax=[float,float]: If this is set, then rng is ignored and the min and max of the color bar 
                          are minmax[0] and minmax[1], respectively
    drng [xmin,xmax,ymin,ymax]: Data range to annotate the y and x axes
    resize(int): Set to resize the img array before plotting using nearest interpolation.
               Use -1 (default) to avoid resizeing.
    intord(int): Order of the spline used for interpolating when resizeing.
    title (string): Main title. Note that a new line showing the spatial mean
                    and stddev is alwways added below title.
                    If img is 3D title must be a string vector with z titles, one per image
    xlable, ylable and cblabel (string):
        To include as x and y axes labels, and color bar label
    fontsz(int): Font size. When resize!=-1, then efective font used is fontsz*resize
    doplot: Set to plot on screen instead of saving to file. If plot is X<0 then after ploting the image, 
            the plot is closed after X seconds.
    pdf:  matplotlib.backends.backend_pdf class to plot to a pdf document

    OUTPUTS
    The plot of img is saved in path
    """

    # CONSTANTS
    DPI = 96.
    MARGIN = [0.2, 0.17]  # Margins to allocate cbar and axes (not symetric)
    NTICKS = 5  # number of ticks
    SCIFMT = '{:8.4e}'  # format string to print sci notation numbers

    bend = mpl.get_backend()
    mpl.use('Agg')

    # MAIN
    if resize != -1:
        oimg = ndimage.zoom(img, resize, order=intord)
        fontsz *= resize
    else:
        oimg = img

    height, width = oimg.shape

    # figure size in inches
    if np.all(figsize) == 0:
        figsize = width / \
            float(DPI) / (1. - MARGIN[0]), height / \
            float(DPI) / (1. - MARGIN[1])

    # img mean and stdddev
    mimg = np.nanmean(oimg)
    stdimg = np.nanstd(oimg)

    # Create a figure of the right size with one axis that takes up the full figure
    fig = plt.figure(figsize=figsize)
    mpl.rcParams.update({'font.size': fontsz})
    ax = fig.add_axes([MARGIN[0] / 3., MARGIN[1] / 3.,
                       (1. - MARGIN[0]), (1. - MARGIN[1])])

    # Display the image.
    if minmax[0] == minmax[1]:
        im = ax.imshow(oimg, cmap='Greys_r', vmin=mimg + rng[0] * stdimg,
                       vmax=mimg + rng[1] * stdimg, interpolation='none')
    else:
        im = ax.imshow(oimg, cmap='Greys_r', vmin=minmax[0],
                       vmax=minmax[1], interpolation='none')
    # Add axis for the color bar and plot it
    box = ax.get_position()
    cax = plt.axes([box.x0 + box.width * 1.01, box.y0,
                    MARGIN[0] / 8., box.height])
    cbar = plt.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_label(cblabel)
    ax.set_title(title + '\n mean=' + SCIFMT.format(mimg) + '; std=' + SCIFMT.format(stdimg)
                 + '\n min=' + SCIFMT.format(np.nanmin(oimg)) + '; max=' + SCIFMT.format(np.nanmax(oimg)))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # add custom ticks
    if drng != [0, 0, 0, 0]:
        ax.set_xticks(np.arange(0, width+width/NTICKS, width/NTICKS))
        bla = np.arange(drng[0], drng[1]+(drng[1]-drng[0]) /
                        NTICKS, (drng[1]-drng[0])/NTICKS)
        ax.set_xticklabels(np.char.mod('%.2f', bla))
        ax.set_yticks(np.arange(0, height+height/NTICKS, height/NTICKS))
        bla = np.arange(drng[2], drng[3]+(drng[3]-drng[2]) /
                        NTICKS, (drng[3]-drng[2])/NTICKS)
        ax.set_yticklabels(np.char.mod('%.2f', bla))
    # Ensure we're displaying with square pixels and the right extent.
    # ax.set(xlim=[0, width], ylim=[height, 0], aspect=1)

    if doplot != 0:
        if doplot < 0:
            plt.ion()
            plt.pause(-doplot)
        else:
            plt.show()
    else:
        if pdf == None:
            plt.savefig(path, dpi=DPI)
        else:
            pdf.savefig()
        plt.close()
    mpl.rcdefaults()
    mpl.use(bend)
    return

def save_movie(path, imgs, step=0, fps=24):  
    """
    Saves a video to path with all images in imgs

    :param: imgs: Numpy array with the images [n,y,x]
    :param: path: Full path to output video file (ext will be added)
    :param: step: If set to any int > 0, the output video will include only one of every step frames
    :param: fps: Frames Per Second
    """ 
    codec = 'DIVX' # 'MP42'
    ext = '.avi'

    if step > 0:
        out_imgs = imgs[::step,:,:]
    else:
        out_imgs = imgs
    sz = np.shape(out_imgs)    
    print('Saving array with shape %s to a movie in %s', sz, path + ext) 
    width = sz[2]
    hieght = sz[1] 
    fourcc = cv2.VideoWriter_fourcc(*codec)      
    video = cv2.VideoWriter(path + ext, fourcc, float(fps), (width, hieght), False)
    for frm in range(sz[0]):
        frm_data = copy.copy(out_imgs[frm,:,:])
        frm_data = ((frm_data - frm_data.min()) * (1/(frm_data.max() - frm_data.min()) * 255)).astype('uint8')    
        frm_data = frm_data.astype('uint8')
        video.write(frm_data)
    video.release()
    return