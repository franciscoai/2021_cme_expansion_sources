#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
##############################################################################################################
This routine read the mass .dat files for the high latitude events. 
The .dat files are located in the server "/gehme/data/pacecraft/mass/detector/date/sector/" and have the format
file = "yyyymmdd_detector_mass_summary.dat"

The routine read the fils and make plots of Mass [g] versus distance and time

Resulting Figures are saved as ".png" in /gehme/data/pacecraft/mass/c2/yyyymmdd/sector/

Created on Fri Mar 22 12:37:42 2024
@author: Fernando Lopez - GEHMe
##############################################################################################################
"""

def dat2csv(file,path,date):
    with open(file) as dat_file, open(path+date+'_file_csv.csv','w') as csv_file:
        csv_writer = csv.writer(csv_file,delimiter=',')
        file=[]
        columns =[]
        for line in dat_file:
            row = [field.strip() for field in line.split(' ')]
            row = [l for l in row if len(l)>0]
            if len(columns) == 0:
                columns = row.copy()
                columns.insert(1,'time')
            else:
                file.append(row)
            csv_writer.writerow(row)
    df = pd.DataFrame(file,columns = columns)
    return(df)

import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
import csv

"Complete the following: "

# date to analyze "yyyymmdd": 
date = ["20101212","20101214","20110317","20110605","20130123","20130129","20130209","20130424","20130502", 
        "20130517","20130527","20130607"]

detectores = ['c2','cor1a','cor1b','cor2a','cor2b']
for date in date: 
    time =[]
    times = []
    dates = []
    mass = []
    dtype = []
    height = []
    aw = []
    for detector in detectores:
        # The path of the file.dat 
        if detector == 'c2':
            mission = 'soho/lasco'
            path = f"/gehme/gehme/data/soho/lasco/mass/{detector}/{date}/sector/"
        else:
            nave = detector[4]
            detector2 = detector[0:4]
            path = f"/gehme/gehme/data/stereo/secchi/mass/{nave}/{detector2}/{date}/sector/"
        file = path + date +"_"+ detector + "_mass_summary.dat"
        # Transform .dat file to csv file to read with pandas 
        df1 = dat2csv(file,path,date)
        df = df1[df1["ROI-Typ"] == "SECTOR"]
                
        times.append(df["time"])
        mass.append(df["Mass"])
        dates.append(df["Date-Obs"])
        height.append(df["Width"])
        aw.append(df["Cen-PA"])

    datec2 = dates[0]
    datecor1a = dates[1]
    datecor1b = dates[2]
    datecor2a = dates[3]
    datecor2b = dates[4]
    
    tc2 = times[0]
    tcor1a = times[1]
    tcor1b = times[2]
    tcor2a = times[3]
    tcor2b = times[4]
    
    mc2 = mass[0]
    mcor1a = mass[1]
    mcor1b = mass[2]
    mcor2a = mass[3]
    mcor2b = mass[4]
    
    hc2 = height[0]
    hcor1a = height[1]
    hcor1b = height[2]
    hcor2a = height[3]
    hcor2b = height[4]
    
    awc2 = aw[0]
    awcor1a = aw[1]
    awcor1b = aw[2]
    awcor2a = aw[3]
    awcor2b = aw[4]
    
    
    timec2 = pd.to_datetime(datec2+' '+tc2, format='%Y/%m/%d %H:%M:%S.%f',utc='True')   
    timecor1a = pd.to_datetime(datecor1a+' '+tcor1a, format='%Y/%m/%d %H:%M:%S.%f',utc='True')   
    timecor1b = pd.to_datetime(datecor1b+' '+tcor1b, format='%Y/%m/%d %H:%M:%S.%f',utc='True')   
    timecor2a = pd.to_datetime(datecor2a+' '+tcor2a, format='%Y/%m/%d %H:%M:%S.%f',utc='True')   
    timecor2b = pd.to_datetime(datecor2b+' '+tcor2b, format='%Y/%m/%d %H:%M:%S.%f',utc='True')   
    
    "###################################################################################################"
    # make the plot Mass vs time
    fig, ax = plt.subplots()
    ax.plot(timec2,np.log10((mc2.astype(float))),'-ro',label ='C2')
    ax.plot(timecor1a,np.log10((mcor1a.astype(float))),'-bo',label = 'COR1-A')
    ax.plot(timecor1b,np.log10((mcor1b.astype(float))),'-go',label = 'COR1-B')
    ax.plot(timecor2a,np.log10((mcor2a.astype(float))),'-bs',label = 'COR2-A')
    ax.plot(timecor2b,np.log10((mcor2b.astype(float))),'-gs',label = 'COR2-B')
    ax.legend(loc ="lower right")
    
    ax.set(xlabel=("Time"), ylabel=("log$_{10}$ (Mass) [g]"),title = ("Event - " + date))
    
    """
    ------------------- Enter the plot section ------------------------
    """
    plotdir = "/gehme/gehme/data/soho/lasco/mass/c2/"+date+"/sector/"
    
    # Plot mass vs time   
    figname = plotdir+date+"_mass_vs_time.png"
    plt.savefig(figname,dpi=300,format='png')
    #plt.show()
    "###################################################################################################"
    # make the plot Mass vs distance
    fig, ax1 = plt.subplots()
    ax1.plot(hc2.astype(float),np.log10((mc2.astype(float))),'-ro',label ='C2')
    ax1.plot(hcor1a.astype(float),np.log10((mcor1a.astype(float))),'-bo',label = 'COR1-A')
    ax1.plot(hcor1b.astype(float),np.log10((mcor1b.astype(float))),'-go',label = 'COR1-B')
    ax1.plot(hcor2a.astype(float),np.log10((mcor2a.astype(float))),'-bs',label = 'COR2-A')
    ax1.plot(hcor2b.astype(float),np.log10((mcor2b.astype(float))),'-gs',label = 'COR2-B')
    ax1.legend(loc ="lower right")
     
    ax1.set(xlabel=("Projected apex height (R$_\odot$)"), ylabel=("log$_{10}$ (Mass) [g]"),title = ("Event - " + date))
     
    #Plot mass vs height   
    figname = plotdir+date+"_mass_vs_height.png"
    plt.savefig(figname,dpi=300,format='png')
    #plt.show()
    #breakpoint()
    "###################################################################################################"
    # make the plot Mass vs AW
    # fig, ax2 = plt.subplots()
    # ax2.plot(awc2.astype(float),np.log10((mc2.astype(float))),'-ro',label ='C2')
    # ax2.plot(awcor1a.astype(float),np.log10((mcor1a.astype(float))),'-bo',label = 'COR1-A')
    # ax2.plot(awcor1b.astype(float),np.log10((mcor1b.astype(float))),'-go',label = 'COR1-B')
    # ax2.plot(awcor2a.astype(float),np.log10((mcor2a.astype(float))),'-bs',label = 'COR2-A')
    # ax2.plot(awcor2b.astype(float),np.log10((mcor2b.astype(float))),'-gs',label = 'COR2-B')
    # ax2.legend(loc ="lower right")
     
    # ax2.set(xlabel=("Angular Width [deg]"), ylabel=("log$_{10}$ (Mass) [g]"),title = ("Event - " + date))
     
    # #Plot mass vs height   
    # figname = plotdir+date+"_mass_vs_aw.png"
    # plt.savefig(figname,dpi=300,format='png')    