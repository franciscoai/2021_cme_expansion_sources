import pandas as pd
import gspread as gs #use google apis to read cmes info from google spreadsheet
import os

#Random tools for wrking with solar physics data.

#read harpnums database
harpnum_db= pd.read_csv("http://jsoc.stanford.edu/doc/data/hmi/harpnum_to_noaa/all_harps_with_noaa_ars.txt",sep=" ")

def noaa2harpnum(noaa):
    #check the database of noaa vs harpnum and returns the harpnum associated with a noaa AR number.
    noaas = list(harpnum_db["NOAA_ARS"])
    harps = list(harpnum_db["HARPNUM"])
    harpnum = []
    i = 0

    while i < len(noaas) - 1:
        if str(noaa) in noaas[i].split(sep=","):
            harpnum.append(harps[i])
        i += 1
    return harpnum

"Dates and Times CMEs"

def google_to_csv(gtitle,fname):
    #Import data from google spreadsheets and save in a csv file
    #using the gspread library
    savepath = os.getcwd() + "/" + fname
    #Start google client by auth Remember to put json auth file in /.config/gspread
    google_client = gs.oauth()

    spreadsheet = google_client.open(gtitle)
    worksheet = spreadsheet.sheet1 #If we only have 1 sheet, need to improve this
    dataframe = pd.DataFrame(worksheet.get_all_values())

    dataframe.to_csv(savepath, sep=",")

