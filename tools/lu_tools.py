import pandas as pd
import gspread as gs #use google apis to read cmes info from google spreadsheet
import os
import drms

#Random tools for working with solar physics data
#Author: MerendaLucianoA.@GEHMe

#read harpnums database
harpnum_db = pd.read_csv("http://jsoc.stanford.edu/doc/data/hmi/harpnum_to_noaa/all_harps_with_noaa_ars.txt",sep=" ")

def noaa2harpnum(noaa):
    #check the database of noaa vs harpnum and returns the harpnum associated with a noaa AR number.
    noaas = list(harpnum_db["NOAA_ARS"])
    harps = list(harpnum_db["HARPNUM"])
    harpnum_count = 0
    harpnuml = []
    i = 0 # aux variable.

    while i < len(noaas) - 1:
        if str(noaa) in noaas[i].split(sep=","):
            harpnuml.append(harps[i])
            harpnum_count += 1
        i += 1

    if harpnum_count > 1:
        print("There's more than 1 harpnum available")
        harpnum = input("Please enter manually which one you need: ")
    else:
        harpnum = str(harpnuml[0])

    return harpnum

def drms_download(start_time,end_time,instrument="hmi",type="sharp",data_series="",segments="",user_email="lucianomerenda3@gmail.com"):

    # vso type  client for downloading data from sdo and soho/mdi.
    # Parameters:
    #   "Time format for start_time and end_time is 2010.12.12_00:00:00_TAI"
    #   instrument can mdi, hmi, aia,
    #   type can be sharp, L1 depends on the instrument
    #   segments depends on the type of data


    # set the drms client
    drms_client = drms.Client(email=user_email, verbose=True)

    #Select data series
    if data_series == "":
        #print data series available of input instrument
        available_series = drms_client.series(instrument.lower(),full=True)
        print(f"Available data series for {instrument}:\n")
        for series in available_series.index:
            print(f"{series}) {available_series.name[series]} : {available_series.note[series]}")

        ds_selected = available_series.name[int(input("\nSelect desired data series: "))]
    else:
        ds_selected = data_series

    #Select segments of the data needed
    if segments == "":
        segments_available = drms_client.info(ds_selected).segments
        for i,segment in enumerate(segments_available.index):
            print(f"{i}) {segment}: {segments_available.note[segment]}")

        segments_selected = "{" + input("Enter segments name separated by a comma: ") + "}"
    else:
        segments_selected = segments

    #Check if we need a harpnum for the data series selected
    if "HARPNUM" in drms_client.pkeys(ds_selected):
        noaa_input = input("Enter NOAA number: ")
        harpnum = "[" + noaa2harpnum(noaa_input) + "]" #only for sharps and smarps
    else:
        harpnum = ""

    #Enter start and end time
    time_window = f"[{start_time}-{end_time}]"

    # Lets create the export request shall we?
    ds_needed = ds_selected + harpnum + time_window + segments_selected
    print("\nDRMS request: ds_needed\n")
    export_request = drms_client.export(ds_needed, method='url', protocol='fits')

    # Wait for the server to prepare requested files
    export_request.wait()

    # Then download files
    # Set the download directory first
    if "sharp" in ds_selected:
        download_dir = "/gehme/data/sdo/hmi/sharps/" + start_time[0:10].replace(".","") + "/"
    elif "smarp" in ds_selected:
        download_dir = "/gehme/data/soho/mdi/smarps/" + start_time[0:10].replace(".","") + "/"
    elif "aia.lev1" in ds_selected:
        download_dir = "/gehme/data/sdo/aia/L1/wavelength/" + start_time[0:10].replace(".","") + "/"
    elif "hmi.M" in ds_selected:
        download_dir = "/gehme/data/sdo/hmi/" + start_time[0:10].replace(".", "") + "/"
    else:
        download_dir = input("Please manually enter the download directory\nRemember to use this convention! ---> /gehme/data/observatory/intrument/../data_type.../date/: \n")

    print(f"Data will be downloaded in {download_dir}")
    export_request.download(download_dir)



def google_to_csv(gtitle,fname):
    #Import data from google spreadsheets and save in a csv file
    #NOTE: before using, read gspread documentation on how to set up the google cloud app
    #      in order to finally be able to use it with your google account
    #IMPORTANT: Can't be used from within the server as far a I understand,
    #           since autentication needs a browser GUI running in the server

    #using the gspread library
    savepath = os.getcwd() + "/" + fname
    #Start google client by auth Remember to put json auth file in /.config/gspread
    google_client = gs.oauth()

    spreadsheet = google_client.open(gtitle)
    worksheet = spreadsheet.sheet1 #If we only have 1 sheet, need to improve this
    dataframe = pd.DataFrame(worksheet.get_all_values())

    dataframe.to_csv(savepath, sep=",")
