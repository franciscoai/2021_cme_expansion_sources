import pandas as pd
import gspread as gs #use google apis to read cmes info from google spreadsheet
import os
import drms
import sunpy.map

#Random tools for working with solar physics data
#Author: MerendaC.LucianoA.@GEHMe

#note to self
#make an object that will have each property as an attribute or method(TBD)
#To make csv properties we will use the same standard
#ID/DateOfEvent, time, prop1_unit1, prop2_unit2, prop3_unit3, ...,  propn_unitn
#then a loader will take the data from the csv a loa into the object.


def sharp_parameter_csvgen(savedir):

    # Files can be a string of the dir where the needed files are, can be used with regex
    #So we take only one type of data segment since we only need the metadata, which is the same for all segments.
    #Levantamos todos los directorios donde hayan sharps

    print("Scanning sharps Dirs")
    sharp_dir = "/gehme/data/sdo/hmi/sharps/"
    dirs = os.listdir(sharp_dir)

    ID = []
    date = []
    usflux = []
    usflux_err = []
    meanjzh = []
    meanjzh_err = []
    totusjh = []
    totusjh_err = []
    absnjzh = []
    absnjzh_err = []
    meanpot = []
    meanpot_err = []
    totpot = []
    totpot_err = []

    #Abrimos carpeta por carpeta y vamos agregando a la listas los datos de los headers de cada directorio
    print("Starting to read files and headers.")
    for dir in dirs:
        #tomamos solo un tipo de archivo(data segments)
        files = sharp_dir + dir + "/*Br.fits"

        #read all the needed segments
        sharps_data = sunpy.map.Map(files,sequence=True)

        #take the metadata from the headers
        headers = sharps_data.all_meta()
        datetime = headers[0].get("T_REC")[0:10]

        #Retrive event id.
        id = date2id(datetime) #FUNCION A DEFINIR ----->>>> que nos tire el id en funcion de la fecha del evento
        #sino preguntar a Fran como hacemos con el tema de los ids.
        for hdr in headers:

            ID.append(id)
            date.append(hdr.get("T_REC"))
            usflux.append(hdr.get("USFLUX"))
            usflux_err.append(hdr.get("ERRVF"))
            meanjzh.append(hdr.get("MEANJZH"))
            meanjzh_err.append(hdr.get("ERRMIH"))
            totusjh.append(hdr.get("TOTUSJH"))
            totusjh_err.append(hdr.get("ERRTUI"))
            absnjzh.append(hdr.get("ABSNJZH"))
            absnjzh_err.append(hdr.get("ERRTAI"))
            meanpot.append(hdr.get("MEANPOT"))
            meanpot_err.append(hdr.get("ERRMPOT"))
            totpot.append(hdr.get("TOTPOT"))
            totpot_err.append(hdr.get("ERRTPOT"))



    print("Saving all SHARPs data in " + savedir)
    parameters = pd.DataFrame({"ID":ID,"T_REC":date,"USFLUX":usflux,"ERRVF":usflux_err,"MEANJZH":meanjzh,"ERRMIH":meanjzh_err,"TOTUSJH":totusjh,"ERRTUI":totusjh_err,"ABSNJZH":absnjzh,"ERRTAI":absnjzh_err,"MEANPOT":meanpot,"ERRMPOT":meanpot_err,"TOTPOT":totpot,"ERRTPOT":totpot_err})

    #Save to csv file in savedir
    parameters.to_csv(savedir + "sharps.csv")

    return parameters


def noaa2harpnum(noaa):

    # read harpnums database
    harpnum_db = pd.read_csv("http://jsoc.stanford.edu/doc/data/hmi/harpnum_to_noaa/all_harps_with_noaa_ars.txt",
                             sep=" ")
    #check the database of noaa vs harpnum and returns the harpnum associated with a noaa AR number.
    noaas = list(harpnum_db["NOAA_ARS"])
    harps = list(harpnum_db["HARPNUM"])
    harpnum_count = 0
    harpnuml = []
    i = 0 # aux variable.

    while i < len(noaas) - 1:
        if any([str(noaa_sub) in noaas[i].split(sep=",") for noaa_sub in noaa.split(sep=",")]):
            harpnuml.append(harps[i])
            harpnum_count += 1
        i += 1

    if harpnum_count > 1:
        print("There's more than 1 harpnum available: ")
        print(harpnuml)
        harpnum = input("Please enter manually which one you need: ")
    else:
        harpnum = str(harpnuml[0])

    return harpnum

def drms_download(start_time,end_time="",time_span="",cadence="",instrument="hmi",data_series="",segments="",noaa="",wavelength="",user_email="lucianomerenda3@gmail.com"):

    # vso type  client for downloading data from sdo and soho/mdi.
    # Parameters:
    #   "Time format for start_time and end_time is 2010.12.12_00:00:00_TAI"
    #   instrument can mdi, hmi, aia,
    #   type can be sharp, L1 depends on the instrument
    #   segments depends on the type of data

    #Check time window
    
    if end_time == "":
        if time_span != "":
            time_window = f"[{start_time}/{time_span}]"
        elif time_span != "" and cadence != "":
            time_window = f"[{start_time}/{time_span}@{cadence}]"
        else:
            raise RuntimeError("Please specify end_time or time_span or cadence & time_span of time")
    else:
        # Enter start and end timee
        time_window = f"[{start_time}-{end_time}]"

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
        segments_selected = "{" + segments + "}"

    #Check if we need a harpnum for the data series selected
    #ir if we need the wavelength for aia data
    if "HARPNUM" in drms_client.pkeys(ds_selected):
        if noaa == "":
            noaa_input = input("Enter NOAA number: ")
            harpnum = "[" + noaa2harpnum(noaa_input) + "]" #only for sharps
        else:
            harpnum = "[" + noaa2harpnum(noaa) + "]"
    else:
        harpnum = ""

    if "aia" in ds_selected:
        if wavelength == "":
            wavelength_s = "[? WAVELNTH = " + input("Enter wavelength in Angstroms needed (e.g. 193: ") + " ?]"
        else:
            wavelength_s = "[? WAVELNTH = " + wavelength + " ?]"
    else:
        wavelength_s = ""

    # Lets create the export request shall we?
    ds_needed = ds_selected + harpnum + time_window + wavelength_s + segments_selected
    print(f"\nDRMS request: {ds_needed}\n")
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
        download_dir = "/gehme/data/sdo/aia/L1/" + wavelength + "/" + start_time[0:10].replace(".","") + "/"
    elif "hmi.M" in ds_selected:
        download_dir = "/gehme/data/sdo/hmi/" + start_time[0:10].replace(".", "") + "/"
    else:
        download_dir = input("Please manually enter the download directory\nRemember to use this convention! ---> /gehme/data/observatory/intrument/../data_type.../date/: \n")

    try:
        os.makedirs(download_dir)
    except FileExistsError:
        pass

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
