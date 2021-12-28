import pandas as pd

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

