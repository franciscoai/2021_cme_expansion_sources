import sunpy.map
import os
import lu_tools as lt
import pandas as pd


if __name__== "__main__":

    # Levantamos todos los directorios donde hayan SHARPs
    # y chequeamos que el directorio al repo git este bien.

    print("Scanning SHARPs Dirs")
    sharp_dir = "/gehme/data/sdo/hmi/sharps/"
    dirs = os.listdir(sharp_dir)
    git_dir_name = "2021_cme_expansion_sources"

    # Chequeamos path al repo git del proyecto.
    if git_dir_name != os.getcwd()[-len(git_dir_name):]:
        git_dir = input(f"Enter {git_dir_name} full path:")
    else:
        git_dir = os.getcwd()

    savedir = git_dir + "/input_data/"

    # Defino listas para almacenar los datos a escribir
    ids = []
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

    # Abrimos carpeta por carpeta y vamos agregando a la listas los datos de los headers de cada directorio
    print("Starting to read files and headers.")
    for directory in dirs:

        # tomamos solo un tipo de archivo(data segments) dado que solo
        # necesitamos la metadata que es la misma para cada segmento.
        files = sharp_dir + directory + "/*Br.fits"

        # read all the needed segments
        sharps_data = sunpy.map.Map(files, sequence=True)

        # take the metadata from the headers
        headers = sharps_data.all_meta()
        datetime = headers[0].get("T_REC")[0:10]

        # Retrive event id from the database.
        id_req = lt.date2id(datetime, ids_dir="/home/lugem/GEHMe/lucho_repo/Dates and Times CMEs.xlsx")

        for hdr in headers:
            ids.append(int(id_req))
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
    parameters = pd.DataFrame(
        {"ID": ids, "T_REC": date, "USFLUX": usflux, "ERRVF": usflux_err, "MEANJZH": meanjzh, "ERRMIH": meanjzh_err,
         "TOTUSJH": totusjh, "ERRTUI": totusjh_err, "ABSNJZH": absnjzh, "ERRTAI": absnjzh_err, "MEANPOT": meanpot,
         "ERRMPOT": meanpot_err, "TOTPOT": totpot, "ERRTPOT": totpot_err})


    # Rename the columns so the units are next to the parameter
    parameters = parameters.rename(columns={"USFLUX": "USFLUX.Mx", "ERRVF": "ERRVF.Mx",
                                            "MEANJZH": "MEANJZH.\\frac{G^{2}}{m}",
                                            "ERRMIH": "ERRMIH.\\frac{G^{2}}{m}", "TOTUSJH": "TOTUSJH.\\frac{G^{2}}{m}",
                                            "ERRTUI": "ERRTUI.\\frac{G^{2}}{m}", "ABSNJZH": "ABSNJZH.\\frac{G^{2}}{m}",
                                            "ERRTAI": "ERRTAI.\\frac{G^{2}}{m}", "MEANPOT": "MEANPOT.\\frac{erg}{cm^{3}}",
                                            "ERRMPOT": "ERRMPOT.\\frac{erg}{cm^{3}}", "TOTPOT": "TOTPOT.\\frac{erg}{cm}",
                                            "ERRTPOT": "ERRTPOT.\\frac{erg}{cm}"})

    # Data preview

    print("Data preview: \n")
    print(parameters.head())

    # Save to csv file in savedir
    parameters.to_csv(savedir + "SHARPs.csv", sep=",", float_format="%04.3e", index=False)
