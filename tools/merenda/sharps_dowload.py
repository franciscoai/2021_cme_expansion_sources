import pandas as pd
import lu_tools as lt

# Script para descargar los datos de SHARPs en el servidor
if __name__ == "__main__":

    # Tomamos los datos del google sheet "Dates and Times CMEs" y los refinames manualmente para obtener
    # los datos pertinentes a los SHARPs.
    sharp_data = pd.read_csv("sharp_download.csv")

    # Parametros comunes para todas las descargas
    sharps_series = "hmi.sharp_cea_720s"
    sharp_segments = "Br,Br_err,Bp,Bp_err,Bt,Bt_err,bitmap,conf_disambig,magnetogram"

    for row in sharp_data.index:

        # tomamos datos fila a fila para cada cme y generemos el download request.
        inicio = sharp_data.loc[row,"Start-fecha"] + "_" + sharp_data.loc[row,"Start-Time"] + "_TAI"
        fin = sharp_data.loc[row,"End-fecha"] + "_" + sharp_data.loc[row,"End-Time"] + "_TAI"
        noaa_number = sharp_data.loc[row,"NOAA AR"]

        lt.drms_download(start_time=inicio,end_time=fin,data_series=sharps_series,segments=sharp_segments,noaa=noaa_number)









