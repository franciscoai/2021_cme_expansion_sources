import numpy as np

import poisson_pil as pp
import sunpy.map

hebe_filament_times = ["2010/12/14 13:37:00","2011/06/05 03:55:00",
                       "2013/01/29 01:18:00","2013/02/09 05:40:00",
                       "2013/05/02 04:58:00","2013/06/07 22:11:00"]

hebe_filaments_data = ["/gehme/data/sdo/hmi/sharps/20101214/hmi.sharp_cea_720s.284.20101214_133600_TAI.magnetogram.fits",
                       "/gehme/data/sdo/hmi/sharps/20110605/hmi.sharp_cea_720s.650.20110605_040000_TAI.magnetogram.fits",
                       "/gehme/data/sdo/hmi/sharps/20130129/hmi.sharp_cea_720s.2414.20130129_012400_TAI.magnetogram.fits",
                       "/gehme/data/sdo/hmi/sharps/20130209/hmi.sharp_cea_720s.2460.20130209_053600_TAI.magnetogram.fits",
                       "/gehme/data/sdo/hmi/sharps/20130502/hmi.sharp_cea_720s.2693.20130502_050000_TAI.magnetogram.fits",
                       "/gehme/data/sdo/hmi/sharps/20130607/hmi.sharp_cea_720s.2790.20130607_221200_TAI.magnetogram.fits"]


for time, data in zip(hebe_filament_times, hebe_filaments_data):

    magmap = sunpy.map.Map(data)
    hdr = magmap.meta
    input_data = magmap.data
    pil_computation = pp.PILOOP(input_data)
    print(f"Filament activation time: {time},", f"Header T_REC: {hdr['T_REC'],}", f"theta0 = {pil_computation['theta0']}")

