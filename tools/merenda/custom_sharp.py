from glob import glob
from sys import exit, argv
from time import sleep
import astropy.units as apu
import matplotlib.pyplot as plt
import pandas as pd
import sunpy.map
from sunpy.coordinates import HeliographicCarrington, propagate_with_solar_surface

from SHARPs import calculate_sharpkeys as cs
from .lu_tools import AreaSelector, cornerl_order

if __name__ == "__main__":

    # Select files to process
    # set the dir to work on, this can be changed as needed
    # will use argv to pass dir as an argument
    try:
        date = argv[1]
    except IndexError:
        date = input("Enter date to work on (YYYYMMDD format):\n")

    sharps_dir = "/gehme/data/sdo/hmi/sharps/"
    date_dir = date + "/"
    work_dir = sharps_dir + date_dir

    # Check in all files, ordered and split in segments.
    br_files = glob(work_dir + "*Br.fits")
    br_files.sort()
    bt_files = glob(work_dir + "*Bt.fits")
    bt_files.sort()
    bp_files = glob(work_dir + "*Bp.fits")
    bp_files.sort()
    br_err_files = glob(work_dir + "*Br_err.fits")
    br_err_files.sort()
    bt_err_files = glob(work_dir + "*Bt_err.fits")
    bt_err_files.sort()
    bp_err_files = glob(work_dir + "*Bp_err.fits")
    bp_err_files.sort()
    confdis_files = glob(work_dir + "*conf_disambig.fits")
    confdis_files.sort()
    bitmap_files = glob(work_dir + "*bitmap.fits")
    bitmap_files.sort()
    los_files = glob(work_dir + "*magnetogram.fits")
    los_files.sort()

    # Check if we have a directory with the exact
    # same amount of each type of file (segment)
    fc = len(br_files)  # files count
    if len(bt_files) == fc and len(bp_files) == fc and len(br_err_files) == fc \
            and len(bt_err_files) == fc and len(bp_err_files) == fc \
            and len(confdis_files) == fc and len(bitmap_files) == fc \
            and len(los_files) == fc:
        pass
    else:
        print("Wrong files count: files quantity must be the "
              "same for all type of files (segments)")
        exit(1)

    # select area to crop
    print("Choose the area to crop from the Br map")
    preview = sunpy.map.Map(br_files[0])
    preview.plot(title="Use mouse to select and press Q after selection is made")
    area_to_select = AreaSelector(plt.gcf(), plt.gca())

    # Give crop coordinates in bottomleft, topright order,
    # in order to use them with sunpy's map submap method.
    # Crop coordinates are selected for first image, all following crop
    # coordinates are corrected for solar rotation
    # To do: Option to select from which image to select crop.
    final_coords = cornerl_order(area_to_select.data_coords, "bltp", "np")

    bottomleft = preview.pixel_to_world(final_coords[0][0] * apu.pix,
                                        final_coords[0][1] * apu.pix)
    topright = preview.pixel_to_world(final_coords[1][0] * apu.pix,
                                      final_coords[1][1] * apu.pix)

    # Make the submap
    crop = preview.submap(bottom_left=bottomleft, top_right=topright)
    # Take a look at it. This is only for testing purposes.
    crop.peek()

    # Initialize arrays to store in SHARPs parameter computations
    # ids = [] Won´t add this for now
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
    bottomleft_arr = []
    topright_arr = []

    # Take all files and calculate sharp parameters using bobra's sharp module
    # and then save them in a csv file.
    for i in range(fc):
        # take the header for
        hdr = sunpy.map.Map(br_files[i]).meta

        # Adjust every frame crop coordinates for solar rotation
        if i == 0:
            new_final_coords = (bottomleft, topright)
        else:
            map_newframe = HeliographicCarrington(obstime=hdr['date-obs'])
            with propagate_with_solar_surface():
                new_final_coords = (bottomleft.transform_to(map_newframe), topright.transform_to(map_newframe))

        bz, by, bx, bz_err, by_err, bx_err, conf_disambig, bitmap, nx, ny, rsun_ref,\
            rsun_obs, cdelt1_arcsec, los, los_err = cs.get_data(br_files[i], bt_files[i],
                                                             bp_files[i],
                                                             br_err_files[i],
                                                             bt_err_files[i],
                                                             bp_err_files[i],
                                                             confdis_files[i],
                                                             bitmap_files[i],
                                                             los_files[i],
                                                             coords=new_final_coords)

        # compute the total unsigned flux and its error
        ausflux, ausflux_err, count_mask = cs.compute_abs_flux(bz, bz_err, conf_disambig,
                                                            bitmap, nx, ny, rsun_ref,
                                                            rsun_obs, cdelt1_arcsec)

        # compute the z current and associated errors and then
        # compute the magnetic current helicity
        current = cs.computeJz(bx, by, bx_err, by_err, conf_disambig, bitmap, nx, ny)
        jz, jz_err, derx, dery = current[0], current[1], current[2], current[3]

        ameanjzh, ameanjzh_err, atotusjh, atotusjh_err, aabsnjzh, aabsnjzh_err = \
            cs.computeHelicity(jz, jz_err, bz, bz_err, conf_disambig, bitmap, nx, ny,
                            rsun_ref, rsun_obs, cdelt1_arcsec)

        # compute the potential field the magnetic free energy
        potential = cs.greenpot(bz, nx, ny)
        bpx, bpy = potential[0], potential[1]

        ameanpot, ameanpot_err, atotpot, atotpot_err = cs.computeFreeEnergy(bx_err, by_err,
                                                                         bx, by, bpx,
                                                                         bpy, nx, ny,
                                                                         rsun_ref,
                                                                         rsun_obs,
                                                                         cdelt1_arcsec,
                                                                         conf_disambig,
                                                                         bitmap)


        #Appending zone
        date.append(hdr['t_rec'])
        usflux.append(ausflux)
        usflux_err.append(ausflux_err)
        meanjzh.append(ameanjzh)
        meanjzh_err.append(ameanjzh_err)
        totusjh.append(atotusjh)
        totusjh_err.append(atotusjh_err)
        absnjzh.append(aabsnjzh)
        absnjzh_err.append(aabsnjzh_err)
        meanpot.append(ameanpot)
        meanpot_err.append(ameanpot_err)
        totpot.append(atotpot)
        totpot_err.append(atotpot_err)
        bottomleft_arr.append([new_final_coords[0].lon.deg, new_final_coords[0].lat.deg])
        topright_arr.append([new_final_coords[1].lon.deg, new_final_coords[1].lat.deg])

    # I need to attach somehow the crop coordinates.
    # Create dataframe and store the data into a csv
    custom_sharp_df = pd.DataFrame({"T_REC": date,
                                    "USFLUX.Mx": usflux,
                                    "ERRVF.Mx": usflux_err,
                                    "MEANJZH.\\frac{G^{2}}{m}": meanjzh,
                                    "ERRMIH.\\frac{G^{2}}{m}": meanjzh_err,
                                    "TOTUSJH.\\frac{G^{2}}{m}": totusjh,
                                    "ERRTUI.\\frac{G^{2}}{m}": totusjh_err,
                                    "ABSNJZH.\\frac{G^{2}}{m}": absnjzh,
                                    "ERRTAI.\\frac{G^{2}}{m}": absnjzh_err,
                                    "MEANPOT.\\frac{erg}{cm^{3}}": meanpot,
                                    "ERRMPOT.\\frac{erg}{cm^{3}}": meanpot_err,
                                    "TOTPOT.\\frac{erg}{cm}": totpot,
                                    "ERRTPOT.\\frac{erg}{cm}": totpot_err,
                                    "BottomLeftCorner[lon,lat].deg": bottomleft_arr,
                                    "TopRightCorner[lon,lat].deg": topright_arr})

    custom_sharp_df.to_csv(work_dir + "croped_sharps.csv", sep=",", float_format="%04.3e", index=False)
