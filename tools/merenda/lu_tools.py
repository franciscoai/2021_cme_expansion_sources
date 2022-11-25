import pandas as pd
import gspread as gs  # use google apis to read cmes info from google spreadsheet
import os
import drms
import sunpy.map
import numpy as np
from astropy.time import Time
from matplotlib.patches import Rectangle
from matplotlib.backend_bases import MouseButton
from matplotlib.backend_bases import MouseEvent
import ipdb

# Random tools for working with solar physics data
# Author: MerendaC.LucianoA.@GEHMe

# note to self
# make an object that will have each property as an attribute or method(TBD)
# To make csv properties we will use the same standard
# ID/DateOfEvent, time, prop1_unit1, prop2_unit2, prop3_unit3, ...,  propn_unitn
# then a loader will take the data from the csv a loa into the object.


class AreaSelector:

    """
        Class to select an area from an image in matplotlib.
    """
    def __init__(self, figure, axes):
        # Initialize the instances
        self.canvas_coords = []
        self.data_coords =[]
        self.bottom_left = 0.0
        self.width = 0.0
        self.height = 0.0
        #We connect the figure to event handlers
        self.axs = axes
        self.fig = figure
        self.press_event = figure.canvas.mpl_connect('button_press_event', self.press)
        self.rel_event = figure.canvas.mpl_connect('button_release_event', self.release)
        self.move_event = figure.canvas.mpl_connect('motion_notify_event', self.move)

    def press(self, event):
        # Register first point of the selected area
        if event.button is MouseButton.LEFT:
            if len(self.canvas_coords) < 1 and len(self.data_coords) < 1:
                self.canvas_coords.append((event.x, event.y))
                self.data_coords.append((event.xdata, event.ydata))
                print(event.x, event.y, event.xdata, event.ydata)
            else:
                # Clear all items in order to start a new rectangle

                self.canvas_coords = []
                self.data_coords = []
                self.canvas_coords.append((event.x, event.y))
                self.data_coords.append((event.xdata, event.ydata))
                print(event.x, event.y, event.xdata, event.ydata)

    def move(self,event):

        # Draw the rectangle as the mouse moves over the image
        if event.button is MouseButton.LEFT:
            if len(self.axs.patches) < 1:
                bottom_left, width, height = cornerl_order(
                   [self.data_coords[0], (event.xdata, event.ydata)], "blwh", "mpl")
                area_to_select = Rectangle(bottom_left, width=width,
                                           height=height, linewidth=1.0, alpha=0.2, facecolor=None,edgecolor='black')
                self.axs.add_patch(area_to_select)
                self.fig.canvas.draw()

            else:
                self.axs.patches.pop()
                self.fig.canvas.draw()
                bottom_left, width, height = cornerl_order(
                    [self.data_coords[0], (event.xdata, event.ydata)], "blwh", "mpl")
                area_to_select = Rectangle(bottom_left, width=width,
                                           height=height, linewidth=1.0, alpha=0.2,
                                           facecolor=None, edgecolor='black')
                self.axs.add_patch(area_to_select)
                self.fig.canvas.draw()


    def release(self, event):
        # Register second point and draw final rectangle.
        if event.button is MouseButton.LEFT:
            self.canvas_coords.append((event.x, event.y))
            self.data_coords.append((event.xdata, event.ydata))

        self.bottom_left, self.width, self.height = cornerl_order(
            self.data_coords, "blwh", "mpl")
        print(self.data_coords,self.width,self.height)
        area_selected = Rectangle(self.bottom_left, width=self.width,
                                      height=self.height, linewidth=1.0, alpha=0.2, facecolor=None,edgecolor='black')
        self.axs.add_patch(area_selected)
        self.fig.canvas.draw()
        print(event.x, event.y,event.xdata,event.ydata)
        print(self.canvas_coords,self.data_coords)


def cornerl_order(coords, order, library):

    """ Takes a list of two tuples that contains the pixel coords of a rectangle and
        returns the same list in the following order[top_left_tuple, bottom_right_tuple]
        according to numpy indexing order or matplotlib too
        input:
                coords = [(x1,y1),(x2,y2)]
                order can be: "tpbr":  TopLeft - BottomRight order to use with opencv's rectangle
                              "bltr": BottomLeft - TopRight order use it with sunpy.map.submap
                              "blwh": BottomLeft width height order to use with matplotlib patches
                library: "np" for numpy and opencv
                         "mpl" for matplotlib
        output: coords in desired format/order
    """

    x1, y1 = coords[0]
    x2, y2 = coords[1]

    if (not x1) or (not y1) or (not x2) or (not y2):
        raise ValueError("Outside of Data")

    if library == 'np':
        if x2 > x1 and y2 < y1:
            top_left = (x1, y2)
            bottom_right = (x2, y1)
            top_right = (x2, y2)
            bottom_left = (x1, y1)
        elif x2 > x1 and y2 > y1:
            top_left = (x1, y1)
            bottom_right = (x1, y2)
            top_right = (x2, y1)
            bottom_left = (x1, y2)
        elif x1 > x2 and y1 > y2:
            top_left = (x2, y2)
            bottom_right = (x1, y1)
            top_right = (x1, y2)
            bottom_left = (x2, y1)
        elif x1 > x2 and y1 < y2:
            top_left = (x2, y1)
            bottom_right = (x1, y2)
            top_right = (x1, y1)
            bottom_left = (x2, y2)
        else:
            print("This isn't a rectangle!")
            if x1 == x2 and (y1 > y2 or y2 > y1):
                top_left = (x2, y1)
                bottom_right = (x2, y2)
                top_right = (x2, y1)
                bottom_left = (x2, y2)
            elif y1 == y2 and (x1 > x2 or x2 > x1):
                top_left = (x1, y2)
                bottom_right = (x2, y2)
                top_right = (x1, y2)
                bottom_left = (x2, y2)


    elif library == 'mpl':
        if x2 > x1 and y2 > y1:
            top_left = (x1, y2)
            bottom_right = (x2, y1)
            top_right = (x2, y2)
            bottom_left = (x1, y1)
        elif x2 > x1 and y2 < y1:
            top_left = (x1, y1)
            bottom_right = (x2, y2)
            top_right = (x2, y1)
            bottom_left = (x1, y2)
        elif x1 > x2 and y1 > y2:
            top_left = (x2, y1)
            bottom_right = (x1, y2)
            top_right = (x1, y1)
            bottom_left = (x2, y2)
        elif x1 > x2 and y1 < y2:
            top_left = (x2, y2)
            bottom_right = (x1, y1)
            top_right = (x1, y2)
            bottom_left = (x2, y1)
        else:
            print("This isn't a rectangle!")
            if x1 == x2 and (y1 > y2 or y2 > y1):
                top_left = (x2, y1)
                bottom_right = (x2, y2)
                top_right = (x2, y1)
                bottom_left = (x2, y2)
            elif y1 == y2 and (x1 > x2 or x2 > x1):
                top_left = (x1, y2)
                bottom_right = (x2, y2)
                top_right = (x1, y2)
                bottom_left = (x2, y2)
    else:
        raise ValueError("Wrong library input")

    if order == "tpbr":
        # This is for opencv's function cv2.rectangle
        return [top_left, bottom_right]
    elif order == "bltp":
        # This is for sunpy.map.submap's method
        return [bottom_left, top_right]
    elif order == "blwh":
        # this is for matplolib's Rectangle patch
        width = bottom_right[0] - bottom_left[0]
        height = top_left[1] - bottom_left[1]
        return (bottom_left,width,height)
    else:
        raise ValueError("Wrong order input")


def uglydateconv(uglydate):

    """
        Will convert just dd/MM/YYYY to YYYY-MM-dd  for now
    """
    dd = uglydate[0:2]
    mm = uglydate[3:5]
    yyyy = uglydate[6:]

    nicedate = yyyy + "-" + mm + "-" + dd

    return nicedate


def date2id(datetime, ids_dir):

    # ids_dir: path to Dates and Times excel downloaded file
    # put datetime in the time format needed
    datetime = datetime.replace(".", "-")

    # get table from google sheets or download it to some dir to get it
    # get all the sheets, select the columns we only need
    # and clean them from nan values added automatically for blank cells

    datesandtimes = pd.read_excel(ids_dir, sheet_name=None)
    id_dates = pd.concat([datesandtimes[key][["ID", "Fecha"]] for key in datesandtimes.keys()])
    id_dates = id_dates[id_dates["ID"].isna() != True]

    ids = list(id_dates["ID"])
    dates = list(id_dates["Fecha"])

    # Get needed id for the csv file
    id_required = [ids[n] for n, fecha in enumerate(dates) if Time(fecha) == Time(datetime)]

    if len(id_required) > 1:
        print(f"Ids found for {datetime}:\n {id_required}")
        raise Exception("\n### 2 or more Ids were found")
    else:
        return id_required[0]


def noaa2harpnum(noaa):

    """
        Returns the harpnum associated with a given noaa Active region number
            input:
                noaa: NOAA active region number
            output:
                harpnum associated.
    """
    # read harpnums database
    harpnum_db = pd.read_csv("http://jsoc.stanford.edu/doc/data/hmi/harpnum_to_noaa/all_harps_with_noaa_ars.txt",
                             sep=" ")
    # check the database of noaa vs harpnum and returns the harpnum associated with a noaa AR number.
    noaas = list(harpnum_db["NOAA_ARS"])
    harps = list(harpnum_db["HARPNUM"])
    harpnum_count = 0
    harpnuml = []
    i = 0  # aux variable.

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


def drms_download(start_time, end_time="", time_span="", cadence="", instrument="hmi", data_series="", segments="",
                  noaa="", wavelength="", user_email="lucianomerenda3@gmail.com"):
    """ vso type  client for downloading data from sdo and soho/mdi.
     Parameters:
       "Time format for start_time and end_time is 2010.12.12_00:00:00_TAI"
       instrument can mdi, hmi, aia,
       type can be sharp, L1 depends on the instrument
       segments depends on the type of data
    """
    # Check time window

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

    # Select data series
    if data_series == "":
        # print data series available of input instrument
        available_series = drms_client.series(instrument.lower(), full=True)
        print(f"Available data series for {instrument}:\n")
        for series in available_series.index:
            print(f"{series}) {available_series.name[series]} : {available_series.note[series]}")

        ds_selected = available_series.name[int(input("\nSelect desired data series: "))]
    else:
        ds_selected = data_series

    # Select segments of the data needed
    if segments == "":
        segments_available = drms_client.info(ds_selected).segments
        for i, segment in enumerate(segments_available.index):
            print(f"{i}) {segment}: {segments_available.note[segment]}")

        segments_selected = "{" + input("Enter segments name separated by a comma: ") + "}"
    else:
        segments_selected = "{" + segments + "}"

    # Check if we need a harpnum for the data series selected
    # ir if we need the wavelength for aia data
    if "HARPNUM" in drms_client.pkeys(ds_selected):
        if noaa == "":
            noaa_input = input("Enter NOAA number: ")
            harpnum = "[" + noaa2harpnum(noaa_input) + "]"  # only for SHARPs
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
        download_dir = "/gehme/data/sdo/hmi/sharps/" + start_time[0:10].replace(".", "") + "/"
    elif "smarp" in ds_selected:
        download_dir = "/gehme/data/soho/mdi/smarps/" + start_time[0:10].replace(".", "") + "/"
    elif "aia.lev1" in ds_selected:
        download_dir = "/gehme/data/sdo/aia/L1/" + wavelength + "/" + start_time[0:10].replace(".", "") + "/"
    elif "hmi.M" in ds_selected:
        download_dir = "/gehme/data/sdo/hmi/" + start_time[0:10].replace(".", "") + "/"
    else:
        download_dir = input(
            "Please manually enter the download directory\nRemember to use this convention: "
            "---> /gehme/data/observatory/intrument/../data_type.../date/: \n")

    try:
        os.makedirs(download_dir)
    except FileExistsError:
        pass

    print(f"Data will be downloaded in {download_dir}")
    export_request.download(download_dir)


def google_to_csv(gtitle, fname):

    """ Import data from google spreadsheets and save in a csv file
            # NOTE: before using, read gspread documentation on how to set up the google
                cloud app in order to finally be able to use it with your google account
            # IMPORTANT: Can't be used from within the server as far as I understand,
            #           since authentication needs a browser GUI running in the server
    """

    # using the gspread library
    savepath = os.getcwd() + "/" + fname
    # Start google client by auth Remember to put json auth file in /.config/gspread
    google_client = gs.oauth()

    spreadsheet = google_client.open(gtitle)
    worksheet = spreadsheet.sheet1  # If we only have 1 sheet, need to improve this
    dataframe = pd.DataFrame(worksheet.get_all_values())

    dataframe.to_csv(savepath, sep=",")
