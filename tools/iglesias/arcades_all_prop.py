# reads arcades_prop_csv and plots some quantitites
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import datetime as dt
import os

prop_to_plot = ['width fit vel [km/s]','length fit vel [km/s]','tilt fit mean value [deg]']
local_path = os.getcwd()
arcades_prop_csv = local_path + '/output_data/arcades/props/all_arcades_props.csv'
#reads csv
df = pd.read_csv(arcades_prop_csv)
# plot props vs event
for prop in prop_to_plot:
    plt.figure()
    plt.plot(df['event'], df[prop], 'o')
    plt.xlabel('event')
    plt.title(prop)
    plt.savefig(local_path + '/output_data/arcades/props/'+prop[0:-7]+'.png')

