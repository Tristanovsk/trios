import base64
import pandas as pd
import numpy as np
import glob
import io
import os
from textwrap import dedent as d
import re
import matplotlib as mpl
import plotly
import plotly.graph_objs as go

import utils.utils as u
from trios.process import *

# import aeronet
# from config import *


# ------------------------------------------------
# above-water data files
dirfig = os.path.abspath('/DATA/OBS2CO/data/trios/fig')
dirout = os.path.abspath('/DATA/OBS2CO/data/trios/above_water')
awrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/aw*idpr*.csv")
iwrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/uw*idpr*.csv")
swrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/Lu0*idpr*.csv")

coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ.csv")[0]
coords = pd.read_csv(coordf, sep=';')
coords
# get idpr numbers
idprs = np.unique([re.findall(r'idpr(\d+)', x)[0] for x in swrfiles])


# loop over idpr
for idpr in idprs:
    #    idpr=idprs[2]
    print(idpr)
    try:
        c = coords[coords.ID_prel == int(idpr)]  # .values[0]
        lat = c['Lat'].values[0]
        lon = c['Lon'].values[0]
        alt = 0 #c['Altitude']
        name = c['ID_lac'].values[0]
    except:
        continue
    dff = pd.DataFrame()

    # -----------------------------------------------
    #   AWR processing
    # -----------------------------------------------

    awr = u.awr_data(idpr, awrfiles)
    if awr.file:

        try:
            df, wl = awr.reader(lat, lon, alt, name)

            Rrs, rho = awr_process(df, wl).call_process(aot=0.05)

            print(rho.mean)
            # outputs saving
            dff = pd.concat([df, Rrs], axis=1)

        except:
            pass


    # -----------------------------------------------
    #   SWR processing
    # -----------------------------------------------

    swr = u.swr_data(idpr, swrfiles)
    if swr.file:
        df, wl = swr.reader(lat, lon, alt)
        Rrs = swr_process(df, wl).call_process()
        #Rrs.to_csv(os.path.join(dirout, 'trios_swr_' + name + '_idpr' + idpr + '.csv'))
        dff = pd.concat([dff, Rrs], axis=1)

    # -----------------------------------------------
    #   IWR processing
    # -----------------------------------------------
    # iwr = u.iwr_data(idpr, iwrfiles)
    # if iwr.file:
    #     df, wl = iwr.reader(c[1], c[2], c[3])
    #     Rrs = iwr_process(df, wl).trios()
    #     dff = pd.concat([dff, Rrs], axis=1)

    # writing output file
    dff.to_csv(os.path.join(dirout, 'trios_awr_' + name + '_idpr' + idpr + '.csv'))

