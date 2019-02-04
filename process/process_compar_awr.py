import base64
import pandas as pd
import numpy as np
import glob
import io
import os
from textwrap import dedent as d
import re
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objs as go
from scipy.interpolate import interp1d

from utils.sunposition import sunpos
import utils.utils as u
import utils.auxdata as ua
from process.process import *

coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ.csv")[0]
coords = pd.read_csv(coordf, sep=';')
dirfig = os.path.abspath('/DATA/OBS2CO/data/trios/fig')

awrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/aw*idpr*.csv")
# awrfiles = glob.glob("/DATA/OBS2CO/data/trios/test_setup/raw/aw*idpr*.csv")
swrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/Lu0*idpr*.csv")

iopw = ua.iopw()
iopw.load_iopw()

def add_curve(ax, x, mean, std, c='red', label=''):
    ax.plot(x, mean, linestyle='solid', c=c, lw=2.5,
            alpha=0.8, label=label)
    ax.fill_between(x,
                    mean - std,
                    mean + std, alpha=0.35, color=c)


idpr = '167'

# get idpr numbers
idprs = np.unique([re.findall(r'idpr(\d+)', x)[0] for x in swrfiles])
#idprs = np.array(['170'])
# loop over idpr
for idpr in idprs:
    c = coords[coords.ID_prel == int(idpr)]  # .values[0]
    lat = c['Lat'].values[0]
    lon = c['Lon'].values[0]
    alt = 0  # c['Altitude']
    name = c['ID_lac'].values[0]

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))
    fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.65)

    # -----------------------------------------------
    #   SWR processing
    # -----------------------------------------------

    uswr = u.swr_data(idpr, swrfiles)
    if uswr.file:
        df, wl_swr = uswr.reader(lat, lon, alt)
        df['sza', ''] = np.nan
        for index, row in df.iterrows():
            # print index
            sza = sunpos(index, lat, lon, alt)[1]
            df.at[index, 'sza'] = sza
        swr = swr_process(df, wl_swr)
        Rrs_swr = swr.call_process()
        add_curve(ax, wl_swr, Rrs_swr.transpose().mean(axis=1), Rrs_swr.transpose().std(axis=1), label='swr', c='black')
        Rrs_swr = swr.call_process(shade_corr=True)
        add_curve(ax, wl_swr, Rrs_swr.transpose().mean(axis=1), Rrs_swr.transpose().std(axis=1), label='swr', c='red')

    # -----------------------------------------------
    #   AWR processing
    # -----------------------------------------------
    azi = 135
    vza = 40
    awr = u.awr_data(idpr, awrfiles)
    if awr.Edf:

        index_idx = [0]

        d = u.data(index_idx)
        Ed, wl_Ed = d.load_csv(awr.Edf)
        Lsky, wl_Lsky = d.load_csv(awr.Lskyf)
        Lt, wl_Lt = d.load_csv(awr.Ltf)

        # ''' interpolate Ed and Lsky data upon Lt wavelength'''
        wl = wl_Lt
        Lt.columns = pd.MultiIndex.from_tuples(zip(['Lt'] * len(wl), wl), names=['param', 'wl'])
        intEd = interp1d(wl_Ed, Ed.values, fill_value='extrapolate')(wl)
        newEd = pd.DataFrame(index=Ed.index,
                             columns=pd.MultiIndex.from_tuples(zip(['Ed'] * len(wl), wl), names=['param', 'wl']),
                             data=intEd)
        intLsky = interp1d(wl_Lsky, Lsky.values, fill_value='extrapolate')(wl)
        newLsky = pd.DataFrame(index=Lsky.index, columns=pd.MultiIndex.from_tuples(zip(['Lsky'] * len(wl), wl),
                                                                                   names=['param', 'wl']), data=intLsky)

        awr = awr_process()
        ws = [2]

        print(azi, vza)

        Lsky = newLsky  # .loc[(newLsky.index.get_level_values(1) ==  vza) & (newLsky.index.get_level_values(2) ==  azi)]
        Ed = newEd  # .loc[(newEd.index.get_level_values(1) ==  vza) & (newEd.index.get_level_values(2) ==  azi)]

        # Lsky_idx = Lsky.index
        # Ed_idx= Ed.index
        # Lt_idx = Lt.index
        # Lsky.reset_index(level=[1,2],inplace=True)
        # Ed.reset_index(level=[1,2],inplace=True)
        # Lt.reset_index(level=[1,2],inplace=True)

        # merge sensor data on time
        df = pd.merge_asof(Lt, Ed, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")
        df = pd.merge_asof(df, Lsky, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")

        # add solar angle data and idpr
        # compute solar angle (mean between fisrt and last aqcuisition time
        df['sza', ''] = np.nan
        for index, row in df.iterrows():
            # print index
            sza = sunpos(index, lat, lon, alt)[1]
            df.at[index, 'sza'] = sza

        rho_h = awr.get_rho_values([df.sza.mean()], [vza], [azi], wl=wl)
        rho15 = awr.get_rho_mobley(awr.rhoM2015, [df.sza.mean()], [vza], [azi], [ws])
        rho99 = awr.get_rho_mobley(awr.rhoM1999, [df.sza.mean()], [vza], [azi], [ws])

        Rrs_h = (df.loc[:, 'Lt'] - rho_h * df.loc[:, 'Lsky']) / df.loc[:, 'Ed']
        Rrs15 = (df.loc[:, 'Lt'] - rho15 * df.loc[:, 'Lsky']) / df.loc[:, 'Ed']

        Rrs99 = (df.loc[:, 'Lt'] - rho99 * df.loc[:, 'Lsky']) / df.loc[:, 'Ed']
        # plt.figure()

        add_curve(ax, wl, Rrs15.transpose().mean(axis=1), Rrs15.transpose().std(axis=1),
                  label='M2015 (' + str(rho15) + ')')
        add_curve(ax, wl, Rrs99.transpose().mean(axis=1), Rrs99.transpose().std(axis=1), c='orange',
                  label='M1999(' + str(rho99) + ')')
        add_curve(ax, wl, Rrs_h.transpose().mean(axis=1), Rrs_h.transpose().std(axis=1), c='grey',
                  label='h(' + str(rho_h.mean()) + ')')

    ax.set_title('azi=' + str(azi) + ', vza=' + str(vza) + ', sza=' + str(sza))

    ax.legend(loc='best', frameon=False)

    ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
    ax.set_xlabel(r'Wavelength (nm)')
    fig.savefig(os.path.join(dirfig, 'trios_awr_' + name + '_idpr' + idpr + '.png'))
    plt.close()

