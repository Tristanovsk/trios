import os, sys
import pandas as pd
import numpy as np
import glob
import io

import matplotlib as mpl
import matplotlib.pyplot as plt

import plotly
import plotly.graph_objs as go
import cmocean

import utils.utils as u
from trios.process import *
from utils.sunposition import sunpos

dir = '/DATA/OBS2CO/data/sabine/data/raw'
odir = '/DATA/OBS2CO/data/sabine/data/awr'
file = os.path.join(dir, 'Ed_Ld_Lu_all_highroc_cruise_oslofjord.dat')

dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')

dff = pd.read_csv(file, sep="\t", date_parser=dateparse, index_col=2)

awr = awr_process()
awr.get_rho_values([30], [40], [135])

groupname = dff.columns.values[1] #'RT$station_area[n]'


def get_param(df, param='Ld', name='Lsky'):
    ind = df.columns.str.contains(param)
    df = df.loc[:, ind]
    wl = df.columns.str.replace(param + '_', '').astype('float')
    df.columns = pd.MultiIndex.from_product([[name], wl], names=['param', 'wl'])
    return df


for group, df_ in dff.groupby(groupname):
    print(group)

    ofile = os.path.join(odir, 'awr_' + group + '.csv')
    Rrs_df = pd.DataFrame()
    Rrs_ = pd.DataFrame()
    for time, df in df_.groupby(df_.index):
        azi = [int(df.iloc[:, 0].unique()[0].split('_')[-1])]
        lat = df.LAT.values.mean()
        lon = df.LON.values.mean()
        ws = [df.WIND_SPEED.values.mean()/2]
        # TODO check altitude
        alt = 0

        Lsky = get_param(df, 'Ld', 'Lsky')
        Lt = get_param(df, 'Lu', 'Lt')
        Ed = get_param(df, 'Ed', 'Ed')

        wl = Lsky.columns.levels[1].values

        # compute solar angle (mean between first and last acquisition time

        sza = sunpos(time, lat, lon, alt)[1]

        Rrs, rho = awr.process(wl, Lt, Lsky.values, Ed.values, sza, ws=ws, azi=azi)

        df.columns = pd.MultiIndex.from_product([df.columns.values, ['']], names=['param', 'wl'])
        Rrs.columns = pd.MultiIndex.from_product([['Rrs(awr)'], wl], names=['param', 'wl'])
        # ['Rrs(awr)'] * Rrs.shape[1]
        Rrs = pd.concat([df.iloc[:, 0:10], Ed,Lsky,Lt,Rrs], axis=1)
        Rrs['sza']=sza
        Rrs['azi']=azi[0]
        Rrs['rho']=rho.mean()
        Rrs_ = pd.concat([Rrs_, Rrs], axis=0)

    Rrs_df = pd.concat([Rrs_df, Rrs_], axis=0)
    Rrs_df.to_csv(ofile)

if False:
    dff = Rrs_df
    # dff.columns = pd.MultiIndex.from_product([['Rrs(awr)'], dff.columns], names=['param', 'wl'])
    # wl = dff.columns.droplevel().values
    trace = []
    i = 0

    for station, df in dff.groupby('RT$station_area[n]'):
        i = i + 1
        for time, df_ in df.iterrows():
            print(time)
            trace.append(go.Scattergl(
                x=wl,  # spectrum,
                y=df_.loc['Rrs(awr)'].values,
                text='',
                name=station,

                legendgroup=station,
                mode='lines',
                line={

                    'color': 'rgba' + str(cmocean.cm.thermal(i * 25)[:-1])}

                # 'line': go.scatter.Line(color='rgba({}, {}, {}, {})'.format(*s_m.to_rgba(parameters[i]).flatten()), width=2),

                , ))

    plotly.offline.plot(trace)
