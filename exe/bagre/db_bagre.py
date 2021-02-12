import os, sys
import pandas as pd
import numpy as np
import xarray as xr
import glob
import io
import matplotlib as mpl
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objs as go
import cmocean
import scipy.optimize as so
from scipy.interpolate import interp1d

from trios.process import *
from trios.utils.sunposition import sunpos

opj = os.path.join
dir = '/DATA/projet/bagre'
dirdata = opj(dir, 'data')
aerosols = ['fine', 'coarse']
aerosol = aerosols[1]
method = 'osoaa_' + aerosol
odir = opj(dirdata, 'L2', method)

files = glob.glob(opj(dirdata, 'raw/trios/csv_QC', 'trios*.csv'))

if not os.path.exists(odir):
    os.makedirs(odir)


# (for raw data: 83b0 = Lu, 853e = Ed, 83ae = Ld40va, 855dLd50va).
# sensors are pointing to a 240º azimuth (cw north).
# TODO, redo compuation for wl <400nm
#  project spectral data on common wavelength set (i.e., wl_common defined in trios.config)
# here we change wl_common to stay within interpolation range (i.e., 410-1025 nm)
# wl_common = wl_common[wl_common>=410]

def format_Rrs(Rrs, df, name='Rrs'):
    Rrs = Rrs.to_pandas()
    Rrs.index = df.index
    Rrs.columns = pd.MultiIndex.from_product([[name], df.Lt.columns, ])
    return Rrs


awr = awr_process(aerosol=aerosol)

vza = 40
ws = 2
aot550 = 0.1
rho = awr.rho.rho.to_xarray()
rho_ = rho.interp(vza=vza, wind=ws, aot=aot550, wl=wl_common)
rho_M99 = awr.rhoM1999.to_xarray().interp(vza=vza, wind=ws)
rho_M99 = rho_M99.interp(sza=np.linspace(0, 80, 81), method='cubic').interp(azi=np.linspace(0, 180, 181),
                                                                            method='cubic')

for file in files:
    basename = os.path.basename(file)
    df = pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)
    print(df.shape)
    lat = df.lat.values[0][0]
    lon = df.lon.values[0][0]
    alt = 0
    df['sza', ''] = np.nan
    df['azi', ''] = np.nan
    geom = sunpos(df.index.to_pydatetime(), lat, lon, alt)[:, 0:2]
    df.at[:, 'sza'] = geom[:, 1]
    relazi = 135  # ( geom[:,0] - azi_sensor)% 360
    # to get values between 0 and 180°
    # since rho values are symmetrical with the principal plane
    # relazi[relazi>180] = 360-relazi[relazi>180]
    df.at[:, 'azi'] = relazi
    print(basename)
    for name, raw in df.resample('1d'):
        print(name)

        # try:
        suff = name.__str__().replace(' ', '_')
        N = len(raw.index)
        if not N:
            print('No data for %s'.format(file))
            continue
        # ------------------
        # filtering
        # ------------------
        # daylight data
        ind = (raw.sza < 80)  & (raw.Lt.mean(axis=1) > 1) & (raw.Lsky.mean(axis=1) < 300 )
        if not any(ind):
            print('No QA checked data for ',file)
            continue
        # TODO add time rolling filters
        # ind = awr.filtering(raw.Lt, raw.Lsky, raw.Ed)
        clean = raw[ind].__deepcopy__()
        Lt, Lsky, Ed, sza, azi = clean.Lt.values, clean.Lsky.values, clean.Ed.values, clean.sza.values, clean.azi.values

        sza_ = xr.DataArray(sza, dims='geom')
        azi_ = xr.DataArray(azi, dims='geom')

        # -----------------------------
        # data processing
        # -----------------------------
        rho_v = rho_.interp(sza=sza_, azi=azi_).T
        rho_M99_ = rho_M99.interp(sza=sza_, azi=azi_).to_array().T.values

        clean['rho', ''] = rho_v.mean(axis=1)
        clean['rho_min', ''] = rho_v.min(axis=1)
        clean['rho_max', ''] = rho_v.max(axis=1)
        clean['rho_M99', ''] = rho_M99_

        Lsurf = (rho_v * Lsky)

        Lsurf_M99 = rho_M99_ * Lsky

        Rrs = format_Rrs((Lt - Lsurf) / clean.Ed, clean, 'Rrs_' + method)
        Rrs_M99 = (Lt - Lsurf_M99) / clean.Ed
        Rrs_M99.columns = pd.MultiIndex.from_product([['Rrs_M99'], Rrs_M99.columns, ])

        clean = pd.concat([Rrs, Rrs_M99, clean], axis=1)
        clean.to_csv(opj(odir, basename.split('.')[0]+'_L2_' + method + '_' + suff + '.csv'))

        # Rrs, rho = awr.process_wrapper(wl_common, raw, raw.sza, ws=ws, vza = [vza] * N, azi=raw.azi)
        #
        # awr.get_rho_values(raw.sza-90,[vza]*N, raw.azi.values, ws=[2])
        #
        # rho = awr.get_rho_mobley(awr.rhoM1999,raw.sza-90,[vza]*N, azi=raw.azi, ws=[2])[0,:,0]
