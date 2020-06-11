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
dir = '/DATA/projet/borges'
dirdata = opj(dir,'data')
method = 'osoaa_coarse'
odir = opj(dirdata,'L2',method)
if not os.path.exists(odir):
    os.makedirs(odir)

# (for raw data: 83b0 = Lu, 853e = Ed, 83ae = Ld40va, 855dLd50va).
# sensors are pointing to a 240º azimuth (cw north).
# TODO, redo compuation for wl <400nm
#  project spectral data on common wavelength set (i.e., wl_common defined in trios.config)
# here we change wl_common to stay within interpolation range (i.e., 410-1025 nm)
wl_common = wl_common[wl_common>=410]

lat, lon, alt = -16.2, -47.32, 850
azi_sensor = 240

def load_csv(file,label=''):
    ''' Load and reproject data on common wavelength set'''
    print(file)
    df = pd.read_csv(file,index_col=0,parse_dates=True,)
    wl = df.columns = df.columns.astype('float')
    df.index.name = 'date'

    raw = interp1d(wl, df.values, fill_value='extrapolate')(wl_common)
    df = pd.DataFrame(index=df.index,columns=pd.MultiIndex.from_tuples(zip([label] * len(wl_common), wl_common),
                             names=['param', 'wl']),data=raw)
    # sort to get data in increasing time order
    df.sort_index(inplace=True)
    return df

Ed = load_csv(opj(dirdata,'raw_data.xlsx_853e.csv'),label='Ed')
Lsky = load_csv(opj(dirdata,'raw_data.xlsx_83ae.csv'),label='Lsky')
#Lsky50,wl = load_csv(opj(dirdata,'raw_data.xlsx_855d.csv'),label='Lsky50')
Lt = load_csv(opj(dirdata,'raw_data.xlsx_83b0.csv'),label='Lt')

# merge sensor data on time
df = pd.merge_asof(Lt, Ed, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                   direction="nearest")
df = pd.merge_asof(df, Lsky, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                   direction="nearest")


# Convert to UTC time
df.index = df.index + pd.Timedelta(hours=3)

df['sza', ''] = np.nan
df['azi', ''] = np.nan
geom = sunpos(df.index.to_pydatetime(), lat, lon, alt)[:,0:2]
df.at[:, 'sza'] = geom[:,1]
relazi =  ( geom[:,0] - azi_sensor)% 360
# to get values between 0 and 180°
# since rho values are symmetrical with the principal plane
relazi[relazi>180] = 360-relazi[relazi>180]
df.at[:, 'azi'] = relazi
awr = awr_process()

vza = 40
ws = 2
aot550 = 0.1
rho = awr.rho.rho.to_xarray()
rho_ = rho.interp(wind=ws,aot=aot550,wl=wl_common)

for name, raw in df.resample('1H'):
    print(name)

    #try:
    suff = name.__str__().replace(' ','_')
    N = len(raw.index)
    if not N:
        continue
     # ------------------
    # filtering
    # ------------------
    # daylight data
    ind = raw.sza < 80
    if not all(ind):
        continue

    #ind = awr.filtering(raw.Lt, raw.Lsky, raw.Ed)
    clean = raw[ind]
    Lt, Lsky, Ed, sza, azi = clean.Lt.values, clean.Lsky.values, clean.Ed.values, clean.sza.values, clean.azi.values
    sza_ = xr.DataArray(sza, dims='geom')
    azi_ = xr.DataArray(azi, dims='geom')

    # -----------------------------
    # data processing
    # -----------------------------
    rho_v = rho_.interp(vza = vza, sza = sza_,  azi = azi_).T

    clean['rho', ''] = rho_v.mean(axis=1)

    Lsurf = (rho_v * Lsky)
    Rrs = (Lt - Lsurf) / clean.Ed
    Rrs = Rrs.to_pandas()
    Rrs.index = clean.index
    Rrs.columns = pd.MultiIndex.from_product([ ['Rrs'], clean.Lt.columns,])
    clean = pd.concat([Rrs,clean],axis=1)
    clean.to_csv(opj(odir,'awr_L2_'+method+'_'+suff+'.csv'))


        # Rrs, rho = awr.process_wrapper(wl_common, raw, raw.sza, ws=ws, vza = [vza] * N, azi=raw.azi)
        #
        #
        #
        # awr.get_rho_values(raw.sza-90,[vza]*N, raw.azi.values, ws=[2])
        #
        # rho = awr.get_rho_mobley(awr.rhoM1999,raw.sza-90,[vza]*N, azi=raw.azi, ws=[2])[0,:,0]
