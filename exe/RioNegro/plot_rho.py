import os
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

plt.ioff()
mpl.rcParams.update({'font.size': 18})
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from scipy import odr as odr
from scipy import stats as scistats

import RTxploitation as rt
iopw = rt.auxdata.iopw().get_iopw

opj= os.path.join
idir = os.path.abspath('/DATA/OBS2CO/data/rogerio')
figdir = opj(idir,'fig')
file = opj(idir,'data/Simulation_Rrs_OSOAA_TSS_COD_aCDOM.xlsx')

data = pd.read_excel(file)
data.sort_values(by="TSS (mg L)", inplace=True)
#data = data.set_index(['Station', 'Date'])

idir = '/DATA/OBS2CO/data/rogerio/data/L2'
stat_file= os.path.join(idir, 'all/rho_osoaa_stat.csv')
df = pd.read_csv(stat_file, index_col=[0], parse_dates=True)
df.sort_values(['date','ID','method','wl'],inplace=True)

df_ave=df.groupby("wl").mean()
df_q05 = df[["wl","0.5"]].groupby("wl").quantile(q=0.05)
df_q95 = df[["wl","0.5"]].groupby("wl").quantile(0.95)

plt.figure(figsize=(10,6))
plt.plot(df_ave.index, df_ave["0.5"], '--k', lw=1.5)
#plt.scatter(df.wl,df['0.5'])
plt.fill_between(df_q05.index, df_q05["0.5"], df_q95["0.5"], alpha=.25, facecolor="r")
plt.ylabel('Reflectance factor')
plt.xlabel('Wavelength (nm)')
plt.tight_layout()
plt.savefig(opj(figdir,'rho_values.png'),dpi=300)