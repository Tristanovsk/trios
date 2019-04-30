import os
import pandas as pd
import numpy as np
import glob
import re
import datetime

from scipy.interpolate import interp1d

from utils.sunposition import sunpos
import utils.utils as u
import utils.auxdata as ua
from trios.process import *

plot=True #False
odir = os.path.abspath('/DATA/OBS2CO/data/trios/above_water')
dirfig = os.path.join(odir,'fig')

awrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/aw*idpr*.csv")

coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ_test.csv")[0]
coords = pd.read_csv(coordf, sep=';')
coords['Date_prel']=pd.to_datetime(coords['Date_prel'])
coords['h_debut']=coords['Date_prel'] + pd.to_timedelta(coords['h_debut'])
coords['h_fin']=coords['Date_prel'] + pd.to_timedelta(coords['h_fin'])
# TODO warning: time is set as start_time + 15 minutes (can be more accurate)
coords['Date_prel']= coords['h_debut']+datetime.timedelta(minutes = 15)

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
idprs = np.unique([re.findall(r'idpr(\d+)', x)[0] for x in awrfiles])
#idprs = np.array(['170'])
# loop over idpr
for idpr in idprs:
    c = coords[coords.ID_prel == int(idpr)]  # .values[0]
    date=c['Date_prel'].dt.strftime('%Y%m%d')
    lat = c['Lat'].values[0]
    lon = c['Lon'].values[0]
    alt = 0  # c['Altitude']
    name = c['ID_lac'].values[0]

    ofile=os.path.join(odir,'Rrs_awr_'+date.values[0]+'_idpr'+idpr+'_'+name+'.csv')
    header = c.stack(dropna=False)
    header.index = header.index.droplevel()
    header.to_csv(ofile, header=None)

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
        Rrs_h = (df.loc[:, 'Lt'] - rho_h * df.loc[:, 'Lsky']) / df.loc[:, 'Ed']

        Rrs_stat = Rrs_h.describe()
        #Rrs_stat.columns=Rrs_stat.columns.droplevel()
        Rrs_stat = Rrs_stat.T
        Rrs_stat.to_csv(ofile,mode='a')

        if plot:
            import matplotlib.pyplot as plt
            plt.rcParams.update({'font.size': 18})
            # -------------------------------
            # for Mobley values :
            rho15 = awr.get_rho_mobley(awr.rhoM2015, [df.sza.mean()], [vza], [azi], [ws])[0]
            rho99 = awr.get_rho_mobley(awr.rhoM1999, [df.sza.mean()], [vza], [azi], [ws])[0]

            Rrs15 = (df.loc[:, 'Lt'] - rho15 * df.loc[:, 'Lsky']) / df.loc[:, 'Ed']
            Rrs99 = (df.loc[:, 'Lt'] - rho99 * df.loc[:, 'Lsky']) / df.loc[:, 'Ed']

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))
            fig.subplots_adjust(left=0.16, right=0.9, hspace=.5, wspace=0.65)
            add_curve(ax, wl, Rrs15.transpose().mean(axis=1), Rrs15.transpose().std(axis=1),
                  label='M2015 (' + str('%.4f' %rho15) + ')')
            add_curve(ax, wl, Rrs99.transpose().mean(axis=1), Rrs99.transpose().std(axis=1), c='orange',
                  label='M1999(' + str('%.4f' %rho99) + ')')
            add_curve(ax, wl, Rrs_h.transpose().mean(axis=1), Rrs_h.transpose().std(axis=1), c='grey',
                  label='h(' + str('%.4f' %rho_h.mean()) + ')')
            info = str(header.Meteo_obs)+' / '+str(header.TRIOS_obs)
            info = info.replace('nan','')
            ax.set_title(info+'\n azi=' + str(azi) + ', vza=' + str(vza) + ', sza=' + str('%.2f' %sza))

            ax.legend(loc='best', frameon=False)
            ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
            ax.set_xlabel(r'$Wavelength (nm)$')
            fig.savefig(os.path.join(dirfig, 'trios_awr_'+date.values[0]+'_'+ name + '_idpr' + idpr + '.png'))
            plt.close()





