import glob
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
from scipy.interpolate import interp1d

from trios.utils.sunposition import sunpos
from trios.utils import utils as u
from trios.process import *

coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ.csv")[0]
coords = pd.read_csv(coordf, sep=';')
dirfig = os.path.abspath('/DATA/OBS2CO/data/trios/fig')

awrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/2018/aw*idpr*.csv")
# awrfiles = glob.glob("/DATA/OBS2CO/data/trios/test_setup/raw/aw*idpr*.csv")
swrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/2018/Lu0*idpr*.csv")

iopw = ua.iopw()
iopw.load_iopw()


def add_curve(ax, x, mean, std=None, c='red', label='',**kwargs):
    ax.plot(x, mean, linestyle='solid', c=c, lw=2.5,
            alpha=0.8, label=label,*kwargs)
    if np.any(std):
        ax.fill_between(x,
                        mean - std,
                        mean + std, alpha=0.35, color=c)


idpr = '167'

# get idpr numbers
idprs = np.unique([re.findall(r'idpr(\d+)', x)[0] for x in swrfiles])
# idprs = np.array(['170'])
# loop over idpr
for idpr in idprs:
    c = coords[coords.ID_prel == int(idpr)]  # .values[0]
    lat = c['Lat'].values[0]
    lon = c['Lon'].values[0]
    alt = 0  # c['Altitude']
    name = c['ID_lac'].values[0]

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
        raw = pd.merge_asof(Lt, Ed, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                            direction="nearest")
        raw = pd.merge_asof(raw, Lsky, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                            direction="nearest")

        # add solar angle data and idpr
        # compute solar angle (mean between fisrt and last aqcuisition time
        raw['sza', ''] = np.nan
        for index, row in raw.iterrows():
            # print index
            sza = sunpos(index, lat, lon, alt)[1]
            raw.at[index, 'sza'] = sza

        # ------------------
        # filtering
        # ------------------
        ind = awr.filtering(raw.Lt, raw.Lsky, raw.Ed)
        clean = raw[ind]
        Lt, Lsky, Ed, sza = clean.Lt.values, clean.Lsky.values, clean.Ed.values, clean.sza.values

        # -----------------------------
        # data processing
        # -----------------------------
        Rrs99, rho99 = awr.process_wrapper(wl, clean, clean.sza, ws=ws, azi=azi)
        Rrs15, rho15 = awr.process_wrapper(wl, clean, clean.sza, ws=ws, azi=azi, method='M15')
        Rrs_h, rho_h = awr.process_wrapper(wl, clean, clean.sza, ws=ws, azi=azi, method='osoaa')
        Rrs_opt, Rrs_opt_std = awr.process_optimization(wl, Lt, Lsky, Ed, sza, azi=azi)
        wl = Rrs99.T.index.get_level_values(1)
        date = Rrs99.index.get_level_values(0).date[0].__str__()

        # ------------------
        # plotting
        # ------------------
        Ltm = Lt.mean(axis=0)
        Edm = Ed.mean(axis=0)

        mpl.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 12))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.45)

        # ---- Ed
        ax = axs[0,0]
        add_curve(ax, wl, Ed.mean(axis=0),
                  label=r'$L_{sky}$',c='red') # just to put the two labels
        add_curve(ax, wl, Ed.mean(axis=0), Ed.std(axis=0),
                  label=r'$E_s$',c='black')
        ax.set_ylabel(r'$E_{d}(0^{+})$')

        # ---- Lsky
        ax2 = ax.twinx()
        add_curve(ax2, wl, Lsky.mean(axis=0), Lsky.std(axis=0),
                  label=r'$L_{sky}$',c='red')
        ax2.set_ylabel(r'$L_{sky}$',color='r')
        ax2.tick_params('y', colors='r')
        ax.set_xlabel(r'Wavelength (nm)')
        ax.legend(loc='best', frameon=False)


        # ---- Lt vs Lsurf
        ax = axs[0,1]
        add_curve(ax, wl, Lt.mean(axis=0), Lt.std(axis=0),
                  label=r'$L_t$',c='black')
        add_curve(ax, wl, Lsky.mean(axis=0)*rho15, Lsky.std(axis=0)*rho15,
                  label='M2015 (' + str(round(rho15,4)) + ')',c='violet')
        add_curve(ax, wl, Lsky.mean(axis=0)*rho99, Lsky.std(axis=0)*rho99, c='orange',
                  label='M1999(' + str(round(rho99,4)) + ')')
        add_curve(ax, wl, Lsky.mean(axis=0)*rho_h, Lsky.std(axis=0)*rho_h, c='green',
                  label='h(' + str(round(rho_h.mean(),4)) + ')')

        ax.set_ylabel(r'$L_t\ or L_{surf}$')
        ax.set_xlabel(r'Wavelength (nm)')

        # ---- Proportion o(Lt - Lsurf ) /Lt
        ax = axs[0,2]


        add_curve(ax, wl, Lsky.mean(axis=0)*rho15/Ltm, Lsky.std(axis=0)*rho15,
                  label='M2015 (' + str(round(rho15,4)) + ')',c='violet')
        add_curve(ax, wl, Lsky.mean(axis=0)*rho99/Ltm, Lsky.std(axis=0)*rho99, c='orange',
                  label='M1999(' + str(round(rho99,4)) + ')')
        add_curve(ax, wl, Lsky.mean(axis=0)*rho_h/Ltm, Lsky.std(axis=0)*rho_h, c='green',
                  label='h(' + str(round(rho_h.mean(),4)) + ')')

        ax.set_ylabel(r'$L_{surf}/L_t$')
        ax.set_xlabel(r'Wavelength (nm)')

        # ---- Lw
        ax = axs[1,0]
        add_curve(ax, wl, Rrs15.mean(axis=0)*Edm, Rrs15.std(axis=0)*Edm,
                  label='M2015 (' + str(round(rho15,4)) + ')',c='violet')
        add_curve(ax, wl, Rrs99.mean(axis=0)*Edm, Rrs99.std(axis=0)*Edm, c='orange',
                  label='M1999(' + str(round(rho99,4)) + ')')
        add_curve(ax, wl, Rrs_h.mean(axis=0)*Edm, Rrs_h.std(axis=0)*Edm, c='green',
                  label='h(' + str(round(rho_h.mean(),4)) + ')')
        add_curve(ax, wl, Rrs_opt*Edm, Rrs_opt_std*Edm, c='blue',
                  label='Optimization')
        ax.set_ylabel(r'$L_{w}\  (sr^{-1})$')
        ax.set_xlabel(r'Wavelength (nm)')

        # ---- Rrs
        ax = axs[1,1]
        add_curve(ax, wl_swr, Rrs_swr.transpose().mean(axis=1), Rrs_swr.transpose().std(axis=1), label='swr', c='black')
        add_curve(ax, wl, Rrs15.transpose().mean(axis=1), Rrs15.transpose().std(axis=1),
                  label='M2015 (' + str(round(rho15,4)) + ')',c='violet')
        add_curve(ax, wl, Rrs99.transpose().mean(axis=1), Rrs99.transpose().std(axis=1), c='orange',
                  label='M1999(' + str(round(rho99,4)) + ')')
        add_curve(ax, wl, Rrs_h.transpose().mean(axis=1), Rrs_h.transpose().std(axis=1), c='green',
                  label='h(' + str(round(rho_h.mean(),4)) + ')')
        add_curve(ax, wl, Rrs_opt, Rrs_opt_std, c='blue',
                  label='Optimization')
        ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
        ax.set_xlabel(r'Wavelength (nm)')
        ax.set_title('azi=' + str(azi) + ', vza=' + str(vza) + ', sza=' + str(round(sza.mean(), 2)))

        # ---- delta Rrs
        ax = axs[1,2]
        Rrs_swr_ = interp1d(wl_swr, Rrs_swr.transpose().mean(axis=1), fill_value='extrapolate')(wl)
        Rrs_swr_[wl > 850]=np.nan
        add_curve(ax, wl, (Rrs15.mean(axis=0)-Rrs_swr_)/Rrs_swr_,
                  label='M2015 (' + str(round(rho15,4)) + ')',c='violet')
        add_curve(ax, wl, (Rrs99.mean(axis=0)-Rrs_swr_)/Rrs_swr_, c='orange',
                  label='M1999(' + str(round(rho99,4)) + ')')
        add_curve(ax, wl, (Rrs_h.mean(axis=0)-Rrs_swr_)/Rrs_swr_, c='green',
                  label='h(' + str(round(rho_h.mean(),4)) + ')')
        add_curve(ax, wl, (Rrs_opt-Rrs_swr_)/Rrs_swr_, c='blue',
                  label='Optimization')
        ax.set_ylabel(r'$\Delta^{rel} R_{rs} $')
        ax.set_xlabel(r'Wavelength (nm)')
        ax.legend(loc='best', frameon=False)

        fig.suptitle('trios_awr ' + name + ' idpr' + idpr, fontsize=16)
        fig.savefig(os.path.join(dirfig, 'trios_awr_' + name + '_idpr' + idpr + '.png'))
        plt.close()
