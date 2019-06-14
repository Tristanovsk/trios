import os, sys
import pandas as pd
import numpy as np
import glob
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmocean
import plotly
import plotly.graph_objs as go
import scipy.optimize as so

import trios.utils.utils as u
import trios.utils.auxdata as ua
from trios.process import *


class fit:
    def __init__(self, N=0, m=2):
        self.popt = np.full([N, m], np.nan)
        self.pcov = np.full([N, m, m], np.nan)


# import aeronet
# from config import *


# ------------------------------------------------
# above-water data files
dirfig = os.path.abspath('/DATA/OBS2CO/data/trios/fig')
odir = os.path.abspath('/DATA/OBS2CO/data/trios/in_water')

iwrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/uw*idpr*.csv")

coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ.csv")[0]
coords = pd.read_csv(coordf, sep=';')
coords['Date_prel'] = pd.to_datetime(coords['Date_prel'])
# get idpr numbers
idprs = np.unique([re.findall(r'idpr(\d+)', x)[0] for x in iwrfiles])


def add_curve(ax, x, mean, std, c='red', label=''):
    ax.plot(x, mean, linestyle='solid', c=c, lw=2.5,
            alpha=0.8, label=label)
    ax.fill_between(x,
                    mean - std,
                    mean + std, alpha=0.35, color=c)


################
# load aux data
iopw = ua.iopw()
iopw.load_iopw()
irr = ua.irradiance()
irr.load_F0()
# TODO check noise values (e.g., NEI from Trios), should it be spectral?
noise = 0.1
idpr = '177'
# loop over idpr
for idpr in idprs:  # [-1:]:
    #    idpr=idprs[2]
    print(idpr)
    try:
        c = coords[coords.ID_prel == int(idpr)]  # .values[0]
        date = c['Date_prel'].dt.strftime('%Y%m%d')
        lat = c['Lat'].values[0]
        lon = c['Lon'].values[0]
        alt = 0  # c['Altitude']
        name = c['ID_lac'].values[0]

        ofile = os.path.join(odir, 'Rrs_swr_' + date.values[0] + '_idpr' + idpr + '_' + name + '.csv')
        header = c.stack()
        header.index = header.index.droplevel()
        header.to_csv(ofile, header=False)

        dff = pd.DataFrame()

        # -----------------------------------------------
        #   IWR processing
        # -----------------------------------------------
        iwr = u.iwr_data(idpr, iwrfiles)
        if iwr.file:
            df, wl_ = iwr.reader(lat, lon, alt) # depth data already corrected for sensor position  , delta_Lu_depth=0.07,delta_Edz_depth=-0.28)
        else:
            continue

        # ---------------------------
        # Data formatting
        # ---------------------------
        # set pandas dataframes with data and parameters
        df['rounded_depth', ''] = df.prof_Edz.round(1)
        mean = df.groupby('rounded_depth').mean()
        median = df.groupby('rounded_depth').median()
        std = df.groupby('rounded_depth').std()
        q25 = df.groupby('rounded_depth').quantile(0.25)
        q75 = df.groupby('rounded_depth').quantile(0.75)
        x = mean.prof_Edz
        depth_ = np.linspace(0, x.max(), 200)

        # ---------------------------
        # Data processing
        # ---------------------------
        # load process object
        process = iwr_process(wl=wl_)

        # process data
        res = process.process(mean, std)
        res_w = process.process(mean, std, mode='log')
        res_med = process.process(median, std)

        # ---------------------------
        # Plotting section
        # ---------------------------
        i = 0

        mpl.rcParams.update({'font.size': 12})

        fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(16, 10))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.25)
        for idx in (28, 37, 51, 71, 91, 105, 130, 140, 170):
            ax = axs.flat[i]
            i += 1
            y = mean.Edz.iloc[:, idx]
            sigma = std.Edz.iloc[:, idx]
            yerr = [median.Edz.iloc[:, idx] - q25.Edz.iloc[:, idx], q75.Edz.iloc[:, idx] - median.Edz.iloc[:, idx]]
            ax.errorbar(x, y, yerr=yerr, ecolor='lightgray', fmt='none', elinewidth=3, capsize=0)
            ax.scatter(x, y,
                       # s = std.Edz.iloc[:,50]/std.Edz.iloc[:,50].mean()+100,
                       c=mean.Ed.iloc[:, idx],
                       alpha=0.6, cmap=cmocean.cm.thermal, label=None
                       )
            Ed_sim = iwr_process().f_Edz(depth_, *res.popt[idx, :])
            ax.plot(depth_, Ed_sim, linestyle='-', c='black', label='mean')
            Ed_sim = iwr_process().f_Edz(depth_, *res_w.popt[idx, :])
            ax.plot(depth_, Ed_sim, linestyle=':', c='black', label='log-space')
            # Ed_sim = iwr_process.f_Edz(depth_, *res_rw[idx][0])
            # ax.plot(depth_, Ed_sim, linestyle=':', c='red', label='mean, relativ weighted')
            Ed_sim = iwr_process().f_Edz(depth_, *res_med.popt[idx, :])
            ax.plot(depth_, Ed_sim, linestyle='--', c='black', label='median')

            ax.semilogy()
            # ax.colorbar()
            ax.set_ylabel(r'$E_d\ ({mW\cdot m^{-2}\cdot nm^{-1}})$')
            ax.set_xlabel('Depth (m)')
            ax.set_title(r'$\lambda = $' + str(round(wl_[idx], 1)) + ' nm, Kd = ' +
                         str(round(res_w.popt[idx, 0], 3)) + '$m^{-1}$, Ed0 =' + str(round(res_w.popt[idx, 1], 1)),
                         fontsize=12)
        ax.legend(loc='best', frameon=False)
        fig.suptitle('trios_iwr ' + name + ' idpr' + idpr, fontsize=16)
        # fig.savefig(os.path.join(dirfig,'trios_iw_idpr' + idpr + '_'+name+'.png'),dpi=600)
        fig.savefig(os.path.join(dirfig, 'trios_iw_idpr' + idpr + '_' + name + '.pdf'))
        plt.close()

        mpl.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.29)

        iparam = 0
        ax = axs.flat[iparam]
        y, std = res.popt[:, iparam], res.pcov[:, iparam, iparam]
        add_curve(ax, wl_, y, std, c='red', label='mean')
        y, std = res_w.popt[:, iparam], res_w.pcov[:, iparam, iparam]
        add_curve(ax, wl_, y, std, c='orange', label='log-space')
        y, std = res_med.popt[:, iparam], res_med.pcov[:, iparam, iparam]
        add_curve(ax, wl_, y, std, c='green', label='median')
        ax.set_ylabel(r'$K_{d}\ ({m^{-1}})$')
        ax.set_xlabel(r'Wavelength (nm)')

        iparam = 1
        ax = axs.flat[iparam]
        y, std = res.popt[:, iparam], np.sqrt(res.pcov[:, iparam, iparam])
        add_curve(ax, wl_, y, std, c='red', label='mean')
        y, std = res_w.popt[:, iparam], np.sqrt(res_w.pcov[:, iparam, iparam])
        add_curve(ax, wl_, y, std, c='orange', label='log-space')
        y, std = res_med.popt[:, iparam], np.sqrt(res_med.pcov[:, iparam, iparam])
        add_curve(ax, wl_, y, std, c='green', label='median')
        add_curve(ax, wl_, mean.Ed.mean(), mean.Ed.std(), c='black', label='Es')

        ax.legend(loc='best', frameon=False)
        ax.set_ylabel(r'$E_{d}(0^{-})\ ({mW\cdot m^{-2}\cdot nm^{-1}})$')
        ax.set_xlabel(r'Wavelength (nm)')
        fig.suptitle('trios_iwr ' + name + ' idpr' + idpr, fontsize=16)
        fig.savefig(os.path.join(dirfig, 'trios_iwr_l2_idpr' + idpr + '_' + name + '.pdf'))
        plt.close()

        # fig, ax = plt.subplots()
        # N=df.Edz.shape[1]
        # for wl, group in df.Edz.iloc[:,range(0,N,10)].iteritems():
        #
        #     ax.scatter(df.prof_Edz, group, marker='o',  c=cmocean.cm.thermal(wl/850*255))
        #
        # ax.set_ylim((0, 500000))

        #
        # writing output file
        # dff.to_csv(os.path.join(dirout, 'trios_awr_' + name + '_idpr' + idpr + '.csv'))

    except:
        plt.close()
        continue
