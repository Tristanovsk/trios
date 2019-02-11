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

import utils.utils as u
import utils.auxdata as ua

from process.process import *



class fit:
    def __init__(self,N=0,m=2):

        self.popt = np.full([N,m],np.nan)
        self.pcov = np.full([N,m,m],np.nan)

# import aeronet
# from config import *


# ------------------------------------------------
# above-water data files
dirfig = os.path.abspath('/DATA/OBS2CO/data/trios/fig')
dirout = os.path.abspath('/DATA/OBS2CO/data/trios/above_water')

iwrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/uw*idpr*.csv")

coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ.csv")[0]
coords = pd.read_csv(coordf, sep=';')
coords
# get idpr numbers
idprs = np.unique([re.findall(r'idpr(\d+)', x)[0] for x in iwrfiles])

################
# load aux data
iopw = ua.iopw()
iopw.load_iopw()
irr = ua.irradiance()
irr.load_F0()
#TODO check noise values (e.g., NEI from Trios), should it be spectral?
noise = 0.1
idpr = '178'
# loop over idpr
for idpr in idprs:#[-1:]:
    #    idpr=idprs[2]
    print(idpr)
    try:
        c = coords[coords.ID_prel == int(idpr)]  # .values[0]
        lat = c['Lat'].values[0]
        lon = c['Lon'].values[0]
        alt = 0  # c['Altitude']
        name = c['ID_lac'].values[0]

        dff = pd.DataFrame()

        # -----------------------------------------------
        #   IWR processing
        # -----------------------------------------------
        iwr = u.iwr_data(idpr, iwrfiles)
        if iwr.file:
            df, wl_ = iwr.reader(lat, lon, alt)
            reflectance = iwr_process(df, wl_).process()
            df = pd.concat([df, reflectance], axis=1)

            df.to_csv(os.path.join(dirout, 'trios_iwr_' + name + '_idpr' + idpr + '.csv'))

        mean = df.groupby('rounded_depth').mean()
        median = df.groupby('rounded_depth').median()
        std = df.groupby('rounded_depth').std()
        q25 = df.groupby('rounded_depth').quantile(0.25)
        q75 = df.groupby('rounded_depth').quantile(0.75)


        N = len(wl_)
        x = mean.prof_Edz #- 0.56
        depth_ = np.linspace(0, x.max(), 200)
        res=fit(N)
        res_w,res_rw,res_med=fit(N),fit(N),fit(N)

        for idx, wl in enumerate(wl_[:-10]):
            aw, bbw = iopw.get_iopw(wl)
            F0 = irr.get_F0(wl)

            y = mean.Edz.iloc[:, idx]
            sigma = std.Edz.iloc[:, idx]
            sigma[sigma<noise]=noise
            sigma.fillna(np.inf,inplace=True)
            res.popt[idx,:], res.pcov[idx,...] = curve_fit(iwr_process.f_Edz, x, y, [1.1 * aw, 100], bounds=([aw, 0], [np.inf, F0]))
            res_w.popt[idx,:], res_w.pcov[idx,...] = curve_fit(iwr_process.f_logEdz, x, np.log(1+y), [1.1 * aw, 100],  bounds=([aw, 0], [np.inf, F0]))#, sigma=sigma, absolute_sigma=True
            #res_rw.append( curve_fit(iwr_process.f_Edz, x, y, [1.1 * aw, 100], sigma=sigma, absolute_sigma=False, bounds=([aw, 0], [np.inf, F0])))
            y = median.Edz.iloc[:, idx]
            res_med.popt[idx,:], res_med.pcov[idx,...] = curve_fit(iwr_process.f_Edz, x, y, [1.1 * aw, 100], bounds=([aw, 0], [np.inf, F0]))

        i=0

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
            Ed_sim = iwr_process.f_Edz(depth_, *res.popt[idx,:])
            ax.plot(depth_, Ed_sim, linestyle='-', c='black', label='mean')
            Ed_sim = iwr_process.f_Edz(depth_, *res_w.popt[idx,:])
            ax.plot(depth_, Ed_sim, linestyle=':', c='black', label='log-space')
            #Ed_sim = iwr_process.f_Edz(depth_, *res_rw[idx][0])
            #ax.plot(depth_, Ed_sim, linestyle=':', c='red', label='mean, relativ weighted')
            Ed_sim = iwr_process.f_Edz(depth_, *res_med.popt[idx,:])
            ax.plot(depth_, Ed_sim, linestyle='--', c='black', label='median')

            ax.semilogy()
            # ax.colorbar()
            ax.set_ylabel(r'$E_d\ ({mW\cdot m^{-2}\cdot nm^{-1}})$')
            ax.set_xlabel('Depth (m)')
            ax.set_title(r'$\lambda = $' + str(round(wl_[idx], 1)) + ' nm, Kd = ' +
                         str(round(res_w.popt[idx,0], 3)) + '$m^{-1}$, Ed0 =' + str(round(res_w.popt[idx,1], 1)), fontsize=12)
        ax.legend(loc='best', frameon=False)
        fig.suptitle('trios_iwr ' + name + ' idpr' + idpr, fontsize=16)
        # fig.savefig(os.path.join(dirfig,'trios_iw_idpr' + idpr + '_'+name+'.png'),dpi=600)
        fig.savefig(os.path.join(dirfig, 'trios_iw_idpr' + idpr + '_' + name + '.pdf'))
        plt.close()


        def add_curve(ax, x, mean, std, c='red', label=''):
            ax.plot(x, mean, linestyle='solid', c=c, lw=2.5,
            alpha=0.8, label=label)
            ax.fill_between(x,
                    mean - std,
                    mean + std, alpha=0.35, color=c)

        mpl.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.29)

        iparam = 0
        ax = axs.flat[iparam]
        y, std = res.popt[:,iparam], res.pcov[:,iparam,iparam]
        add_curve(ax, wl_, y, std, c='red', label='mean')
        y, std = res_w.popt[:,iparam], res_w.pcov[:,iparam,iparam]
        add_curve(ax, wl_, y, std, c='orange', label='log-space')
        y, std = res_med.popt[:,iparam], res_med.pcov[:,iparam,iparam]
        add_curve(ax, wl_, y, std, c='green', label='median')
        ax.set_ylabel(r'$K_{d}\ ({m^{-1}})$')
        ax.set_xlabel(r'Wavelength (nm)')

        iparam = 1
        ax = axs.flat[iparam]
        y, std = res.popt[:,iparam], np.sqrt(res.pcov[:,iparam,iparam])
        add_curve(ax, wl_, y, std, c='red', label='mean')
        y, std = res_w.popt[:,iparam], np.sqrt(res_w.pcov[:,iparam,iparam])
        add_curve(ax, wl_, y, std, c='orange', label='log-space')
        y, std = res_med.popt[:,iparam], np.sqrt(res_med.pcov[:,iparam,iparam])
        add_curve(ax, wl_, y, std, c='green', label='median')
        add_curve(ax, wl_, mean.Ed.mean(), mean.Ed.std(), c='black', label='Es')

        ax.legend(loc='best', frameon=False)
        ax.set_ylabel(r'$E_{d}(0^{-})\ ({mW\cdot m^{-2}\cdot nm^{-1}})$')
        ax.set_xlabel(r'Wavelength (nm)')
        fig.suptitle('trios_iwr ' + name + ' idpr' + idpr, fontsize=16)
        fig.savefig(os.path.join(dirfig, 'trios_iwr_l2_idpr' + idpr + '_'+name+ '.pdf'))
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
