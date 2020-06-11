
import os, sys
import pandas as pd
import numpy as np
import glob
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d

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

iwrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/2018/uw*idpr*.csv")

coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ.csv")[0]
coords = pd.read_csv(coordf, sep=';')
coords['Date_prel'] = pd.to_datetime(coords['Date_prel'])
# get idpr numbers
idprs = np.unique([re.findall(r'idpr(\d+)', x)[0] for x in iwrfiles])


def add_curve(ax, x, mean, std=None, c='red', label='',**kwargs):
    ax.plot(x, mean, linestyle='solid', c=c, lw=2.5,
            alpha=0.8, label=label,*kwargs)
    if np.any(std):
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
idpr = '150'

idx_list_for_plot=(28, 37, 51, 71, 91, 105, 130, 140, 170)
#idx_list_for_plot=(165,170,175,180,185,190,195,200,205)
# loop over idpr
for idpr in idprs[0:]:
    #    idpr=idprs[2]
    print(idpr)
    try:
        c = coords[coords.ID_prel == int(idpr)]  # .values[0]
        date = c['Date_prel'].dt.strftime('%Y%m%d')
        lat = c['Lat'].values[0]
        lon = c['Lon'].values[0]
        alt = 0  # c['Altitude']
        name = c['ID_lac'].values[0]

        # -----------------------------------------------
        #   write data header
        # -----------------------------------------------
        ofile = os.path.join(odir, 'iwr_' + date.values[0] + '_idpr' + idpr + '_' + name + '.csv')
        header = c.stack()
        header.index = header.index.droplevel()
        header.to_csv(ofile, header=False)

        dff = pd.DataFrame()

        # -----------------------------------------------
        #   IWR processing
        # -----------------------------------------------
        iwr = u.iwr_data(idpr, iwrfiles)
        # if iwr.file:
        df, wl_ = iwr.reader(lat, lon,
                             alt)  # depth data already corrected for sensor position  , delta_Lu_depth=0.07,delta_Edz_depth=-0.28)
        # else:
        #    continue

        # ---------------------------
        # Data formatting
        # ---------------------------
        # set pandas dataframes with data and parameters
        df['rounded_depth', ''] = df.prof_Edz.round(1)

        # data filtering
        df[df.rounded_depth < -0.05] = np.nan
        #df[df == 0] = np.nan

        mean = df.groupby('rounded_depth').mean()
        median = df.groupby('rounded_depth').median()
        std = df.groupby('rounded_depth').std()
        df_=df.drop(df.columns[df.dtypes=='object'],axis=1)
        q25 = df_.groupby('rounded_depth').quantile(0.25)
        q75 = df_.groupby('rounded_depth').quantile(0.75)

        # ---------------------------
        # Data processing
        # ---------------------------
        # load process object

        process = iwr_process(wl=wl_)

        # process data
        res = process.process(mean, std)
        res_w = process.process(mean, std, mode='log')
        res_med = process.process(median, std)
        res_lsq = process.process(mean, std, mode='lsq')

        Kd = res_lsq.popt[:, 0]
        Edm = res_lsq.popt[:, 1]
        Es = mean.Ed.mean()
        Lwm = res_lsq.popt[:, 3]
        rrs = Lwm / Edm
        Es_from_Edm = 0.96 * Edm
        Lwp = 0.541 * Lwm
        Rrs_Edp = Lwp / Es
        Rrs_Edm = Lwp / Es_from_Edm

        Rrs_df = pd.DataFrame(Rrs_Edp)
        Rrs_df.columns = ['Rrs']
        Rrs_df['Kd']=Kd

        # -----------------------------------------------
        #   write data header
        # -----------------------------------------------
        Rrs_df.to_csv(ofile,mode='a')


        # Open SWR files for comparison
        swr = pd.read_csv(
            '/DATA/OBS2CO/data/trios/surface_water/Rrs_swr_' + date.values[0] + '_idpr' + idpr + '_' + name + '.csv',
            skiprows=57)

        # ---------------------------
        # Plotting section
        # ---------------------------

        pdf = PdfPages(os.path.join(dirfig, 'trios_iwr_all_idpr' + idpr + '_' + name + '.pdf'))

        # -------------------------------------
        # Profile plotting
        # -------------------------------------
        # Ed profile
        i = 0
        x = mean.prof_Edz
        depth_ = np.linspace(0, x.max(), 200)

        mpl.rcParams.update({'font.size': 12})

        fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(16, 10))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.25)

        for idx in idx_list_for_plot:
            ax = axs.flat[i]
            i += 1
            y = mean.Edz.iloc[:, idx]
            sigma = std.Edz.iloc[:, idx]
            yerr = [median.Edz.iloc[:, idx] - q25.Edz.iloc[:, idx], q75.Edz.iloc[:, idx] - median.Edz.iloc[:, idx]]
            ax.errorbar(x, y, yerr=yerr, ecolor='lightgray', alpha=0.43, fmt='none', elinewidth=3, capsize=0, label='std')
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
        # fig.savefig(os.path.join(dirfig, 'trios_iwr_Ed_profile_idpr' + idpr + '_' + name + '.pdf'))
        pdf.savefig()
        plt.close()

        # Lu profile
        i = 0
        x = mean.prof_Luz
        depth_ = np.linspace(0, x.max(), 200)
        mpl.rcParams.update({'font.size': 12})

        fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(16, 10))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.25)

        for idx in idx_list_for_plot:
            ax = axs.flat[i]
            i += 1
            y = mean.Luz.iloc[:, idx]
            sigma = std.Luz.iloc[:, idx]
            yerr = [y - q25.Luz.iloc[:, idx], q75.Luz.iloc[:, idx] - y]
            yerr = sigma
            ax.errorbar(x, y, yerr=yerr, ecolor='lightgray', alpha=0.43, fmt='none', elinewidth=3, capsize=0, label='std')
            ax.scatter(x, y,
                       # s = std.Luz.iloc[:,50]/std.Luz.iloc[:,50].mean()+100,
                       c=mean.Ed.iloc[:, idx],
                       alpha=0.6, cmap=cmocean.cm.thermal, label=None
                       )
            Lu_sim = iwr_process().f_Lu(depth_, *res_lsq.popt[idx, -2:])
            ax.plot(depth_, Lu_sim, linestyle='-', c='black', label='mean')
            # Lu_sim = iwr_process().f_Luz(depth_, *res_w.popt[idx, :])
            # ax.plot(depth_, Lu_sim, linestyle=':', c='black', label='log-space')
            # # Lu_sim = iwr_process.f_Luz(depth_, *res_rw[idx][0])
            # # ax.plot(depth_, Lu_sim, linestyle=':', c='red', label='mean, relativ weighted')
            # Lu_sim = iwr_process().f_Luz(depth_, *res_med.popt[idx, :])
            # ax.plot(depth_, Lu_sim, linestyle='--', c='black', label='median')

            ax.semilogy()
            # ax.colorbar()
            ax.set_ylabel(r'$L_u\ ({mW\cdot m^{-2}\cdot sr^{-1}\cdot nm^{-1}})$')
            ax.set_xlabel('Depth (m)')
            ax.set_title(r'$\lambda = $' + str(round(wl_[idx], 1)) + r'$ nm, K_{Lu} = $' +
                         str(round(res_lsq.popt[idx, 2], 3)) + '$m^{-1}$, Lu0 =' + str(round(res_lsq.popt[idx, 3], 1)),
                         fontsize=12)
        ax.legend(loc='best', frameon=False)
        fig.suptitle('trios_iwr ' + name + ' idpr' + idpr, fontsize=16)
        # fig.savefig(os.path.join(dirfig,'trios_iw_idpr' + idpr + '_'+name+'.png'),dpi=600)
        # fig.savefig(os.path.join(dirfig, 'trios_iwr_Lu_profile_idpr' + idpr + '_' + name + '.pdf'))
        pdf.savefig()

        plt.close()

        # -------------------------------------
        # Sprectra plotting
        # -------------------------------------
        mpl.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 12))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.29)

        iparam = 0
        ax = axs.flat[iparam]
        y, sigma = res.popt[:, iparam], res.pcov[:, iparam, iparam]
        add_curve(ax, wl_, y, sigma, c='red', label='mean')
        y, sigma = res_w.popt[:, iparam], res_w.pcov[:, iparam, iparam]
        add_curve(ax, wl_, y, sigma, c='orange', label='log-space')
        y, sigma = res_med.popt[:, iparam], res_med.pcov[:, iparam, iparam]
        add_curve(ax, wl_, y, sigma, c='green', label='median')
        y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam])
        add_curve(ax, wl_, y, sigma, c='blue', label='lsq')
        ax.set_ylabel(r'$K_{d}\ ({m^{-1}})$')
        ax.set_xlabel(r'Wavelength (nm)')

        iparam = 1
        ax = axs.flat[iparam]
        y, sigma = res.popt[:, iparam], np.sqrt(res.pcov[:, iparam, iparam])
        add_curve(ax, wl_, y, sigma, c='red', label='mean')
        y, sigma = res_w.popt[:, iparam], np.sqrt(res_w.pcov[:, iparam, iparam])
        add_curve(ax, wl_, y, sigma, c='orange', label='log-space')
        y, sigma = res_med.popt[:, iparam], np.sqrt(res_med.pcov[:, iparam, iparam])
        add_curve(ax, wl_, y, sigma, c='green', label='median')
        y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam])
        add_curve(ax, wl_, y, sigma, c='blue', label='lsq')
        add_curve(ax, wl_, mean.Ed.mean(), mean.Ed.std(), c='black', label='Es')
        ax.legend(loc='best', frameon=False)
        ax.set_ylabel(r'$E_{d}(0^{-})\ ({mW\cdot m^{-2}\cdot nm^{-1}})$')
        ax.set_xlabel(r'Wavelength (nm)')

        iparam = 2
        ax = axs.flat[iparam]
        y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam]) * 0
        sigma[sigma > 10 * np.nanmedian(sigma)] = np.nan
        add_curve(ax, wl_, y, sigma, c='blue', label='lsq')
        ax.set_ylabel(r'$K_{Lu}\ ({m^{-1}})$')
        ax.set_xlabel(r'Wavelength (nm)')

        iparam = 3
        ax = axs.flat[iparam]

        y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam]) * 0
        sigma[sigma > 10 * np.nanmedian(sigma)] = np.nan
        add_curve(ax, wl_, y, sigma, c='blue', label='lsq')
        ax.set_ylabel(r'$L_{w}(0^{-})\ ({mW\cdot m^{-2}\cdot sr^{-1}\cdot nm^{-1})$')
        ax.set_xlabel(r'Wavelength (nm)')

        fig.suptitle('trios_iwr ' + name + ' idpr' + idpr, fontsize=16)
        # fig.savefig(os.path.join(dirfig, 'trios_iwr_l2_idpr' + idpr + '_' + name + '.pdf'))
        pdf.savefig()

        plt.close()

        # Edm = res_lsq.popt[:,1]
        # Es = mean.Ed.mean()
        # Lwm = res_lsq.popt[:,3]
        # rrs = Lwm / Edm
        # Es_from_Edm = 0.96 * Edm
        # Lwp = 0.541 * Lwm
        # Rrs_Edp = Lwp / Es
        # Rrs_Edm = Lwp / Es_from_Edm

        mpl.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 12))
        fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.29)

        ax = axs.flat[0]
        add_curve(ax, wl_, Es / Edm, c='red', label='Es over Ed0-')
        ax.set_xlabel(r'Wavelength (nm)')
        ax.legend(loc='best', frameon=False)
        ax = axs.flat[1]
        add_curve(ax, wl_, rrs, c='red', label='rrs')
        ax.set_xlabel(r'Wavelength (nm)')
        ax.legend(loc='best', frameon=False)
        ax = axs.flat[2]
        add_curve(ax, wl_, Rrs_Edm, c='red', label='Rrs (from Ed0-)')
        add_curve(ax, wl_, Rrs_Edp, c='blue', label='Rrs (from Ed0+)')
        add_curve(ax, swr.wl, swr['mean'], swr['std'], c='black', label='Rrs (from SWR)')
        ax.set_xlabel(r'Wavelength (nm)')
        ax.legend(loc='best', frameon=False)
        ax = axs.flat[3]
        Rrs_swr = interp1d(swr.wl, swr['mean'], fill_value='extrapolate')(wl_)
        diff = Rrs_Edp - Rrs_swr
        add_curve(ax, wl_, diff, c='black', label='Rrs(IWR) - Rrs(SWR)')

        fig.suptitle('trios_iwr ' + name + ' idpr' + idpr, fontsize=16)
        pdf.savefig()
        plt.close()
        pdf.close()

    except:
        plt.close()
        continue
