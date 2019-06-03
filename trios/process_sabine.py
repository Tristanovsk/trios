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
import scipy.optimize as so

from trios.process import *
from trios.utils.sunposition import sunpos

dir = '/DATA/OBS2CO/data/sabine/data/raw'
odir = '/DATA/OBS2CO/data/sabine/data/awr'
dirfig = '/DATA/OBS2CO/data/sabine/fig'
file = os.path.join(dir, 'Ed_Ld_Lu_all_highroc_cruise_oslofjord.dat')

dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')

dff = pd.read_csv(file, sep="\t", date_parser=dateparse, index_col=2)

awr = awr_process()
awr.get_rho_values([30], [40], [135])

groupname = dff.columns.values[1]  # 'RT$station_area[n]'


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

    # -----------------------------
    # data formatting
    # -----------------------------
    raw_ = pd.DataFrame()
    for time, df in df_.groupby(df_.index):
        azi = [int(df.iloc[:, 0].unique()[0].split('_')[-1])]
        lat = df.LAT.values.mean()
        lon = df.LON.values.mean()
        ws = [df.WIND_SPEED.values.mean() / 2]
        # TODO check altitude
        alt = 0

        Lsky = get_param(df, 'Ld', 'Lsky')
        Lt = get_param(df, 'Lu', 'Lt')
        Ed = get_param(df, 'Ed', 'Ed')

        wl = Lsky.columns.levels[1].values

        # compute solar angle (mean between first and last acquisition time

        sza = sunpos(time, lat, lon, alt)[1]
        sza = np.repeat(sza, Lt.__len__())

        # Rrs, rho = awr.process(wl, Lt, Lsky, Ed, sza, ws=ws, azi=azi)

        # keep metadata from de and reshape dataframe
        df.columns = pd.MultiIndex.from_product([df.columns.values, ['']], names=['param', 'wl'])

        # ['Rrs(awr)'] * Rrs.shape[1]
        raw = pd.concat([df.iloc[:, 0:10], Ed, Lsky, Lt], axis=1)
        raw['sza'] = sza.mean()
        raw['azi'] = azi[0]
        # raw['rho']=rho.mean()
        raw_ = pd.concat([raw_, raw], axis=0)


    def filtering(Lt, Lsky, Ed, **kargs):
        '''

        :param Lt:
        :param Lsky:
        :param Ed:
        :param kargs:
        :return:
        '''

        ind_Ed, notused = calc.spectra_median_filter(Ed, kargs)
        ind_sky, notused = calc.spectra_median_filter(Lsky, kargs)
        ind = ind_Ed & ind_sky
        return ind


    def process_wrapper(wl, df, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):

        print(sza, vza, azi, ws, aot, method)
        Rrs, rho = process(wl, df.Lt, df.Lsky.values, df.Ed.values, sza, vza, azi, ws, aot, method)

        Rrs.columns = pd.MultiIndex.from_product([['Rrs'], wl], names=['param', 'wl'])

        return Rrs, rho


    def process(wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):
        '''

        :param wl:
        :param Lt:
        :param Lsky:
        :param Ed:
        :param sza:
        :param vza:
        :param azi:
        :param ws:
        :param aot:
        :param method: 'M99, 'M15, 'osoaa'
        :return:
        '''

        # -----------------------------
        # standard data processing
        # -----------------------------

        if method == 'osoaa':
            rho = awr.get_rho_values(np.median(sza), vza, azi, wl=wl, ws=ws, aot=aot)
        elif method == 'M99':
            rho = awr.get_rho_mobley(awr.rhoM1999, [np.median(sza)], [vza], [azi], [ws])
        elif method == 'M15':
            rho = awr.get_rho_mobley(awr.rhoM2015, [np.median(sza)], [vza], [azi], [ws])
        else:
            return print('ERROR: no method for rho factor')

        Rrs = (Lt - rho * Lsky) / Ed
        # Rrs.columns = pd.MultiIndex.from_product([['Rrs'], wl], names=['param', 'wl'])

        return Rrs, rho.mean()


    def cost_func(x, param, meas, Rrs_bar):
        sza, vza, azi = param
        Lt, Lsky, Ed = meas
        ws = x

        rho = awr.get_rho_mobley(awr.rhoM1999, [sza], [vza], [azi], [ws])

        Rrs = (Lt - rho * Lsky) / Ed

        return Rrs - Rrs_bar

        # x_ave = x.mean()
        # return np.sum(np.abs(x - x_ave))


    def process_optimization(wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):

        plt.figure()

        # ------------------------------
        # initialization of mean/median values
        # ------------------------------
        Rrs, rho = process(wl, Lt, Lsky, Ed, sza, ws=ws, azi=azi)
        Rrs_bar = Rrs.mean(axis=0)

        # -----------------------------
        # non-linear optimization
        # -----------------------------

        for j in range(10):
            x_est = []
            res = []
            Rrs_est = []
            rho_est = []
            plt.plot(wl, Rrs_bar)
            for i in range(len(Lt)):
                geom = [sza[i], vza, azi]
                meas = [Lt[i], Lsky[i], Ed[i]]
                x0 = ws
                res_lsq = so.least_squares(cost_func, x0, args=(geom, meas, Rrs_bar))
                res.append(res_lsq)
                x_est.append(res_lsq.x[0])
                Rrs, rho = process(wl, Lt[i], Lsky[i], Ed[i], sza[i], ws=res_lsq.x[0], azi=azi)
                Rrs_est.append(Rrs)
                rho_est.append(rho)
                print(res_lsq.x, res_lsq.cost)
            Rrs_bar = np.mean(Rrs_est,axis=0)
            Rrs_std = np.std(Rrs_est,axis=0)
        return Rrs_bar, Rrs_std


    # ------------------
    # filtering
    # ------------------
    ind = filtering(raw_.Lt, raw_.Lsky, raw_.Ed)
    clean = raw_[ind]
    Lt, Lsky, Ed, sza = clean.Lt.values, clean.Lsky.values, clean.Ed.values, clean.sza.values
    # -----------------------------
    # data processing
    # -----------------------------

    Rrs, rho = process_wrapper(wl, raw_, raw_.sza, ws=ws, azi=azi)
    Rrs_15, rho = process_wrapper(wl, raw_, raw_.sza, ws=ws, azi=azi, method='M15')
    Rrs_h, rho = process_wrapper(wl, raw_, raw_.sza, ws=ws, azi=azi, method='osoaa')
    Rrs_opt, Rrs_opt_std = process_optimization(wl, Lt, Lsky, Ed, sza, azi=azi)

    wl = Rrs.T.index.get_level_values(1)
    date = Rrs.index.get_level_values(0).date[0].__str__()


    def add_curve(ax, x, mean, std, c='red', label=''):
        ax.plot(x, mean, linestyle='solid', c=c, lw=2.5,
                alpha=0.8, label=label)
        ax.fill_between(x,
                        mean - std,
                        mean + std, alpha=0.35, color=c)


    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))
    fig.subplots_adjust(left=0.16, right=0.9, hspace=.5, wspace=0.65)
    add_curve(ax, wl, Rrs[~ind].transpose().mean(axis=1), Rrs[~ind].transpose().std(axis=1), label='QC-excluded',
              c='red')
    add_curve(ax, wl, Rrs[ind].transpose().mean(axis=1), Rrs[ind].transpose().std(axis=1), label='M99, QC-passed',
              c='black')
    add_curve(ax, wl, Rrs_h[ind].transpose().mean(axis=1), Rrs_h[ind].transpose().std(axis=1), label='H, QC-passed',
              c='green')
    add_curve(ax, wl, Rrs_15[ind].transpose().mean(axis=1), Rrs_15[ind].transpose().std(axis=1), label='M15, QC-passed',
              c='violet')
    add_curve(ax, wl, Rrs_opt, Rrs_opt_std, label='Optimization',
              c='blue')
    ax.legend(loc='best', frameon=False)
    ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
    ax.set_xlabel(r'$Wavelength\ (nm)$')
    fig.savefig(os.path.join(dirfig, 'trios_awr_' + date + '_' + group + '.png'))
    plt.close()

    raw_df = pd.concat([raw_df, raw_], axis=0)

    Rrs_df = pd.concat([Rrs_df, Rrs], axis=1)
    Rrs_df.to_csv(ofile)

if False:
    dff = Rrs  # _df
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
