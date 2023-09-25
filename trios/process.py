# standard lib
import pandas as pd
import xarray as xr
from scipy import interpolate, integrate
import scipy.optimize as so

# for plotting pruposes
import cmocean
import plotly.graph_objs as go
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})

# package
from trios.utils.utils import plot as up
import trios.utils.utils as uu
from trios.utils.utils import reshape as r
import trios.utils.auxdata as ua
from trios.config import *


class awr_process:
    def __init__(self, df=None, wl=None, aerosol = 'fine', name="", idpr=""):
        self.df = df
        self.wl = wl
        self.name = name
        self.idpr = idpr

        self.rhosoaa_fine_file = rhosoaa_fine_file
        self.rhosoaa_coarse_file = rhosoaa_coarse_file
        self.M1999_file = M1999_file
        self.M2015_file = M2015_file
        self.load_rho_lut()

        self.rho = self.rhosoaa_coarse
        if aerosol == 'fine':
            self.rho = self.rhosoaa_fine


    def load_rho_lut(self):
        self.rhosoaa_fine = pd.read_csv(self.rhosoaa_fine_file, index_col=[0, 1, 2, 3, 4, 5])
        self.rhosoaa_coarse = pd.read_csv(self.rhosoaa_coarse_file, index_col=[0, 1, 2, 3, 4, 5])
        self.rhoM1999 = pd.read_csv(self.M1999_file, skiprows=7, index_col=[0, 1, 2, 3])
        self.rhoM2015 = pd.read_csv(self.M2015_file, skiprows=8, index_col=[0, 1, 2, 3])
        self.rhoM1999.dropna(inplace=True)
        self.rhoM2015.dropna(inplace=True)

    def get_rho_values_old(self, sza, vza, azi, ws=[2], aot=[0.1], wl=[550], sunglint=False):
        '''
        Interpolate the rho factor values from tabulated data

        :param sza: solar zenith angle in deg, array-like
        :param vza: view zenith angle in deg, array-like
        :param azi: relative azimuth in deg (=0 when looking at Sun), array-like
        :param ws: wind speed, m/s,(based on Cox-Munk parametrization of surface roughness) array-like
        :param aot: aerosol optical thickness at 550 nm, array-like
        :param wl: wavelength in nm, array-like
        :param sunglint: add sunglint component in rho calculation if True
        :return:
        '''

        grid = self.rho.rho.index.levels

        # convert pandas dataframe into 6D array of the tabulated rho values for interpolation
        rhoname = 'rho'
        if sunglint:
            rhoname = 'rho_g'

        rho_6d = r().df2ndarray(self.rho, rhoname)

        rho_ = calc().spline_2d(grid[-2:], rho_6d, (azi, vza))

        rho_wl = calc().spline_4d(grid[:-2], rho_, (ws, aot, wl, sza))

        return rho_wl.squeeze()

    def get_rho_values(self, sza, vza, azi, ws=[2], aot=[0.1], wl=[550], sunglint=False):
        '''
        Interpolate the rho factor values from tabulated data

        :param sza: solar zenith angle in deg, array-like
        :param vza: view zenith angle in deg, array-like
        :param azi: relative azimuth in deg (=0 when looking at Sun), array-like
        :param ws: wind speed, m/s,(based on Cox-Munk parametrization of surface roughness) array-like
        :param aot: aerosol optical thickness at 550 nm, array-like
        :param wl: wavelength in nm, array-like
        :param sunglint: add sunglint component in rho calculation if True
        :return:
        '''
        rhoname = 'rho'
        if sunglint:
            rhoname = 'rho_g'
        rho = self.rho[rhoname].to_xarray()
        rho = rho.interp(vza=vza, azi = azi)
        rho = rho.interp(wind=ws,aot=aot)

        return rho.interp(sza = sza, wl=wl).squeeze().T.values

    def get_rho_mobley(self, sza, vza, azi, ws,method='M99'):
        '''
        Get the Mobley rho factor from cubic interpolation in the tabulated values


        :param sza:
        :param vza:
        :param azi:
        :param ws:
        :param method: LUT
        :return:
        '''
        # Warning SZA set to 90 if Sun below the horizon
        # sza[sza>90]=90
        if method == 'M99':
            rhodf = self.rhoM1999
        elif method == 'M15':
            rhodf = self.rhoM2015
        else:
            return print('Requested method does not exists, please try again!' )
        rhodf = rhodf.query('sza<75 & vza >0')
        rhodf.index = rhodf.index.remove_unused_levels()

        # grid {wind, sza, vza, azi}
        grid = rhodf.index.levels

        rho_ = r().df2ndarray(rhodf, 'rho')
        rho_mobley = calc().spline_4d(grid, rho_, (ws, sza, vza, azi))
        return rho_mobley

    def call_process(self, method='M99', ofile="", vza=40, azi=135, ws=2, aot=0.1, plot_file=""):

        wl = self.wl
        vza, azi, ws, aot = [vza], [azi], [ws], [aot]  # formatting for interpolation functions

        # ------------------
        # filtering
        # ------------------
        ind = self.filtering(self.df.Lt, self.df.Lsky, self.df.Ed)
        clean = self.df  # [ind]
        Lt, Lsky, Ed, sza = clean.Lt.values, clean.Lsky.values, clean.Ed.values, clean.sza.values

        # -----------------------------
        # data processing
        # -----------------------------
        Rrs, rho = self.process_wrapper(wl, clean, clean.sza, vza=vza, azi=azi, ws=ws, aot=aot, method=method)
        # if method == 'M99':
        #
        # elif method == 'M15':
        #     Rrs, rho = self.process_wrapper(wl, clean, clean.sza, vza=vza, azi=azi, ws=ws, aot=aot, method=method)
        # elif method == 'osoaa':
        #     Rrs, rho = self.process_wrapper(wl, clean, clean.sza, vza=vza, azi=azi, ws=ws, aot=aot, method=method)
        # elif method == 'temp_opt':
        #     Rrs, rho, Rrs_opt, Rrs_opt_std = self.process_optimization(wl, Lt, Lsky, Ed, sza, vza=vza, azi=azi)

        self.Rrs = Rrs.__deepcopy__()

        # self.Rrs['rho'] = rho
        # if method == 'temp_opt':
        #     rho = np.nanmedian(rho)

        if ofile:
            Rrs_df = pd.concat([self.df, self.Rrs,rho], axis=1)
            Rrs_df.to_csv(ofile)
            Rrs_stat = Rrs.describe()
            Rrs_stat.columns = Rrs_stat.columns.droplevel()
            Rrs_stat = Rrs_stat.T
            # Rrs_stat.to_csv(ofile)

        if plot_file:
            # ------------------
            # plotting envelope
            # ------------------
            Ltm = Lt.mean(axis=0)
            Edm = Ed.mean(axis=0)

            Lsurf = Lsky * rho.values

            def add_envelope(ax, wl, values, label='', **kwargs):
                up.add_curve(ax, wl, values.mean(axis=0), label=label, c='black', **kwargs)
                ax.fill_between(wl, np.quantile(values, 0.05, axis=0), np.quantile(values, 0.95, axis=0), alpha=0.25,
                                color='grey')
                ax.fill_between(wl, np.quantile(values, 0.25, axis=0), np.quantile(values, 0.75, axis=0), alpha=0.35,
                                color='red')

            mpl.rcParams.update({'font.size': 18})
            fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(20, 12))
            fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.45)

            # ---- Ed
            ax = axs[0, 0]
            add_envelope(ax, wl, Ed, label=r'$E_s$')
            # up.add_curve(ax, wl, Ed.mean(axis=0), Ed.std(axis=0), label=r'$E_s$', c='black')
            ax.set_ylabel(r'$E_{d}(0^{+})\ (mW\ m^{-2}\ nm^{-1})$')
            ax.set_xlabel(r'Wavelength (nm)')
            ax.legend(loc='best', frameon=False)

            # ---- Lsky
            ax = axs[0, 1]  # ax2 = ax.twinx()
            add_envelope(ax, wl, Lsky, label=r'$L_{sky}$')

            ax.set_ylabel(r'$L_{sky}\ (mW\ m^{-2}\ nm^{-1}\ sr^{-1})$')  # , color='r')
            # ax.tick_params('y', colors='r')
            ax.set_xlabel(r'Wavelength (nm)')
            ax.legend(loc='best', frameon=False)

            # ---- Lt vs Lsurf
            ax = axs[0, 2]
            add_envelope(ax, wl, Lt, label=r'$L_t$')
            # up.add_curve(ax, wl, Lt.mean(axis=0), Lt.std(axis=0), label=r'$L_t$', c='black')
            print('-----',rho)
            add_envelope(ax, wl, Lsurf, label=r'$L_{surf}$', linestyle=':')
            # up.add_curve(ax, wl, Lsky.mean(axis=0) * rho, Lsky.std(axis=0) * rho,
            #              label=method + ' (' + str(round(rho, 4)) + ')', c='violet')
            ax.set_ylabel(r'$L_t\ or\ L_{surf}\ (mW\ m^{-2}\ nm^{-1}\ sr^{-1})$')
            ax.set_xlabel(r'Wavelength (nm)')
            ax.legend(loc='best', frameon=False)

            # ---- Proportion o(Lt - Lsurf ) /Lt
            ax = axs[1, 0]
            add_envelope(ax, wl, Lsurf / Ltm, label=r'$L_{surf}/L_t$')
            # up.add_curve(ax, wl, Lsky.mean(axis=0) * rho / Ltm, Lsky.std(axis=0) * rho,
            #              label=method + ' (' + str(round(rho, 4)) + ')', c='violet')
            ax.set_ylabel(r'$L_{surf}/L_t$')
            ax.set_xlabel(r'Wavelength (nm)')

            # ---- Lw
            ax = axs[1, 1]
            add_envelope(ax, wl, Rrs * Edm)
            # up.add_curve(ax, wl, Rrs.mean(axis=0) * Edm, Rrs.std(axis=0) * Edm,
            #              label=method + ' (' + str(round(rho, 4)) + ')', c='blue')

            ax.set_ylabel(r'$L_{w}\ (mW\ m^{-2}\ nm^{-1}\ sr^{-1})$')
            ax.set_xlabel(r'Wavelength (nm)')

            # ---- Rrs
            ax = axs[1, 2]
            add_envelope(ax, wl, Rrs)
            # up.add_curve(ax, wl, Rrs.transpose().mean(axis=1), Rrs.transpose().std(axis=1),
            #              label=method + ' (' + str(round(rho, 4)) + ')', c='blue')
            ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
            ax.set_xlabel(r'Wavelength (nm)')
            ax.set_title('azi=' + str(azi) + ', vza=' + str(vza) + ', sza=' + str(round(sza.mean(), 2)))

            fig.suptitle('trios_awr ' + self.name + ' idpr' + self.idpr, fontsize=18)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            fig.savefig(plot_file)
            plt.close()

        return self.Rrs

    def process_wrapper(self, wl, df, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):
        '''
        Wrapper to call standard processing upon pandas multiindex dataframe format
        :param wl:
        :param df:
        :param sza:
        :param vza:
        :param azi:
        :param ws:
        :param aot:
        :param method:
        :return:
        '''

        print(sza, vza, azi, ws, aot, method)
        if method == 'temp_opt':
            Rrs, rho, Rrs_opt, Rrs_opt_std = \
                self.process_optimization(wl, df.Lt.values, df.Lsky.values, df.Ed.values, sza, vza=vza, azi=azi)
        else:
            Rrs, rho = self.process(wl, df.Lt, df.Lsky.values, df.Ed.values, sza, vza, azi, ws, aot, method)

        if 'osoaa' in method:
            rho = pd.DataFrame(data=rho)
            rho.columns = pd.MultiIndex.from_product([['rho'], wl], names=['param', 'wl'])

        elif method == 'temp_opt':
            rho = pd.DataFrame(data=rho)
            rho.columns = pd.MultiIndex.from_arrays([['rho'], [550]], names=['param', 'wl'])

        else:
            rho = pd.DataFrame(data=[rho]*Rrs.shape[0])

            rho.columns = pd.MultiIndex.from_arrays([['rho'], [550]], names=['param', 'wl'])
        rho.index = df.Lt.index
        print(rho)
        Rrs.columns = pd.MultiIndex.from_product([['Rrs'], wl], names=['param', 'wl'])
        Rrs.index = df.Lt.index
        return Rrs, rho

    def process(self, wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):
        '''
        Standard processing based on estimation of Lsurf based on rho-factor and Lsky
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
            rho = self.get_rho_values(sza, vza, azi, wl=wl, ws=ws, aot=aot)
        elif method == 'M99':
            rho = self.get_rho_mobley( [np.median(sza)], [vza], [azi], [ws],method=method)
        elif method == 'M15':
            rho = self.get_rho_mobley( [np.median(sza)], [vza], [azi], [ws],method=method)
        else:
            return print('ERROR: no method for rho factor')
        print('rhoooooooooooo',rho.shape,Lsky.shape)
        self.Rrs = (Lt - rho * Lsky) / Ed
        # Rrs.columns = pd.MultiIndex.from_product([['Rrs'], wl], names=['param', 'wl'])

        return self.Rrs, rho #.mean()

    def cost_func(self, x, param, meas, Rrs_bar):
        sza, vza, azi = param
        Lt, Lsky, Ed = meas
        ws = x

        rho = self.get_rho_mobley([sza], [vza], [azi], [ws],method='M99')

        Rrs = (Lt - rho * Lsky) / Ed

        return Rrs - Rrs_bar

        # x_ave = x.mean()
        # return np.sum(np.abs(x - x_ave))

    def process_optimization(self, wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):

        # ------------------------------
        # initialization of mean/median values
        # ------------------------------
        Rrs, rho = self.process(wl, Lt, Lsky, Ed, sza, ws=ws, azi=azi)
        Rrs_bar = np.nanmedian(Rrs, axis=0)

        # -----------------------------
        # non-linear optimization
        # -----------------------------

        for j in range(1):
            x_est = []
            res = []
            Rrs_est = []
            rho_est = []

            for i in range(len(Lt)):
                geom = [sza[i], vza, azi]
                meas = [Lt[i], Lsky[i], Ed[i]]
                x0 = ws
                res_lsq = so.least_squares(self.cost_func, x0, bounds=(0, 15), args=(geom, meas, Rrs_bar))
                res.append(res_lsq)
                x_est.append(res_lsq.x[0])
                Rrs, rho = self.process(wl, Lt[i], Lsky[i], Ed[i], sza[i], ws=res_lsq.x[0], azi=azi)
                Rrs_est.append(Rrs)
                rho_est.append(rho[0])
                print(res_lsq.x, res_lsq.cost)

            Rrs_bar = np.mean(Rrs_est, axis=0)
            Rrs_std = np.std(Rrs_est, axis=0)
            Rrs_est = pd.DataFrame(Rrs_est)
            Rrs_est.columns = pd.MultiIndex.from_product([['Rrs'], wl], names=['param', 'wl'])

        return Rrs_est, rho_est, Rrs_bar, Rrs_std  # np.mean(rho_est, axis=0), Rrs_bar, Rrs_std

    @classmethod
    def filtering(cls, Lt, Lsky, Ed, **kargs):
        '''

        :param Lt:
        :param Lsky:
        :param Ed:
        :param kargs:
        :return:
        '''

        ind_Ed, notused = filters.spectra_median_filter(Ed, *kargs)
        ind_sky, notused = filters.spectra_median_filter(Lsky, *kargs)

        ind_surf, notused = filters.spectra_minimum_filter(Lt)
        print(ind_Ed.shape, ind_sky.shape, ind_surf.shape)
        ind = ind_Ed & ind_sky & ind_surf
        return ind


class swr_process:
    def __init__(self, df=None, wl=None, ):
        self.df = df
        self.wl = wl

    def call_process(self, ofile="", shade_corr=False):
        '''

        :param ofile: if ofile is given, Rrs results are written in ofile
        :param shade_corr:
        :return:
        '''
        wl = self.wl
        Lu = self.df.loc[:, ("Lu0+")]
        Ed = self.df.loc[:, ("Ed")]
        sza = self.df.loc[:, ("sza")].values.mean()
        Rrs = self.process(Lu, Ed, sza, wl, shade_corr=shade_corr)
        Rrs.columns = pd.MultiIndex.from_product([['Rrs(swr)'], Rrs.columns], names=['param', 'wl'])
        self.Rrs = Rrs

        if ofile:
            Rrs_df = pd.concat([self.df, Rrs], axis=1)
            Rrs_df.to_csv(ofile)
            Rrs_stat = Rrs.describe()
            Rrs_stat.columns = Rrs_stat.columns.droplevel()
            Rrs_stat = Rrs_stat.T
            # Rrs_stat.to_csv(ofile)
        return Rrs

    def process(self, Lu, Ed, sza, wl, R=0.05, shade_corr=False):
        Rrs = Lu / Ed
        ang_w = calc().angle_w(sza)

        iopw = ua.iopw()
        iopw.load_iopw()
        iopw.get_iopw(wl)
        a, bb = iopw.aw, iopw.bbw
        # TODO add particulate and dissolved component to a and bb values
        # a,bb = aux.get_iop(..., withwater=True)
        acdom = ua.cdom(0.5, wl).get_acdom()
        a = a + acdom
        bb = bb
        if shade_corr:
            Rrs = self.shade_corr(Rrs, R, ang_w, a, bb, wl)
        # Rrs.columns = pd.MultiIndex.from_product([['Rrs(swr)'], Rrs.columns], names=['param', 'wl'])
        self.Rrs = Rrs
        self.a = a
        self.bb = bb
        self.acdom = acdom
        return self.Rrs

    def epsilon(self, K, R, ang_w):
        '''
        epsilon from Shang et al, 2017, Applied Optics
        :param K:
        :param R:
        :param ang_w: Sun zenith angle below surface (in deg)
        :return:
        '''

        self.eps = np.array(1 - np.exp(-K * R / np.tan(np.radians(ang_w))))
        return self.eps

    def K(self, a, bb, ang_w):
        '''
        K (sum attenuation coef. of Lu in and outside the shade) from Shang et al, 2017, Applied Optics
        :param a: total absorption coefficient (m-1)
        :param bb: total backscattering coefficient (m-1)
        :param ang_w: Sun zenith angle below surface (in deg)
        :return:
        '''
        sin_ang_w = np.sin(np.radians(ang_w))
        self.K_ = (3.15 * sin_ang_w + 1.15) * a * np.exp(-1.57 * bb) \
                  + (5.62 * sin_ang_w - 0.23) * bb * np.exp(-0.5 * a)
        return self.K_

    def shade_corr(self, Rrs, R, ang_w, a, bb, wl, wl_cutoff=900):
        '''
        Correction of shading error from Shang et al, 2017, Applied Optics
        :param Rrs:
        :param R:
        :param ang_w:
        :param a:
        :param bb:
        :return:
        '''

        K = self.K(a, bb, ang_w)
        eps = self.epsilon(K, R, ang_w)
        eps[wl > wl_cutoff] = 0
        self.Rrs = Rrs / (1 - eps)
        return self.Rrs


class iwr_process:
    def __init__(self, df=None, wl=None, name="", idpr=""):
        self.df = df
        self.wl = wl
        self.name = name
        self.idpr = idpr

        ################
        # load aux data
        iopw = ua.iopw()
        iopw.load_iopw()
        irr = ua.irradiance()
        irr.load_F0()
        self.aw, self.bbw = iopw.get_iopw(wl)
        self.F0 = irr.get_F0(wl)


    def call_process(self, ofile="", plot_file=""):

        # ---------------------------
        # Data formatting
        # ---------------------------
        df = self.df
        wl = self.wl

        # set pandas dataframes with data and parameters
        df['rounded_depth', ''] = df.prof_Edz.round(1)

        # data filtering
        df[df.rounded_depth < -0.05] = np.nan
        # df[df == 0] = np.nan

        mean = df.groupby('rounded_depth').mean()
        median = df.groupby('rounded_depth').median()
        std = df.groupby('rounded_depth').std()
        df_ = df.drop(df.columns[df.dtypes == 'object'], axis=1)
        q25 = df_.groupby('rounded_depth').quantile(0.25)
        q75 = df_.groupby('rounded_depth').quantile(0.75)

        # ---------------------------
        # Data processing
        # ---------------------------
        res, res_Lu = self.process(mean, std)
        res_w, void = self.process(mean, std, mode='log')
        res_med, void = self.process(median, std)
        res_lsq, void = self.process(mean, std, mode='lsq')

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
        Rrs_df['Kd'] = Kd

        if ofile:
            # -----------------------------------------------
            #   write result data
            # -----------------------------------------------
            Rrs_df.to_csv(ofile)

        if plot_file:
            # ---------------------------
            # Plotting section
            # ---------------------------

            pdf = PdfPages(plot_file)

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

            print(std.prof_Edz)

            for idx in idx_list_for_plot:
                ax = axs.flat[i]
                i += 1
                y = mean.Edz.iloc[:, idx]
                sigma = std.Edz.iloc[:, idx]
                yerr = [median.Edz.iloc[:, idx] - q25.Edz.iloc[:, idx], q75.Edz.iloc[:, idx] - median.Edz.iloc[:, idx]]

                xerr = std.prof_Edz.values  # [mean.prof_Edz - std.prof_Edz, mean.prof_Edz + std.prof_Edz]
                yerr = [mean.Edz.iloc[:, idx] - std.Edz.iloc[:, idx], mean.Edz.iloc[:, idx] + std.Edz.iloc[:, idx]]

                ax.errorbar(x, y, xerr=xerr, yerr=yerr, ecolor='gray', alpha=0.73, fmt='none', elinewidth=2, capsize=0,
                            label='std')
                ax.scatter(x, y,
                           # s = std.Edz.iloc[:,50]/std.Edz.iloc[:,50].mean()+100,
                           c=mean.Ed.iloc[:, idx],
                           alpha=0.6, cmap=cmocean.cm.thermal, label=None
                           )
                Ed_sim = self.f_Edz(depth_, *res.popt[idx, :])
                ax.plot(depth_, Ed_sim, linestyle='-', c='black', label='mean')
                Ed_sim = self.f_Edz(depth_, *res_w.popt[idx, :])
                ax.plot(depth_, Ed_sim, linestyle=':', c='black', label='log-space')
                # Ed_sim = self.f_Edz(depth_, *res_rw[idx][0])
                # ax.plot(depth_, Ed_sim, linestyle=':', c='red', label='mean, relativ weighted')
                Ed_sim = self.f_Edz(depth_, *res_med.popt[idx, :])
                ax.plot(depth_, Ed_sim, linestyle='--', c='black', label='median')

                ax.semilogy()
                # ax.semilogx()
                # ax.colorbar()
                ax.set_ylabel(r'$E_d\ ({mW\cdot m^{-2}\cdot nm^{-1}})$')
                ax.set_xlabel('Depth (m)')
                ax.set_title(r'$\lambda = $' + str(round(wl[idx], 1)) + ' nm, Kd = ' +
                             str(round(res_w.popt[idx, 0], 3)) + '$m^{-1}$, Ed0 =' + str(round(res_w.popt[idx, 1], 1)),
                             fontsize=12)
            ax.legend(loc='best', frameon=False)
            fig.suptitle('trios_iwr ' + self.name + ' idpr' + self.idpr, fontsize=16)
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
                ax.errorbar(x, y, yerr=yerr, ecolor='lightgray', alpha=0.43, fmt='none', elinewidth=3, capsize=0,
                            label='std')
                ax.scatter(x, y,
                           # s = std.Luz.iloc[:,50]/std.Luz.iloc[:,50].mean()+100,
                           c=mean.Ed.iloc[:, idx],
                           alpha=0.6, cmap=cmocean.cm.thermal, label=None
                           )
                Lu_sim = self.f_Lu(depth_, *res_lsq.popt[idx, -3:-1])
                ax.plot(depth_, Lu_sim, linestyle='-', c='black', label='mean')
                # Lu_sim = self.f_Luz(depth_, *res_w.popt[idx, :])
                # ax.plot(depth_, Lu_sim, linestyle=':', c='black', label='log-space')
                # # Lu_sim = self.f_Luz(depth_, *res_rw[idx][0])
                # # ax.plot(depth_, Lu_sim, linestyle=':', c='red', label='mean, relativ weighted')
                # Lu_sim = self.f_Luz(depth_, *res_med.popt[idx, :])
                # ax.plot(depth_, Lu_sim, linestyle='--', c='black', label='median')

                ax.semilogy()
                # ax.semilogx()

                # ax.colorbar()
                ax.set_ylabel(r'$L_u\ ({mW\cdot m^{-2}\cdot sr^{-1}\cdot nm^{-1}})$')
                ax.set_xlabel('Depth (m)')
                ax.set_title(r'$\lambda = $' + str(round(wl[idx], 1)) + r'$ nm, K_{Lu} = $' +
                             str(round(res_lsq.popt[idx, 2], 3)) + '$m^{-1}$, Lu0 =' + str(
                    round(res_lsq.popt[idx, 3], 1)),
                             fontsize=12)
            ax.legend(loc='best', frameon=False)
            fig.suptitle('trios_iwr ' + self.name + ' idpr' + self.idpr, fontsize=16)
            # fig.savefig(os.path.join(dirfig,'trios_iw_idpr' + idpr + '_'+name+'.png'),dpi=600)
            # fig.savefig(os.path.join(dirfig, 'trios_iwr_Lu_profile_idpr' + idpr + '_' + name + '.pdf'))
            pdf.savefig()

            plt.close()

            # -------------------------------------
            # Spectra plotting
            # -------------------------------------
            mpl.rcParams.update({'font.size': 18})
            fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15, 12))
            fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.29)

            iparam = 0
            ax = axs.flat[iparam]
            y, sigma = res.popt[:, iparam], res.pcov[:, iparam, iparam]
            up.add_curve(ax, wl, y, sigma, c='red', label='mean')
            y, sigma = res_w.popt[:, iparam], res_w.pcov[:, iparam, iparam]
            up.add_curve(ax, wl, y, sigma, c='orange', label='log-space')
            y, sigma = res_med.popt[:, iparam], res_med.pcov[:, iparam, iparam]
            up.add_curve(ax, wl, y, sigma, c='green', label='median')
            y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam])
            up.add_curve(ax, wl, y, sigma, c='blue', label='lsq')
            ax.set_ylabel(r'$K_{d}\ ({m^{-1}})$')
            ax.set_xlabel(r'Wavelength (nm)')

            iparam = 1
            ax = axs.flat[iparam]
            y, sigma = res.popt[:, iparam], np.sqrt(res.pcov[:, iparam, iparam])
            up.add_curve(ax, wl, y, sigma, c='red', label='mean')
            y, sigma = res_w.popt[:, iparam], np.sqrt(res_w.pcov[:, iparam, iparam])
            up.add_curve(ax, wl, y, sigma, c='orange', label='log-space')
            y, sigma = res_med.popt[:, iparam], np.sqrt(res_med.pcov[:, iparam, iparam])
            up.add_curve(ax, wl, y, sigma, c='green', label='median')
            y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam])
            up.add_curve(ax, wl, y, sigma, c='blue', label='lsq')
            up.add_curve(ax, wl, mean.Ed.mean(), mean.Ed.std(), c='black', label='Es')
            ax.legend(loc='best', frameon=False)
            ax.set_ylabel(r'$E_{d}(0^{-})\ ({mW\cdot m^{-2}\cdot nm^{-1}})$')
            ax.set_xlabel(r'Wavelength (nm)')

            iparam = 2
            ax = axs.flat[iparam]
            y, sigma = res_Lu.popt[:, 0], np.sqrt(res_Lu.pcov[:, 0,0])
            up.add_curve(ax, wl, y, sigma, c='red', label='mean')
            y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam]) * 1e-3
            sigma[sigma > 10 * np.nanmedian(sigma)] = np.nan
            up.add_curve(ax, wl, y, sigma, c='blue', label='lsq')
            ax.set_ylabel(r'$K_{Lu}\ ({m^{-1}})$')
            ax.set_xlabel(r'Wavelength (nm)')

            iparam = 3
            ax = axs.flat[iparam]

            y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam]) * 1e-3
            sigma[sigma > 10 * np.nanmedian(sigma)] = np.nan
            up.add_curve(ax, wl, y, sigma, c='blue', label='lsq')
            y, sigma = res_Lu.popt[:, 1], np.sqrt(res_Lu.pcov[:, 1, 1])
            up.add_curve(ax, wl, y, sigma, c='red', label='mean')
            ax.legend(loc='best', frameon=False)
            ax.set_ylabel(r'$L_{w}(0^{-})\ (mW\cdot m^{-2}\cdot sr^{-1}\cdot nm^{-1})$')
            ax.set_xlabel(r'Wavelength (nm)')

            fig.suptitle('trios_iwr ' + self.name + ' idpr' + self.idpr, fontsize=16)
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
            up.add_curve(ax, wl, Es / Edm, c='red', label='Es over Ed0-')
            ax.set_xlabel(r'Wavelength (nm)')
            ax.legend(loc='best', frameon=False)
            ax = axs.flat[1]
            up.add_curve(ax, wl, rrs, c='red', label='rrs')
            ax.set_xlabel(r'Wavelength (nm)')
            ax.legend(loc='best', frameon=False)
            ax = axs.flat[2]
            up.add_curve(ax, wl, Rrs_Edm, c='red', label='Rrs (from Ed0-)')
            up.add_curve(ax, wl, Rrs_Edp, c='blue', label='Rrs (from Ed0+)')

            ax.set_xlabel(r'Wavelength (nm)')
            ax.legend(loc='best', frameon=False)

            ax = axs.flat[3]
            iparam = 0
            y, sigma = res.popt[:, iparam], res.pcov[:, iparam, iparam]
            ax.plot(wl, np.gradient(np.gradient(y, wl), wl), c='red', label='mean')
            y, sigma = res_w.popt[:, iparam], res_w.pcov[:, iparam, iparam]
            ax.plot(wl, np.gradient(np.gradient(y, wl), wl),  c='orange', label='log-space')
            y, sigma = res_med.popt[:, iparam], res_med.pcov[:, iparam, iparam]
            ax.plot(wl, np.gradient(np.gradient(y, wl), wl),  c='green', label='median')
            y, sigma = res_lsq.popt[:, iparam], np.sqrt(res_lsq.pcov[:, iparam, iparam])
            ax.plot(wl, np.gradient(np.gradient(y, wl), wl),  c='blue', label='lsq')
            ax.set_xlabel(r'Wavelength (nm)')
            ax.set_ylabel(r'$\delta K_{d}/\delta\lambda\ ({m^{-1}})$')


            fig.suptitle('trios_iwr ' + self.name + ' idpr' + self.idpr, fontsize=16)
            pdf.savefig()
            plt.close()
            pdf.close()

    def process(self, meas, std, mode='linear'):
        wl_ = self.wl

        # TODO check noise values (e.g., NEI from Trios), should it be spectral?
        noise = 0.1

        N = len(wl_)
        x = meas.prof_Edz  # - 0.56
        res = uu.fit(N)
        res_Lu = uu.fit(N)
        if mode == 'lsq':
            res = uu.fit(N, 5)
        first =True
        for idx, wl in enumerate(wl_[:-10]):
            aw, bbw = self.aw[idx], self.bbw[idx]
            F0 = self.F0[idx]

            y = meas.Edz.iloc[:, idx]
            sigma = std.Edz.iloc[:, idx]
            sigma[sigma < noise] = noise
            sigma.fillna(np.inf, inplace=True)

            if mode == 'linear':
                res.popt[idx, :], res.pcov[idx, ...] = so.curve_fit(self.f_Edz, x, y, [1.1 * aw, 0.9*F0],
                                                                    bounds=([aw, 0], [np.inf, F0]))
                xx = meas.prof_Luz
                yy = meas.Luz.iloc[:, idx]
                res_Lu.popt[idx, :], res_Lu.pcov[idx, ...] = so.curve_fit(self.f_Lu, xx, yy, [1.1 * aw, 0.005],
                                                                bounds=([aw/2, 0], [np.inf, np.inf]))
            elif mode == 'log':
                res.popt[idx, :], res.pcov[idx, ...] = so.curve_fit(self.f_logEdz, x, np.log(1 + y),
                                                                    [1.1 * aw, 100], bounds=(
                        [aw, 0], [np.inf, F0]))  # , sigma=sigma, absolute_sigma=True

            elif mode == 'lsq':

                z = (meas.prof_Edz, meas.prof_Luz)
                y = (meas.Edz.iloc[:, idx], meas.Luz.iloc[:, idx], meas.Ed.iloc[:, idx].mean())

                sig_Edz = self.format_sigma(std.Edz.iloc[:, idx], meas.Edz.iloc[:, idx], 0.1)
                sig_Luz = self.format_sigma(std.Luz.iloc[:, idx], meas.Luz.iloc[:, idx], 1e-3)

                sigma = (sig_Edz, sig_Luz)
                sigma = (meas.Ed.iloc[:, idx].values,meas.Ed.iloc[:, idx].values)
                #sigma = (1, 1)

                x0 = [1.1 * aw, min(F0, meas.Ed.iloc[:, idx].mean()), 1.1 * aw, meas.Luz.iloc[0, idx], 0]
                print('x0', wl, F0, aw, x0)
                print(([aw, 0, aw / 2, 0,-0.1], [np.inf, F0, np.inf, np.inf,0.1]))
                lsq = so.least_squares(self.cost_func, x0, args=(z, y, sigma),
                                       bounds=([0, 0, 0, 0,-0.1], [np.inf, F0, np.inf, np.inf,0.1]))
                print('x', wl, F0, aw, lsq.x)
                cost = 2 * lsq.cost  # res.cost is half sum of squares!
                res.popt[idx, :], res.pcov[idx, ...] = lsq.x, calc().cov_from_jac(lsq.jac, cost)

                # x0 = [ 1.1 * aw, meas.Luz.iloc[0, idx]]
                # llsq = so.least_squares(self.f_Lu, x0, args=(z, y),
                #                        bounds=([ aw / 2, 0], [ np.inf, np.inf]))
                # print(llsq.x,llsq.x[1]-lsq.x[1])
            # TODO add lmfit methods and constraints on Z and Es ~ 0.96 * Ed0-

            #if mode == 'lsq':
                # TODO formalize and do more clever things for Quality Control
                # discard retrieval if error covariance > threshold error covariance median
                QC_idx = res.pcov[:, 3, 3] > 20 * np.nanmedian(res.pcov[:, 3, 3])

                res.popt[QC_idx, 3] = np.nan

        return res, res_Lu

    def format_sigma(self, sigma, rescale=1, noise=0.1):
        '''

        :param sigma:
        :return:
        '''

        sigma = (sigma + noise) / rescale

        return sigma.fillna(np.inf)

    # @staticmethod
    def f_Edz(self, depth, Kd, Ed0):
        '''simple Edz model for homogeneous water column'''
        return Ed0 * np.exp(-Kd * depth)

    # @staticmethod
    def f_logEdz(self, depth, Kd, Ed0):
        '''simple Edz model for homogeneous water column'''
        return np.log(1 + self.f_Edz(depth, Kd, Ed0))  # Ed0) -Kd*depth

    def f_Lu(self, depth, KLu, Lw0minus):
        '''simple Edz model for homogeneous water column'''
        return Lw0minus * np.exp(-KLu * depth)

    def cost_func(self, x, z, mes, sigma):

        z_Edz = z[0]
        z_Lu = z[1]
        Edz = mes[0]
        Lu = mes[1]
        Es = mes[2]
        sig_Edz = sigma[0]
        sig_Luz = sigma[1]
        sig_cons = Es #sigma[2]

        cost_f1 = (Edz - self.f_Edz(z_Edz+x[4], x[0], x[1])) / sig_Edz
        cost_f2 = (Lu - self.f_Lu(z_Lu+x[4], x[2], x[3])) / sig_Luz

        constraint = [(Es - 0.96 * self.f_Edz(0+x[4], x[0], x[1]))/ sig_cons,
                      x[4]]

        return np.concatenate([cost_f1, cost_f2, constraint])

    def cost_func_multi_wl(self, x, z, mes, sigma):
        # TODO write and try out multi-wavelength approach with smoothness constraints (e.g., second derivatives)
        z_Edz = z[0]
        z_Lu = z[1]

        Edz = mes[0]
        Lu = mes[1]
        Es = mes[2]

        sig_Edz = sigma[0]
        sig_Luz = sigma[1]
        sig_cons = Es #sigma[2]

        cost_f1 = (Edz - self.f_Edz(z_Edz+x[4], x[0], x[1])) / sig_Edz
        cost_f2 = (Lu - self.f_Lu(z_Lu+x[4], x[2], x[3])) / sig_Luz

        constraint = [(Es - 0.96 * self.f_Edz(0+x[4], x[0], x[1]))/ sig_cons,
                      x[4]]

        return np.concatenate([cost_f1, cost_f2, constraint])
    def Kd(self, depth, Edz):
        Kd = np.diff(Edz) / np.diff(depth)

    def plot_raw(self, x='Luz', y='prof_Luz'):
        trace = go.Scattergl(
            x=self.df[x].values,
            y=self.df[y].values,

            text=self.df.index.get_level_values(0),
            hoverinfo="text",
            marker={
                'size': 7,
                'opacity': 0.5,
                # 'color': 'rgba({}, {}, {}, {})'.format(*s_m.to_rgba(parameters[i]).flatten()),
                # x.unique(),#color': df.index.get_level_values(0),
                'line': {'width': 0.5, 'color': 'white'},
            },
            # error_y=ErrorY(
            #     type='data',
            #     array=df['Error'],
            #     thickness=1.5,
            #     width=2,
            #     color='#B4E8FC'
            # ),

        )

        layout = go.Layout(
            height=450,
            xaxis=dict(
                range=[0, 200],
                showgrid=False,
                showline=False,
                zeroline=False,
                fixedrange=True,
                tickvals=[0, 50, 100, 150, 200],
                ticktext=['200', '150', '100', '50', '0'],
                title=''
            ),
            yaxis=dict(
                range=[min(-5, min(self.df[y])),
                       max(0, max(self.df[y]))],
                showline=False,
                fixedrange=True,
                zeroline=False,
                # nticks=max(6, round(df['Speed'].iloc[-1]/10))
            ),
            margin=go.Margin(
                t=45,
                l=50,
                r=50
            )
        )

        return go.Figure(data=[trace], layout=layout)


class self_shading:
    def __init__(self):
        '''GordonDing 1992 values for epsilon'''

        self.ang = np.linspace(0, 90, 10)
        self.eps_dir_LuZ = [2.17, 2.17, 2.23, 2.23, 2.29, 2.37, 2.41, 2.45, 2.45, 2.45]
        self.eps_dir_EuZ = [3.14, 3.14, 3.05, 2.94, 2.80, 2.64, 2.47, 2.33, 2.33, 2.33]
        self.eps_dif_LuZ = 4.61,
        self.eps_dif_EuZ = 2.70

    def epsilon(self, sza):
        eps = interpolate.interp1d(self.ang, self.eps_dif_EuZ)(sza)
        return eps


class filters:
    def __init__(self):
        pass

    @classmethod
    def spectra_median_filter(cls, spectra, threshold=0.1):
        '''
        Get indices of the spectra whose values are close to median in the acquisition sequence
        :param spectra: 2D pandas data frame
        :param threshold: defines the interval of retained spectra (median +/- threhold * median)
        :return: boolean indices, array of data within interval median +/- threshold
        '''

        spec = spectra.sum(axis=1)
        med = spec.median()
        ind = np.abs(1 - spec / med) < threshold
        return ind, spectra[ind]

    @classmethod
    def spectra_minimum_filter(cls, spectra, idxs=None, cutoff=0.2, threshold=0.1):
        '''
        Get indices of the spectra whose values are minimum in the acquisition sequence
        :param spectra: 2D pandas data frame
        :param idxs: indices on which the filter is applied (e.g., select only IR wavelengths)
        :param cutoff: percentage of minimum spectra used to get the median value
        :param threshold: defines the interval of retained spectra (median +/- threhold * median)
        :return: boolean indices, array of data within interval minimum to minimum * (1 + threshold)
        '''
        print(spectra.shape)
        if idxs == None:
            idxs = np.arange(spectra.shape[1])

        spec = np.nansum(spectra.values[:, idxs], axis=1)
        N = spec.shape[0]
        idx_num = int(cutoff * N)
        idx_min = spec.argsort()[:idx_num]

        med = np.nanmedian(spec[idx_min])
        ind = np.abs(1 - spec / med) < threshold

        return ind, spectra[ind]

# TODO add time rolling filters


class calc:
    def __init__(self):
        pass

    def PAR(self, wl, Ed):
        '''
        Compute instantaneous PAR from Ed spectrum.
        PAR in mW m-2
        PAR_quanta in µmol photon m-2 s-1
        :param wl:
        :param Ed:
        :return:
        '''
        # ---------------------------------------------
        #      PARAMETERS
        # Planck constant in J s or W s2
        h = 6.6260695729e-3  # d-34
        # light speed in m s-1
        c = 2.99792458e0  # d8
        # Avogadro Number in mol-1
        Avogadro = 6.0221412927e0  # d23
        hc = Avogadro * h * c
        # ---------------------------------------------

        idx_par = (wl >= 400) & (wl <= 700)
        wl = wl[idx_par]
        Ed = Ed[idx_par]
        par = integrate.trapz(Ed, wl)
        par_quanta = integrate.trapz(np.multiply(wl, Ed), wl) / hc
        return par, par_quanta

    def earth_sun_correction(self, dayofyear):
        '''
        Earth-Sun distance correction factor for adjustment of mean solar irradiance

        :param dayofyear:
        :return: correction factor
        '''
        theta = 2. * np.pi * dayofyear / 365
        d2 = 1.00011 + 0.034221 * np.cos(theta) + 0.00128 * np.sin(theta) + \
             0.000719 * np.cos(2 * theta) + 0.000077 * np.sin(2 * theta)
        return d2

    def bidir(self, sza, vza, azi):

        bidir = 1

        return bidir

    def angle_w(self, angle_air, n=1.334):
        '''
        convert above surface angle (angle_air) into sub-surface angle
        :param angle_air in deg
        :param n: refractive index of water
        :return: sub-surface angle in deg
        '''
        return np.degrees(np.arcsin(np.sin(np.radians(angle_air)) / n))

    def spline_2d(self, gin, arr, gout):
        '''
        Interpolation of a 6D array (arr) with bicubic splines on a 2D grid
        corresponding to the 5th and 6th dimensions of arr.
        Return 4D array interpolated on gout.

        :param gin: regular 2D grid of the tabulated data (tuple/array/list of arrays)
        :param arr: tabulated data (N dimensions, interpolation on N-1 and N)
        :param gout: new 2D grid on which data are interpolated (with dims 2 and 3 of the same length);
                    (tuple/array/list of arrays)
        :return: Interpolated data (1D or 3D array depending on the dimension shapes of gout
        '''

        N = arr.shape
        interp = np.zeros(N[:-2])

        for i in range(N[0]):
            for j in range(N[1]):
                for k in range(N[2]):
                    for l in range(N[3]):
                        interp[i, j, k, l] = interpolate.RectBivariateSpline(gin[0], gin[1], arr[i, j, k, l, ...])(
                            gout[0], gout[1], grid=False)

        return interp

    def spline_4d(self, gin, lut, gout):
        '''
        Interpolation with two successive bicubic splines on a regular 4D grid.
        Designed for interpolation in radiative transfer look-up tables with the two last dimensions
        (i.e., wavelength and solar zenith angle) of the same length.
        Those dimensions are then reduced/merged to a single one to get interpolated data on a 3D grid.

        :param gin: regular 4D grid of the tabulated data (tuple/array/list of arrays)
        :param lut: tabulated data
        :param gout: new 4D grid on which data are interpolated (with dims 2 and 3 of the same length);
                    (tuple/array/list of arrays)
        :return: Interpolated data (1D or 3D array depending on the dimension shapes of gout
        '''

        N = gin[0].__len__(), gin[1].__len__(), gin[2].__len__(), gin[3].__len__()
        Nout = gout[0].__len__(), gout[1].__len__(), gout[2].__len__()
        tmp = np.zeros([N[0], N[1], Nout[2]])

        for i in range(N[0]):
            for j in range(N[1]):
                tmp[i, j, :] = interpolate.RectBivariateSpline(gin[2], gin[3], lut[i, j, :, :])(gout[2], gout[3],
                                                                                                grid=False)
        if Nout[0] == Nout[1] == 1:
            interp = np.ndarray(Nout[2])
            for iband in range(Nout[2]):
                interp[iband] = interpolate.RectBivariateSpline(gin[0], gin[1], tmp[:, :, iband])(gout[0], gout[1],
                                                                                                  grid=False)
        else:
            interp = np.ndarray([Nout[0], Nout[1], Nout[2]])
            for iband in range(Nout[2]):
                interp[:, :, iband] = interpolate.RectBivariateSpline(gin[0], gin[1], tmp[:, :, iband])(gout[0],
                                                                                                        gout[1],
                                                                                                        grid=True)

        return interp

    def cov_from_jac(self, jac, cost):
        '''
        Compute covariance from jacobian matrix computed in optimization processes
        Use Moore-Penrose inverse discarding zero singular values.
        :param jac: jacobian
        :param cost: cost residual
        :return: pcov
        '''

        from scipy.linalg import svd

        M, N = jac.shape
        _, s, VT = svd(jac, full_matrices=False)
        threshold = np.finfo(float).eps * max(jac.shape) * s[0]
        s = s[s > threshold]
        VT = VT[:s.size]
        pcov = np.dot(VT.T / s ** 2, VT)

        if pcov is None:
            # indeterminate covariance
            pcov = np.zeros((N, N), dtype=float)
            pcov.fill(np.inf)
        else:
            s_sq = cost / (M - N)
            pcov = pcov * s_sq
        return pcov
