import pandas as pd
from scipy import interpolate, integrate
import scipy.optimize as so

import plotly.graph_objs as go

import trios.utils.utils as uu
from trios.utils.utils import reshape as r
import trios.utils.auxdata as ua
from trios.config import *


class awr_process:
    def __init__(self, df=None, wl=None):
        self.df = df
        self.aot = 0.1
        self.ws = 2
        self.wl = wl
        self.rhosoaa_fine_file = rhosoaa_fine_file
        self.rhosoaa_coarse_file = rhosoaa_coarse_file
        self.M1999_file = M1999_file
        self.M2015_file = M2015_file
        self.load_rho_lut()
        self.rho = self.rhosoaa_fine

    def load_rho_lut(self):
        self.rhosoaa_fine = pd.read_csv(self.rhosoaa_fine_file, index_col=[0, 1, 2, 3, 4, 5])
        self.rhosoaa_coarse = pd.read_csv(self.rhosoaa_coarse_file, index_col=[0, 1, 2, 3, 4, 5])
        self.rhoM1999 = pd.read_csv(self.M1999_file, skiprows=7, index_col=[0, 1, 2, 3])
        self.rhoM2015 = pd.read_csv(self.M2015_file, skiprows=8, index_col=[0, 1, 2, 3])
        self.rhoM1999.dropna(inplace=True)
        self.rhoM2015.dropna(inplace=True)

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

        grid = self.rho.rho.index.levels

        # convert pandas dataframe into 6D array of the tabulated rho values for interpolation
        rhoname = 'rho'
        if sunglint:
            rhoname = 'rho_g'

        rho_6d = r().df2ndarray(self.rho, rhoname)

        rho_ = calc().spline_2d(grid[-2:], rho_6d, (azi, vza))

        rho_wl = calc().spline_4d(grid[:-2], rho_, (ws, aot, wl, sza))

        return rho_wl.squeeze()

    def get_rho_mobley(self, rhodf, sza, vza, azi, ws):
        '''
        Get the Mobley rho factor from cubic interpolation in the tabulated values

        :param rhodf:
        :param sza:
        :param vza:
        :param azi:
        :param ws:
        :return:
        '''

        rhodf = rhodf.query('sza<75 & vza >0')
        rhodf.index = rhodf.index.remove_unused_levels()

        # grid {wind, sza, vza, azi}
        grid = rhodf.index.levels

        rho_ = r().df2ndarray(rhodf, 'rho')
        rho_mobley = calc().spline_4d(grid, rho_, (ws, sza, vza, azi))
        return rho_mobley

    # def call_process(self, vza=[40], azi=[135], ws=2, aot=0.1):
    #     wl = self.wl
    #     Lt = self.df.loc[:, ("Lt")]
    #     Lsky = self.df.loc[:, ("Lsky")]
    #     Ed = self.df.loc[:, ("Ed")]
    #     sza = self.df.loc[:, ("sza")].values.mean()
    #     Rrs = self.process(wl, Lt, Lsky, Ed, sza, vza, azi, ws, aot)
    #     Rrs.columns = pd.MultiIndex.from_product([['Rrs(awr)'], self.Rrs.columns], names=['param', 'wl'])
    #     self.Rrs = Rrs
    #     return self.Rrs
    #
    # def process(self, wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=[2], aot=[0.1]):
    #     '''
    #
    #     :param wl:
    #     :param Lt:
    #     :param Lsky:
    #     :param Ed:
    #     :param sza:
    #     :param vza:
    #     :param azi:
    #     :param ws:
    #     :param aot:
    #     :return:
    #     '''
    #
    #     # ------------------
    #     # filtering
    #     # ------------------
    #     ind_Ed, notused = calc.spectra_median_filter(Ed)
    #     ind_sky, notused = calc.spectra_median_filter(Lsky)
    #     ind = ind_Ed & ind_sky
    #     Lt, Lsky, Ed, sza = Lt[ind], Lsky[ind], Ed[ind], sza[ind]
    #
    #     rho = self.get_rho_values(np.median(sza), vza, azi, wl=wl, ws=ws, aot=aot)
    #     self.Rrs = (Lt - rho * Lsky.values) / Ed.values
    #     self.Rrs.columns = pd.MultiIndex.from_product([['Rrs(awr)'], wl], names=['param', 'wl'])
    #     return self.Rrs, rho

    @staticmethod
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

    def process_wrapper(self, wl, df, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):

        print(sza, vza, azi, ws, aot, method)
        Rrs, rho = self.process(wl, df.Lt, df.Lsky.values, df.Ed.values, sza, vza, azi, ws, aot, method)

        Rrs.columns = pd.MultiIndex.from_product([['Rrs'], wl], names=['param', 'wl'])

        return Rrs, rho

    def process(self, wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):
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
            rho = self.get_rho_values(np.median(sza), vza, azi, wl=wl, ws=ws, aot=aot)
        elif method == 'M99':
            rho = self.get_rho_mobley(self.rhoM1999, [np.median(sza)], [vza], [azi], [ws])
        elif method == 'M15':
            rho = self.get_rho_mobley(self.rhoM2015, [np.median(sza)], [vza], [azi], [ws])
        else:
            return print('ERROR: no method for rho factor')

        self.Rrs = (Lt - rho * Lsky) / Ed
        # Rrs.columns = pd.MultiIndex.from_product([['Rrs'], wl], names=['param', 'wl'])

        return self.Rrs, rho.mean()

    def cost_func(self, x, param, meas, Rrs_bar):
        sza, vza, azi = param
        Lt, Lsky, Ed = meas
        ws = x

        rho = self.get_rho_mobley(self.rhoM1999, [sza], [vza], [azi], [ws])

        Rrs = (Lt - rho * Lsky) / Ed

        return Rrs - Rrs_bar

        # x_ave = x.mean()
        # return np.sum(np.abs(x - x_ave))

    def process_optimization(self, wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=[2], aot=[0.1], method='M99'):

        # ------------------------------
        # initialization of mean/median values
        # ------------------------------
        Rrs, rho = self.process(wl, Lt, Lsky, Ed, sza, ws=ws, azi=azi)
        Rrs_bar = Rrs.mean(axis=0)

        # -----------------------------
        # non-linear optimization
        # -----------------------------

        for j in range(10):
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
                rho_est.append(rho)
                print(res_lsq.x, res_lsq.cost)
            Rrs_bar = np.mean(Rrs_est, axis=0)
            Rrs_std = np.std(Rrs_est, axis=0)
        return Rrs_bar, Rrs_std


class swr_process:
    def __init__(self, df=None, wl=None, ):
        self.df = df
        self.wl = wl

    def call_process(self, shade_corr=False):
        wl = self.wl
        Lu = self.df.loc[:, ("Lu0+")]
        Ed = self.df.loc[:, ("Ed")]
        sza = self.df.loc[:, ("sza")].values.mean()
        Rrs = self.process(Lu, Ed, sza, wl, shade_corr=shade_corr)
        Rrs.columns = pd.MultiIndex.from_product([['Rrs(swr)'], Rrs.columns], names=['param', 'wl'])
        self.Rrs = Rrs
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
    def __init__(self, df=None, wl=None, ):
        self.df = df
        self.wl = wl

    def process(self, meas, std, mode='linear'):
        wl_ = self.wl

        ################
        # load aux data
        iopw = ua.iopw()
        iopw.load_iopw()
        irr = ua.irradiance()
        irr.load_F0()
        # TODO check noise values (e.g., NEI from Trios), should it be spectral?
        noise = 0.1

        N = len(wl_)
        x = meas.prof_Edz  # - 0.56
        res = uu.fit(N)
        if mode == 'lsq':
            res = uu.fit(N, 4)

        for idx, wl in enumerate(wl_[:-10]):
            aw, bbw = iopw.get_iopw(wl)
            F0 = irr.get_F0(wl)

            y = meas.Edz.iloc[:, idx]
            sigma = std.Edz.iloc[:, idx]
            sigma[sigma < noise] = noise
            sigma.fillna(np.inf, inplace=True)
            if mode == 'linear':
                res.popt[idx, :], res.pcov[idx, ...] = so.curve_fit(self.f_Edz, x, y, [1.1 * aw, 100],
                                                                    bounds=([aw, 0], [np.inf, F0]))
            elif mode == 'log':
                res.popt[idx, :], res.pcov[idx, ...] = so.curve_fit(self.f_logEdz, x, np.log(1 + y),
                                                                    [1.1 * aw, 100], bounds=(
                        [aw, 0], [np.inf, F0]))  # , sigma=sigma, absolute_sigma=True
            elif mode == 'lsq':
                z = (meas.prof_Edz, meas.prof_Luz)
                y = (meas.Edz.iloc[:, idx], meas.Luz.iloc[:, idx])

                sig_Edz = self.format_sigma(std.Edz.iloc[:, idx],meas.Edz.iloc[:, idx], 0.1)
                sig_Luz = self.format_sigma(std.Luz.iloc[:, idx],meas.Luz.iloc[:, idx], 1e-3)

                sigma = (sig_Edz, sig_Luz)
                sigma = (1,1)
                x0 = [1.1 * aw, meas.Ed.iloc[:, idx].mean(), 1.1 * aw, meas.Luz.iloc[0, idx]]

                lsq = so.least_squares(self.cost_func, x0, args=(z, y, sigma),
                                       bounds=([aw, 0, aw/2, 0], [np.inf, F0, np.inf, np.inf]))
                cost = 2 * lsq.cost  # res.cost is half sum of squares!
                res.popt[idx, :], res.pcov[idx, ...] = lsq.x, calc().cov_from_jac(lsq.jac, cost)

            if mode == 'lsq':
                # TODO formalize and do more clever things for Quality Control
                # discard retrieval if error covariance > threshold error covariance median
                QC_idx = res.pcov[:, 3,3] > 20 *  np.nanmedian(res.pcov[:, 3,3])

                res.popt[QC_idx, 3] = np.nan

        return res

    def format_sigma(self, sigma, rescale = 1, noise=0.1):
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
        sig_Edz = sigma[0]
        sig_Luz = sigma[1]

        cost_f1 = (Edz - self.f_Edz(z_Edz, x[0], x[1])) / sig_Edz
        cost_f2 = (Lu - self.f_Lu(z_Lu, x[2], x[3])) / sig_Luz

        return np.append(cost_f1, cost_f2)

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

    @staticmethod
    def spectra_median_filter(spectra, threshold=0.1):
        '''
        
        :param series: pandas object
        :param threshold: relative value of median 
        :return: boolean indices, array of data within interval median +/- threshold
        '''
        spec = spectra.sum(axis=1)
        med = spec.median()
        ind = np.abs(1 - spec / med) < 0.1
        return ind, spectra[ind]

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
