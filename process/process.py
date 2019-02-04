import numpy as np
import pandas as pd
from scipy import interpolate, integrate
from scipy.optimize import curve_fit

import plotly.plotly as py
# import plotly.graph_objs as go
from plotly.graph_objs import *

from utils.utils import reshape as r
import utils.auxdata as ua
from config import *


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

    def call_process(self, vza=[40], azi=[135], ws=2, aot=0.1):
        wl = self.wl
        Lt = self.df.loc[:, ("Lt")]
        Lsky = self.df.loc[:, ("Lsky")]
        Ed = self.df.loc[:, ("Ed")]
        sza = self.df.loc[:, ("sza")].values.mean()
        Rrs = self.process(wl, Lt, Lsky, Ed, sza, vza, azi, ws, aot)
        Rrs.columns = pd.MultiIndex.from_product([['Rrs(awr)'], self.Rrs.columns], names=['param', 'wl'])
        self.Rrs = Rrs
        return self.Rrs

    def process(self, wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=[2], aot=[0.1]):
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
        :return:
        '''

        rho = self.get_rho_values([sza], vza, azi, wl=wl, ws=ws, aot=aot)
        self.Rrs = (Lt - rho * Lsky) / Ed

        return self.Rrs, rho


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
        a = a + acdom + 0.4
        bb = bb + 0.05
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

    def process(self):
        wl = self.wl
        df = self.df

        reflectance = df.loc[:, ("Luz")] / df.loc[:, ("Edz")]
        reflectance.columns = pd.MultiIndex.from_product([['reflectance'], reflectance.columns], names=['param', 'wl'])
        self.reflectance = reflectance

        df['rounded_depth', ''] = df.prof_Edz.round(1)
        df.groupby('rounded_depth').mean()

        return self.reflectance

    @staticmethod
    def f_Edz(depth, Kd, Ed0):
        '''simple Edz model for homogeneous water column'''
        return Ed0 * np.exp(-Kd*depth)


    def Kd(self, depth, Edz):
        Kd = np.diff(Edz) / np.diff(depth)

    def plot_raw(self,x='Luz',y='prof_Luz'):
        trace = Scattergl(
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

        layout = Layout(
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
            margin=Margin(
                t=45,
                l=50,
                r=50
            )
        )

        return Figure(data=[trace], layout=layout)


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
        par_quanta = integrate.trapz(np.multiply(wl,Ed), wl) / hc
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
                tmp[i, j, :] = interpolate.RectBivariateSpline(gin[2], gin[3], lut[i, j, :, :])(gout[2], gout[3], grid=False)
        if Nout[0] == Nout[1] == 1:
            interp = np.ndarray(Nout[2])
            for iband in range(Nout[2]):
                interp[iband] = interpolate.RectBivariateSpline(gin[0], gin[1], tmp[:, :, iband])(gout[0], gout[1], grid=False)
        else:
            interp = np.ndarray([Nout[0], Nout[1], Nout[2]])
            for iband in range(Nout[2]):
                interp[:, :, iband] = interpolate.RectBivariateSpline(gin[0], gin[1], tmp[:, :, iband])(gout[0], gout[1],
                                                                                               grid=True)

        return interp
