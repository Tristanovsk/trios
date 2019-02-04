import numpy as np
import pandas as pd
import scipy.interpolate as si
import plotly.plotly as py
# import plotly.graph_objs as go
from plotly.graph_objs import *

from utils.utils import reshape as r
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

    def get_rho_values(self, sza, vza, azi, ws=2, aot=0.1, wl=[550]):
        '''

        :param sza:
        :param vza:
        :param azi:
        :param ws:
        :param aot:
        :param wl:
        :return:
        '''



        grid = self.rho.rho.index.levels[2:]
        # convert pandas dataframe into 6D array of the tabulated rho values for interpolation
        rho_ = r().df2ndarray(self.rho, 'rho')

        rho_wl = calc().spline_4d(grid, rho_[1, 1, ...], (wl, sza, azi, vza))
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
        sza = self.df.loc[:, ("sza")]
        return self.process(wl, Lt, Lsky, Ed, sza, vza, azi, ws, aot)

    def process(self, wl, Lt, Lsky, Ed, sza, vza=[40], azi=[135], ws=2, aot=0.1):
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

        rho = self.get_rho_values([sza.values.mean()],vza,azi,wl=wl )
        self.Rrs = (Lt - rho * Lsky) / Ed
        self.Rrs.columns = pd.MultiIndex.from_product([['Rrs(awr)'], self.Rrs.columns], names=['param', 'wl'])

        return self.Rrs, rho


class swr_process:
    def __init__(self, df=None, wl=None, ):
        self.df = df
        self.wl = wl

    def process(self):
        wl = self.wl
        df = self.df

        Rrs = df.loc[:, ("Lu0+")] / df.loc[:, ("Ed")]
        Rrs.columns = pd.MultiIndex.from_product([['Rrs(swr)'], Rrs.columns], names=['param', 'wl'])
        self.Rrs = Rrs
        return self.Rrs


class iwr_process:
    def __init__(self, df=None, wl=None, ):
        self.df = df
        self.wl = wl

    def process(self):
        wl = self.wl
        df = self.df

        Rrs = df.loc[:, ("Lu0+")] / df.loc[:, ("Ed")]
        Rrs.columns = pd.MultiIndex.from_product([['Rrs(swr)'], Rrs.columns], names=['param', 'wl'])
        self.Rrs = Rrs
        return self.Rrs

    def Kd(self, depth, Edz):
        Kd = np.diff(Edz) / np.diff(depth)

    def plot_raw(self):
        trace = Scattergl(
            x=self.df['Luz'].values,
            y=self.df['depth_Luz'].values,

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
                range=[min(-5, min(self.df['depth_Luz'])),
                       max(0, max(self.df['depth_Luz']))],
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
        eps = si.interp1d(self.ang, self.eps_dif_EuZ)(sza)
        return eps


class calc:
    def __init__(self):
        pass

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
        import scipy.interpolate as si

        N = gin[0].__len__(), gin[1].__len__(), gin[2].__len__(), gin[3].__len__()
        Nout = gout[0].__len__(), gout[1].__len__(), gout[2].__len__()
        tmp = np.zeros([N[0], N[1], Nout[2]])

        for i in range(N[0]):
            for j in range(N[1]):
                tmp[i, j, :] = si.RectBivariateSpline(gin[2], gin[3], lut[i, j, :, :])(gout[2], gout[3], grid=False)
        if Nout[0] == Nout[1] == 1:
            interp = np.ndarray(Nout[2])
            for iband in range(Nout[2]):
                interp[iband] = si.RectBivariateSpline(gin[0], gin[1], tmp[:, :, iband])(gout[0], gout[1], grid=False)
        else:
            interp = np.ndarray([Nout[0], Nout[1], Nout[2]])
            for iband in range(Nout[2]):
                interp[:, :, iband] = si.RectBivariateSpline(gin[0], gin[1], tmp[:, :, iband])(gout[0], gout[1],
                                                                                               grid=True)

        return interp
