import pandas as pd
import datetime as dt
from scipy.interpolate import interp1d

from trios.utils.sunposition import sunpos
from trios.config import *


class awr_data:
    '''
    Above-water radiometry
    '''

    def __init__(self, idpr=None, files=None, Edf=None, Lskyf=None, Ltf=None):
        # ''' get file names for Ed, Lsky and Lt data'''
        if not files is None:
            self.file = list(filter(lambda x: 'idpr' + idpr in x, files))
            file = self.file
            self.Edf = list(filter(lambda x: 'Ed' in x, file))
            self.Lskyf = list(filter(lambda x: 'Lsky' in x, file))
            self.Ltf = list(filter(lambda x: 'Lt' in x, file))

        elif all(v is not None for v in [Edf, Lskyf, Ltf]):
            self.Edf = Edf
            self.Lskyf = Lskyf
            self.Ltf = Ltf
        else:
            raise SyntaxError('ERROR: must specify `files` or `(Edf, Lskyf, Ltf)` variables')

        self.idpr = idpr

    def reader(self, lat, lon, alt=0, name='', index_idx=[0], utc_conv=0, file_format='csv'):
        '''
        Read above-water data files for a given acquisition series (idpr),
        merge the different data types:
          - by interpolating over wavelengths on a common band set (from those of Lt sensor)
          - by searching the nearest neighbor in time
        compute solar zenith angle
        return full data frame

        :param Edf: file path of irradiance data
        :param Lskyf: file pat of sky radiance data
        :param Ltf:  file path of water radiance data
        :param lat: latitude (decimal)
        :param lon: longitude (decimal)
        :param alt: altitude (m)
        :param idpr: ID of the acquisition series
        :param utc_conv: decimal hours added to convert local time into UTC
        :return:
        '''


        # ''' read files with pandas format '''

        d = data(index_idx, file_type=file_format)

        Ed, wl_Ed = d.load_file(self.Edf, utc_conv=utc_conv)
        Lsky, wl_Lsky = d.load_file(self.Lskyf, utc_conv=utc_conv)
        Lt, wl_Lt = d.load_file(self.Ltf, utc_conv=utc_conv)

        # ''' interpolate Ed, Lt and Lsky data upon common wavelength'''
        wl = wl_common

        intEd = interp1d(wl_Ed, Ed.values, fill_value='extrapolate')(wl)
        newEd = pd.DataFrame(index=Ed.index, columns=pd.MultiIndex.from_tuples(zip(['Ed'] * len(wl), wl),
                                                                               names=['param', 'wl']), data=intEd)

        intLt = interp1d(wl_Lt, Lt.values, fill_value='extrapolate')(wl)
        newLt = pd.DataFrame(index=Lt.index, columns=pd.MultiIndex.from_tuples(zip(['Lt'] * len(wl), wl),
                                                                               names=['param', 'wl']), data=intLt)

        intLsky = interp1d(wl_Lsky, Lsky.values, fill_value='extrapolate')(wl)
        newLsky = pd.DataFrame(index=Lsky.index, columns=pd.MultiIndex.from_tuples(zip(['Lsky'] * len(wl), wl),
                                                                                   names=['param', 'wl']), data=intLsky)

        # merge sensor data on time
        df = pd.merge_asof(newLt, newEd, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")
        df = pd.merge_asof(df, newLsky, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")
        # Convert to UTC time
        df.index = df.index  # + pd.Timedelta(hours=3)

        # add solar angle data and idpr
        # compute solar angle (mean between fisrt and last aqcuisition time
        df['sza', ''] = np.nan
        for index, row in df.iterrows():
            # print index
            sza = sunpos(index, lat, lon, alt)[1]
            df.at[index, 'sza'] = sza

        df['idpr', ''] = self.idpr
        df['name', ''] = name

        return df, wl


class iwr_data:
    '''
    In-water radiometry
    '''

    def __init__(self, idpr, files):
        # ''' get file names for Ed, Lsky and Lt data'''
        self.file = list(filter(lambda x: 'idpr' + idpr in x, files))
        file = self.file
        self.Edf = list(filter(lambda x: 'Ed_' in x, file))
        self.Edzf = list(filter(lambda x: 'Edz' in x, file))
        self.Luzf = list(filter(lambda x: 'Luz' in x, file))
        self.idpr = idpr

    def reader(self, lat, lon, alt=0, name='', delta_Lu_depth=0, delta_Edz_depth=0):
        '''
        Read above-water data files for a given acquisition series (idpr),
        merge the different data types:
          - by interpolating over wavelengths on a common band set (from those of Lt sensor)
          - by searching the nearest neighbor in time
        compute solar zenith angle
        return full data frame

        :param Edf: file path of irradiance data
        :param Edzf: file path of downward in-water irradiance data
        :param Luzf:  file path of upward in-water radiance data
        :param lat: latitude (decimal)
        :param lon: longitude (decimal)
        :param alt: altitude (m)
        :param delta_Lu_depth: adjustment of actual depth for Lu sensor (distance from depth sensor);
                               in meters for depth counted positively
        :param delta_Edz_depth: similar to delta_Lu_depth for Edz sensor
        :param idpr: ID of the acquisition series
        :return:
        '''

        # ''' read files with pandas format '''
        d = data([1, 0])

        Ed, wl_Ed = d.load_csv(self.Edf)
        Edz, wl_Edz = d.load_csv(self.Edzf)
        Luz, wl_Luz = d.load_csv(self.Luzf)

        # mask negative values TODO save number of discarded data
        Ed[Ed < 0] = 0  # .mask(Ed<0,inplace=True)
        Edz[Edz < 0] = 0  # .mask(Edz<0,inplace=True)
        Luz[Luz < 0] = 0  # .mask(Luz<0,inplace=True)

        # copy depth data to Ed frame on date index
        # Ed.index = Ed.index.droplevel(level=1)

        # ''' interpolate Ed, Edz and Luz data upon common wavelength'''
        wl = wl_common
        intEd = interp1d(wl_Ed, Ed.values, fill_value='extrapolate')(wl)
        newEd = pd.DataFrame(index=Ed.index.get_level_values(0),
                             columns=pd.MultiIndex.from_tuples(list(zip(['Ed'] * len(wl), wl)), names=['param', 'wl']),
                             data=intEd)

        intEdz = interp1d(wl_Edz, Edz.values, fill_value='extrapolate')(wl)
        newEdz = pd.DataFrame(index=Edz.index, columns=pd.MultiIndex.from_tuples(list(zip(['Edz'] * len(wl), wl)),
                                                                                 names=['param', 'wl']), data=intEdz)

        intLuz = interp1d(wl_Luz, Luz.values, fill_value='extrapolate')(wl)
        newLuz = pd.DataFrame(index=Luz.index, columns=pd.MultiIndex.from_tuples(list(zip(['Luz'] * len(wl), wl)),
                                                                                 names=['param', 'wl']), data=intLuz)

        print('read merge ok')
        # correct depth data for sensor to sensor distance
        newLuz.reset_index(level=1, inplace=True)
        newLuz.iloc[:, 0] = newLuz.iloc[:, 0] + delta_Lu_depth
        # newEd.reset_index(level=1,inplace=True)

        newEdz.reset_index(level=1, inplace=True)
        newEdz.iloc[:, 0] = newEdz.iloc[:, 0] + delta_Edz_depth

        # merge sensor data on time
        df = pd.merge_asof(newLuz, newEd, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")
        df = pd.merge_asof(df, newEdz, left_index=True, right_index=True, suffixes=('_Luz', '_Edz'),
                           tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")  # by="depth",

        # add solar angle data and idpr
        # compute solar angle (mean between fisrt and last acquisition time
        df['sza', ''] = np.nan
        for index, row in df.iterrows():
            # print index
            sza = sunpos(index, lat, lon, alt)[1]
            df.at[index, 'sza'] = sza

        df['idpr', ''] = self.idpr
        df['name', ''] = name

        return df, wl

    # def load_csv(self, file):
    #
    #     dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    #     if len(file) > 1:
    #         print('Warning! Multiple files found but only one expected, trios first file of the list:')
    #         print(file)
    #     file = file[0]
    #     df = pd.read_csv(file, sep=';', index_col=[1, 0], na_values=['-NAN'])
    #     df = df.dropna(axis=1, how='all').dropna(axis=0, how='all')
    #     df.index = df.index.set_levels([pd.to_datetime(df.index.levels[0]), df.index.levels[1]])
    #     df.columns = df.columns.astype('float')  # str.extract('(\d+)',expand=False).astype('float')
    #     # resort to get data in increasing time order
    #     df.sort_index(inplace=True, level=0)
    #     wl = df.columns
    #
    #     return df, wl


class swr_data:
    '''
    Surface-water radiometry
    '''

    def __init__(self, idpr, files):
        # ''' get file names for Ed, Lsky and Lt data'''
        self.file = list(filter(lambda x: 'idpr' + idpr in x, files))
        file = self.file
        self.Edf = list(filter(lambda x: '_Ed' in x, file))
        self.Lu0f = list(filter(lambda x: '_Lu0+' in x, file))
        self.idpr = idpr

    def reader(self, lat=None, lon=None, alt=0):
        '''
        Read above-water data files for a given acquisition series (idpr),
        merge the different data types:
          - by interpolating over wavelengths on a common band set (from those of Lt sensor)
          - by searching the nearest neighbor in time
        compute solar zenith angle
        return full data frame

        :param Edf: file path of irradiance data
        :param Lu0f:  file path of upward in-water radiance data
        :param lat: latitude (decimal)
        :param lon: longitude (decimal)
        :param alt: altitude (m)
        :param idpr: ID of the acquisition series
        :return:
        '''

        df = pd.DataFrame()

        # ''' read files with pandas format '''
        Ed, wl_Ed = data().load_csv(self.Edf)
        Lu0, wl_Lu0 = data().load_csv(self.Lu0f)

        # ''' interpolate Ed and Lsky data upon common wavelengths'''
        wl = wl_common
        intEd = interp1d(wl_Ed, Ed.values, fill_value='extrapolate')(wl)
        newEd = pd.DataFrame(index=Ed.index,
                             columns=pd.MultiIndex.from_tuples(zip(['Ed'] * len(wl), wl), names=['param', 'wl']),
                             data=intEd)
        intLu0 = interp1d(wl_Lu0, Lu0.values, fill_value='extrapolate')(wl)
        newLu0 = pd.DataFrame(index=Lu0.index, columns=pd.MultiIndex.from_tuples(zip(['Lu0+'] * len(wl), wl),
                                                                                 names=['param', 'wl']), data=intLu0)

        # merge sensor data on time
        df = pd.merge_asof(newLu0, newEd, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")

        # add solar angle data and idpr
        # compute solar angle (mean between fisrt and last aqcuisition time
        df['sza', ''] = np.nan
        for index, row in df.iterrows():
            # print index
            sza = sunpos(index, lat, lon, alt)[1]
            df.at[index, 'sza'] = sza

        df['idpr', ''] = self.idpr

        return df, wl


class fit:
    def __init__(self, N=0, m=2):
        self.popt = np.full([N, m], np.nan)
        self.pcov = np.full([N, m, m], np.nan)


class data:
    def __init__(self, index_idx=[0], file_type='csv'):
        # first position should be datetime index
        # followed by the other parameters used for indexing (e.g. azimuth, view angle)
        self.index_idx = index_idx
        self.file_type = file_type
        pass

    def load_file(self, file, utc_conv=0):

        if self.file_type == 'csv':
            return self.load_csv(file, utc_conv=utc_conv)
        elif self.file_type == 'mlb':
            return self.load_mlb(file, utc_conv=utc_conv)

    def load_csv(self, file, utc_conv=0):
        print(file)
        # dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') + pd.to_timedelta(utc_conv, 'h')
        if len(file) > 1 or not isinstance(file, str):
            print('Warning! Multiple files found but only one expected, trios first file of the list:')
            print(file)
            file_ = file[0]
        else:
            file_ = file
        # df = pd.read_csv(file, date_parser=dateparse, sep=';', index_col=0, na_values=['-NAN'])
        df = pd.read_csv(file_, sep=';|,', na_values=['-NAN'], engine='python')

        # get list of indexes
        col = df.columns.values[self.index_idx]

        # local to UTC conversion
        df[col[0]] = pd.to_datetime(df[col[0]]) + pd.to_timedelta(utc_conv, 'h')

        df.set_index(col.tolist(), inplace=True)
        df = df.dropna(axis=1, how='all').dropna(axis=0, how='all')
        df.columns = df.columns.astype('float')  # str.extract('(\d+)',expand=False).astype('float')
        # resort to get data in increasing time order
        df.sort_index(inplace=True)
        wl = df.columns
        return df, wl

    def load_mlb(self, file, utc_conv=0):

        print(file)
        # dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') + pd.to_timedelta(utc_conv, 'h')
        if len(file) > 1 and not isinstance(file, str):
            print('Warning! Multiple files found but only one expected, trios first file of the list:')
            print(file)
            file_ = file[0]
        else:
            file_ = file
        # df = pd.read_csv(file, date_parser=dateparse, sep=';', index_col=0, na_values=['-NAN'])
        header = pd.read_csv(file_, sep='\s+', skiprows=19, nrows=1)
        df = pd.read_csv(file_, sep='\s+', na_values=['-NAN'], engine='python', skiprows=21, header=None)
        df.columns = header.values[0]

        # get list of indexes
        col = self.index_idx[0]
        # local to UTC conversion
        df.iloc[:, col] = pd.TimedeltaIndex(df.iloc[:, col], unit='d') + dt.datetime(1899, 12, 30) \
                          + pd.to_timedelta(utc_conv, 'h')
        df.set_index(df.iloc[:, col], inplace=True)
        # keep spectra/radiometric data only:
        df = df.loc[:, df.columns.notnull()]

        df = df.dropna(axis=1, how='all').dropna(axis=0, how='all')
        df.columns = df.columns.astype('float')  # str.extract('(\d+)',expand=False).astype('float')
        # resort to get data in increasing time order
        df.sort_index(inplace=True)
        wl = df.columns
        return df, wl


class reshape:
    def __init__(self):
        pass

    def ndarray2df(self, arr, grid, names):
        arr = np.column_stack(list(map(np.ravel, np.meshgrid(*grid))) + [arr.ravel()])
        df = pd.DataFrame(arr, columns=names)  # e.g., names=['wind','aot','wl','sza','azi','vza','rho','rho_g'])
        return df

    def df2ndarray(self, df, name):
        shape = list(map(len, df.index.levels))
        arr = np.full(shape, np.nan)
        # fill it using Numpy's advanced indexing
        arr[tuple(df.index.codes)] = df[name].values.flat
        return arr


class plot:
    def __init__(self):
        pass

    @staticmethod
    def add_curve(ax, x, mean, std=None, c='red', label='', **kwargs):
        ax.plot(x, mean,  c=c, lw=2.5,
                alpha=0.8, label=label, **kwargs)
        if np.any(std):
            ax.fill_between(x,
                            mean - std,
                            mean + std, alpha=0.35, color=c)
