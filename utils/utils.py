import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

from utils.sunposition import sunpos


class awr_data:
    '''
    Above-water radiometry
    '''

    def __init__(self, idpr=None, files=None):
        # ''' get file names for Ed, Lsky and Lt data'''
        self.file = list(filter(lambda x: 'idpr' + idpr in x, files))
        file = self.file
        self.Edf = list(filter(lambda x: 'Ed' in x, file))
        self.Lskyf = list(filter(lambda x: 'Lsky' in x, file))
        self.Ltf = list(filter(lambda x: 'Lt' in x, file))
        self.idpr = idpr

    def reader(self, lat, lon, alt=0, name='', index_idx=[0]):
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
        :return:
        '''

        df = pd.DataFrame()

        # ''' read files with pandas format '''
        d = data(index_idx)
        Ed, wl_Ed = d.load_csv(self.Edf)
        Lsky, wl_Lsky = d.load_csv(self.Lskyf)
        Lt, wl_Lt = d.load_csv(self.Ltf)

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
        # merge sensor data on time
        df = pd.merge_asof(Lt, newEd, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")
        df = pd.merge_asof(df, newLsky, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")

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
        :param Edzf: file pat of downward in-water irradiance data
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

        df = pd.DataFrame()

        # ''' read files with pandas format '''
        d = data([1, 0])

        Ed, wl_Ed = d.load_csv(self.Edf)
        Edz, wl_Edz = d.load_csv(self.Edzf)
        Luz, wl_Luz = d.load_csv(self.Luzf)

        #mask negative values TODO save number of discarded data
        # Ed.mask(Ed<0,inplace=True)
        # Edz.mask(Edz<0,inplace=True)
        # Luz.mask(Luz<0,inplace=True)



        # copy depth data to Ed frame on date index
        # Ed.index = Ed.index.droplevel(level=1)

        #''' interpolate Ed and Lsky data upon Lt wavelength'''
        wl = wl_Luz
        Luz.columns = pd.MultiIndex.from_tuples(zip(['Luz'] * len(wl), wl), names=['param', 'wl'])
        intEd = interp1d(wl_Ed, Ed.values, fill_value='extrapolate')(wl)
        newEd = pd.DataFrame(index=Ed.index.get_level_values(0),
                             columns=pd.MultiIndex.from_tuples(zip(['Ed'] * len(wl), wl), names=['param', 'wl']),
                             data=intEd)
        intEdz = interp1d(wl_Edz, Edz.values, fill_value='extrapolate')(wl)
        newEdz = pd.DataFrame(index=Edz.index, columns=pd.MultiIndex.from_tuples(zip(['Edz'] * len(wl), wl),
                                                                                 names=['param', 'wl']), data=intEdz)

        # correct depth data for sensor to sensor distance
        Luz.reset_index(level=1, inplace=True)
        Luz.iloc[:, 0] = Luz.iloc[:, 0] + delta_Lu_depth
        # newEd.reset_index(level=1,inplace=True)

        newEdz.reset_index(level=1, inplace=True)
        newEdz.iloc[:, 0] = newEdz.iloc[:, 0] + delta_Edz_depth

        # merge sensor data on time
        df = pd.merge_asof(Luz, newEd, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
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
    #         print('Warning! Multiple files found but only one expected, process first file of the list:')
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

        # ''' interpolate Ed and Lsky data upon Lt wavelength'''
        wl = wl_Lu0
        Lu0.columns = pd.MultiIndex.from_tuples(zip(['Lu0+'] * len(wl), wl), names=['param', 'wl'])
        intEd = interp1d(wl_Ed, Ed.values, fill_value='extrapolate')(wl)
        newEd = pd.DataFrame(index=Ed.index,
                             columns=pd.MultiIndex.from_tuples(zip(['Ed'] * len(wl), wl), names=['param', 'wl']),
                             data=intEd)

        # merge sensor data on time
        df = pd.merge_asof(Lu0, newEd, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
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


class data:
    def __init__(self, index_idx=[0]):
        # first position should be datetime index
        # followed by the other parameters used for indexing (e.g. azimuth, view angle)
        self.index_idx = index_idx
        pass

    def load_csv(self, file):
        print(file)
        dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        if len(file) > 1:
            print('Warning! Multiple files found but only one expected, process first file of the list:')
            print(file)
        file_ = file[0]
        # df = pd.read_csv(file, date_parser=dateparse, sep=';', index_col=0, na_values=['-NAN'])
        df = pd.read_csv(file_, sep=';', na_values=['-NAN'])

        # get list of indexes
        col = df.columns.values[self.index_idx]
        df[col[0]] = pd.to_datetime(df[col[0]])

        df.set_index(col.tolist(), inplace=True)
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
        arr[tuple(df.index.labels)] = df[name].values.flat
        return arr
