import os
import numpy as np
import pandas as pd
import pandas_access as mdb
import datetime as dt
import glob
from scipy.interpolate import interp1d

from trios.config import *

wl = wl_common

opj = os.path.join

suffix = 'bagre'
DBdir = os.path.abspath('/DATA/projet/bagre/data/raw/trios')
odir = opj(DBdir, 'csv')
# DBfiles = glob.glob(opj(DBdir, '*.mdb'))
DBfile = opj(DBdir, 'trios_3107_pt3b.db')
DBfiles = glob.glob(opj(DBdir, 'Rrs*.mdb'))
DBfile = DBfiles[0]

coords = pd.read_csv(opj(DBdir, '../coord_bagre.csv'))


# isolate measurement methods
def reformat(_, wl_common):
    record = _['Data'].transform(lambda x: x.split('\r\n '))
    # transpose measurement values per wavelength
    wls = np.array([x.split(' ')[0] for x in record.values[0] if x.split(' ')[0] not in (u'', u'0')])
    data = np.array([x.split(' ')[1] for x in record.values[0] if x.split(' ')[0] in wls])
    wl = wls.astype('float')
    data = data.astype('float')

    return interp1d(wl, data, fill_value='extrapolate')(wl_common)


for DBfile in DBfiles:

    mdb.read_schema(DBfile)
    df = mdb.read_table(DBfile, 'tblData')
    query = (df.IDDataType == 'SPECTRUM') & (df.IDDataTypeSub1 == 'CALIBRATED')
    df = df[query]
    data, dates = [], []
    for date, df_ in df.groupby('DateTime'):

        print(date)
        # (here filter out in-water measurements)
        if not df_['Comment'].str.contains('deck').all():
            print(date, df_['Comment'])
            continue

        # get full acquisition set (3 meas: Lt, Lsky, Ed)
        if df_.shape[0] != 3:
            continue
        dates.append(date)
        Lt_ = reformat(df_[df_['Comment'] == 'Lu_deck'], wl_common)
        Lsky_ = reformat(df_[df_['Comment'] == 'Ld_deck'], wl_common)
        Ed_ = reformat(df_[df_['Comment'] == 'Ed_deck'], wl_common)

        data.append(np.concatenate([Lt_, Lsky_, Ed_]))

    columns = pd.MultiIndex.from_product([['Lt'], wl_common]).append(
        pd.MultiIndex.from_product([['Lsky'], wl_common])).append(
        pd.MultiIndex.from_product([['Ed'], wl_common])).set_names(['param', 'wl'])
    trios = pd.DataFrame(data, columns=columns, index=pd.Index(dates, name='date'))
    trios[trios == '-NAN'] = np.nan
    # trios= trios.dropna(axis=1)

    ID = os.path.basename(DBfile).split('-')[-1].replace('.mdb', '')
    date_ = dt.datetime.strptime(date, '%m/%d/%y %H:%M:%S').strftime('%Y%m%d')

    # reformat ID to be generic 'Pt00.'
    import re

    ID = ID.capitalize()
    ID__ = ID_ = re.search('(\d+)', ID).group(1)
    if len(ID_) == 1:
        ID_ = '0' + ID_
    ID = ID.replace(ID__, ID_)

    trios['ID'] = ID
    trios['lat'] = coords[coords.ID == ID].lat.values[0]
    trios['lon'] = coords[coords.ID == ID].lon.values[0]

    ofile = opj(odir, 'trios_' + date_ + '_' + ID + '_' + suffix + '.csv')
    trios.to_csv(ofile)
