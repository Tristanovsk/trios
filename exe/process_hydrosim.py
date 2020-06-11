import os
import numpy as np
import pandas as pd
import datetime as dt
import fiona
import geopandas as gpd

# Enable fiona driver
gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
import glob

import subprocess
from multiprocessing import Pool

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.size': 18})
import plotly.graph_objects as go
import plotly.express as px
import plotly.offline as po

opb = os.path.basename
opj = os.path.join

idir = '/DATA/projet/hydrosim/data/trios'
odir = opj(idir, 'L2/awr')
figdir = opj(idir, 'fig/L2')

coordf = os.path.join(idir, '../metadata', 'Coordenadas_Rrs__2016_2019.xlsx')

metaf = os.path.join(idir, 'Metadata_hydroSim_awr.xlsx')
infofile = os.path.join(idir, 'triosdatainfo.csv')
L1dir = os.path.join(idir, 'L1')

metas = pd.read_excel(metaf)

idirs = glob.glob(L1dir + '/2*[0-9]')
idirs.sort()


def get_meta(file):
    meta = pd.read_csv(file, sep='=', nrows=18, header=None).set_index(0).T
    meta.columns = meta.columns.str.replace(' ', '')
    return meta


def generate_database_info(idirs, infofile):
    meta_attrs = ['%IDDevice', '%IntegrationTime', '%IDDataBack', '%IDDataCal', '%CalFactor']

    df = pd.DataFrame(columns=np.concatenate([['ID', 'date', 'lat', 'lon', 'wind', 'vza', 'azi',
                                               'Ed_file', 'Lsky_file', 'Lt_file'],
                                              ['Ed' + s for s in meta_attrs],
                                              ['Lsky' + s for s in meta_attrs],
                                              ['Lt' + s for s in meta_attrs]]))

    # loop on directories (for each date)
    for i, meta in metas.iterrows():
        # print(i)
        date_dir = meta[0].strftime('%Y%m%d')
        full_date = dt.datetime.combine(meta[0], meta[1])
        ID, site, comment, lat, lon, wind, vza, azi = meta[2:10]
        idir_ = os.path.join(idir, 'L1', date_dir)

        date = opb(idir_)
        print(idir_ + '/' + ID.lower() + ID.replace('id', ''))

        files = pd.Series(
            np.unique(
                np.append(glob.glob(idir_ + '/awr/' + ID.lower() + '_*.mlb'),
                          glob.glob(idir_ + '/awr/ID' + ID.replace('id', '').replace('ID', '') + '_*.mlb'))))

        if len(files) == 0:
            continue

        print('nb files ', len(files))

        for file in files:
            name = opb(file).replace('.mlb', '')
            print(name.split('_'))
            ID = name.split('_')[0]
        Edf = opj(idir_, 'awr', ID + '_Ed.mlb')
        Lskyf = opj(idir_, 'awr', ID + '_Ld.mlb')
        Ltf = opj(idir_, 'awr', ID + '_Lu.mlb')

        # read metadata from files
        Edmeta, Lskymeta, Ltmeta = get_meta(Edf), get_meta(Lskyf), get_meta(Ltf)

        # save into dataframe df
        info = np.concatenate([[ID, date, lat, lon, wind, vza, azi, opb(Edf), opb(Lskyf), opb(Ltf)], \
                               Edmeta[meta_attrs].values[0], Lskymeta[meta_attrs].values[0],
                               Ltmeta[meta_attrs].values[0]])
        df.loc[i] = info

    df.to_csv(infofile, index=False)


def plot_map( idir, datafile=''):
    ''' Plot for band 'wl' and 'method' '''
    wl = 659
    #wl = 860
    method = 'M99'

    if datafile == '':
        datafile = opj(idir, 'Rrs_stat.csv')

    df = pd.read_csv(datafile, parse_dates=True)
    df['datetime'] = df['date']

    df['date'] = pd.DatetimeIndex(df.date).normalize()

    df.set_index(['date', 'ID'], inplace=True)
    df['latitude'] = np.nan
    df['longitude'] = np.nan
    info = pd.read_csv(infofile, index_col=[1, 0], parse_dates=True)

    for idx, info_ in info.iterrows():
        print(idx)
        df.loc[idx, 'longitude'] = info_.lon
        df.loc[idx, 'latitude'] = info_.lat

    df.reset_index(inplace=True)
    df['month'] = df['date'].dt.month
    df['date'] = df['date'].dt.strftime('%Y-%m-%d')

    df_ = df.loc[(df.wl == wl) & (df.method == method)]
    df_['size'] = 8

    fig = px.scatter_mapbox(df_, lat="latitude", lon="longitude", color='0.5', hover_name="ID",
                            hover_data=["date", 'method'],
                            size='size',
                            range_color=[0, df_['0.5'].max() * 0.7],
                            animation_frame='month',
                            color_continuous_scale=px.colors.diverging.RdYlBu_r,
                            zoom=5.5, opacity=0.65, height=900, width=1100,
                            title='Sample points; Rrs at ' + str(wl) + ' nm')
    # fig.update_traces(marker=dict(line_width=2),selector=dict(mode='markers'))
    fig.update_traces(

        line=dict(
            width=0.5,

        ), )
    fig.update_layout(mapbox_style="open-street-map")
    fig.update_layout(title_font_size=24,
                      # mapbox=go.layout.Mapbox(zoom=6,pitch=10),
                      mapbox_style="white-bg",
                      mapbox_layers=[
                          {
                              "below": 'traces',
                              "sourcetype": "raster",
                              "source": [
                                  "https://services.arcgisonline.com/arcgis/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
                              ]
                          },

                      ])
    # fig.update_layout(margin={"r":10,"t":0,"l":10,"b":0})
    po.plot(fig, filename=opj(figdir, 'map_Rrs' + str(wl) + '.html'))


def call(command):
    print(command)
    cp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if cp.returncode != 0:
        print('ERROR for ', command.split(' ')[1:3])
        with open('err.txt', 'a') as ferr:
            ferr.write(command)

    return


def process(infofile, method='M99', ncore=10):
    info = pd.read_csv(infofile)
    utc_conv = '4'

    commands = []

    for i, info_ in info.iterrows():

        ID = str(info_.ID)
        date = str(info_.date)

        print(date, ID)

        odir_ = os.path.join(odir, method, date)
        if not os.path.exists(odir_):
            os.makedirs(odir_)

        lat, lon,vza,azi,wind = str(info_.lat), str(info_.lon),str(info_.vza), str(info_.azi),str(info_.wind)

        Edf, Lskyf, Ltf = info_.Ed_file, info_.Lsky_file, info_.Lt_file
        command = 'trios_processing ' + os.path.join(L1dir, date,'awr') + ' ' + ID + \
                  ' awr --lat ' + lat + ' --lon ' + lon + ' --name _' + method + \
                  ' --vza ' + vza + ' --azi ' + azi + ' --ws ' + wind + \
                  ' --data_files "' + Edf + ' ' + Lskyf + ' ' + Ltf + \
                  '" --odir ' + odir_ + ' --utc_conv=' + utc_conv + \
                  ' --method ' + method + ' --plot --figdir ' + figdir + ' --format mlb --no_clobber'
        commands.append(command)
        print(command)
    with Pool(processes=ncore) as pool:

        print(pool.map(call, commands))
        pool.close


def post_process(idir,figdir,create_stat_file=False):

    stat_file = os.path.join(idir, 'all/Rrs_stat.csv')

    if create_stat_file:
        odf = pd.DataFrame()
        for method in ['M99', 'osoaa', 'temp_opt']:
            files = glob.glob(os.path.join(idir, method) + '/2*[0-9]/Rrs*csv')
            for file in files:
                info = file.split('_')
                ID = info[-2].replace('idpr', '')
                if method == 'temp_opt':
                    ID = info[-3].replace('idpr', '')
                df = pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)
                date = df.index.mean()
                Rrs = df.Rrs.quantile([0.25, 0.5, 0.75]).T
                wl = Rrs.index.values
                Rrs.index = pd.MultiIndex.from_product([[date], [ID], [method], wl],
                                                       names=['date', 'ID', 'method', 'wl'])

                odf = pd.concat([odf, Rrs])

        # save data file
        odf.to_csv(stat_file)

    df = pd.read_csv(stat_file, index_col=[0], parse_dates=True)
    df.sort_values(['date', 'ID', 'method', 'wl'], inplace=True)

    import matplotlib.pyplot as plt
    plt.ioff()
    ncols = 4

    obj = df.groupby(['date', 'ID'])
    N = obj.groups.__len__()
    rows = N // ncols + (1 if N % ncols else 0)
    aspect_ratio = 1.3 * ncols
    fig, axs = plt.subplots(nrows=rows, ncols=ncols, figsize=(25, rows * aspect_ratio))
    fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.45)
    c = dict(M99='red', osoaa='black', temp_opt='blue')

    for (idx, group), ax in zip(obj, axs.flatten()):
        print(group.ID.unique())

        for method, g in group.groupby('method'):
            ax.plot(g.wl, g['0.5'], label=method, color=c[method])
            ax.fill_between(g.wl, g['0.25'], g['0.75'], alpha=0.35, color=c[method])
        ax.legend(loc='best', frameon=False)
        ax.set_xlabel(r'Wavelength (nm)')
        ax.set_title('ID ' + idx[1] + ', ' + idx[0].date().__str__())
    plt.tight_layout()
    plt.savefig(figdir, bbox_inches='tight',dpi=300)
    plt.close()


# ----------------------------
# to execute
# ----------------------------

# plot_map(coords)
# generate_database_info(idirs,infofile)
method = 'M99'
process(infofile, method=method, ncore=10)

for method in ['M99', 'osoaa', 'temp_opt']:
    process(infofile, method=method, ncore=10)

post_process(odir,figdir,create_stat_file=True)

# ----
# END
