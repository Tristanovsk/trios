import glob
import re
import datetime

from trios.utils.sunposition import sunpos
from trios.utils import utils as u
from trios.process import *

plot=True #False
odir = os.path.abspath('/DATA/OBS2CO/data/trios/surface_water')
dirfig = os.path.join(odir,'fig')

swrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/2018/Lu0*idpr*.csv")

coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ_test.csv")[0]
coords = pd.read_csv(coordf, sep=';')
coords['Date_prel']=pd.to_datetime(coords['Date_prel'])
coords['h_debut']=coords['Date_prel'] + pd.to_timedelta(coords['h_debut'])
coords['h_fin']=coords['Date_prel'] + pd.to_timedelta(coords['h_fin'])
# TODO warning: time is set as start_time + 15 minutes (can be more accurate)
coords['Date_prel']= coords['h_debut']+datetime.timedelta(minutes = 15)

iopw = ua.iopw()
iopw.load_iopw()

def add_curve(ax, x, mean, std, c='red', label=''):
    ax.plot(x, mean, linestyle='solid', c=c, lw=2.5,
            alpha=0.8, label=label)
    ax.fill_between(x,
                    mean - std,
                    mean + std, alpha=0.35, color=c)


idpr = '167'

# get idpr numbers
idprs = np.unique([re.findall(r'idpr(\d+)', x)[0] for x in swrfiles])
#idprs = np.array(['170'])
# loop over idpr
for idpr in idprs:
    c = coords[coords.ID_prel == int(idpr)]  # .values[0]
    date=c['Date_prel'].dt.strftime('%Y%m%d')
    lat = c['Lat'].values[0]
    lon = c['Lon'].values[0]
    alt = 0  # c['Altitude']
    name = c['ID_lac'].values[0]

    ofile=os.path.join(odir,'Rrs_swr_'+date.values[0]+'_idpr'+idpr+'_'+name+'.csv')
    header = c.stack(dropna=False)
    header.index = header.index.droplevel()
    header.to_csv(ofile, header=None)

    # -----------------------------------------------
    #   SWR processing
    # -----------------------------------------------

    uswr = u.swr_data(idpr, swrfiles)
    if uswr.file:
        df, wl_swr = uswr.reader(lat, lon, alt)
        df['sza', ''] = np.nan
        for index, row in df.iterrows():
            # print index
            sza = sunpos(index, lat, lon, alt)[1]
            df.at[index, 'sza'] = sza
        swr = swr_process(df, wl_swr)
        Rrs_swr = swr.call_process()
        Rrs_stat = Rrs_swr.describe()
        Rrs_stat.columns=Rrs_stat.columns.droplevel()
        Rrs_stat = Rrs_stat.T
        Rrs_stat.to_csv(ofile,mode='a')
        if plot:
            import matplotlib.pyplot as plt
            plt.rcParams.update({'font.size': 18})
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))
            fig.subplots_adjust(left=0.16, right=0.9, hspace=.5, wspace=0.65)
            add_curve(ax, wl_swr, Rrs_swr.transpose().mean(axis=1), Rrs_swr.transpose().std(axis=1), label='swr', c='black')
            Rrs_swr = swr.call_process(shade_corr=True)
            add_curve(ax, wl_swr, Rrs_swr.transpose().mean(axis=1), Rrs_swr.transpose().std(axis=1), label='swr', c='red')
            ax.legend(loc='best', frameon=False)
            ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
            ax.set_xlabel(r'$Wavelength (nm)$')
            fig.savefig(os.path.join(dirfig, 'trios_swr_'+date.values[0]+'_'+ name + '_idpr' + idpr + '.png'))
            plt.close()