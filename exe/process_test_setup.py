import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from trios.utils.sunposition import sunpos
from trios.utils import utils as u
from trios.process import *


coordf = glob.glob("/DATA/OBS2CO/data/info/mesures_in_situ.csv")[0]
coords = pd.read_csv(coordf, sep=';')
awrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/aw*idpr*.csv")

awrfiles = glob.glob("/DATA/OBS2CO/data/trios/test_setup/raw/aw*idpr*.csv")
swrfiles = glob.glob("/DATA/OBS2CO/data/trios/raw/Lu0*idpr*.csv")

idpr='167'

c = coords[coords.ID_prel == int(idpr)]  # .values[0]
lat = c['Lat'].values[0]
lon = c['Lon'].values[0]
alt = 0 #c['Altitude']
name = c['ID_lac'].values[0]

# -----------------------------------------------
#   SWR processing
# -----------------------------------------------

swr = u.swr_data(idpr, swrfiles)
if swr.file:
    df, wl = swr.reader(lat, lon, alt)
    Rrs_swr = swr_process(df, wl).process()

# -----------------------------------------------
#   AWR processing
# -----------------------------------------------
awr = u.awr_data(idpr, awrfiles)

index_idx=[2,0,1]

d=u.data(index_idx)
Ed, wl_Ed = d.load_csv(awr.Edf)
Lsky, wl_Lsky = d.load_csv(awr.Lskyf)
Lt0, wl_Lt = d.load_csv(awr.Ltf)

# ''' interpolate Ed and Lsky data upon Lt wavelength'''
wl = wl_Lt
Lt0.columns = pd.MultiIndex.from_tuples(zip(['Lt'] * len(wl), wl), names=['param', 'wl'])
intEd = interp1d(wl_Ed, Ed.values, fill_value='extrapolate')(wl)
newEd = pd.DataFrame(index=Ed.index,
                     columns=pd.MultiIndex.from_tuples(zip(['Ed'] * len(wl), wl), names=['param', 'wl']),
                     data=intEd)
intLsky = interp1d(wl_Lsky, Lsky.values, fill_value='extrapolate')(wl)
newLsky = pd.DataFrame(index=Lsky.index, columns=pd.MultiIndex.from_tuples(zip(['Lsky'] * len(wl), wl),
                        names=['param', 'wl']), data=intLsky)

awr = awr_process()
ws=[2]
fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(16, 10))
fig.subplots_adjust(left=0.1, right=0.9, hspace=.5, wspace=0.65)

i=0
for azi, Lt1 in Lt0.groupby(level=2):
    for vza,Lt in Lt1.groupby(level=1):
        ax = axs.flat[i]
        i=i+1
        print(azi,vza)

        Lsky = newLsky.loc[(newLsky.index.get_level_values(1) ==  vza) & (newLsky.index.get_level_values(2) ==  azi)]
        Ed = newEd.loc[(newEd.index.get_level_values(1) ==  vza) & (newEd.index.get_level_values(2) ==  azi)]

        Lsky_idx = Lsky.index
        Ed_idx= Ed.index
        Lt_idx = Lt.index
        Lsky.reset_index(level=[1,2],inplace=True)
        Ed.reset_index(level=[1,2],inplace=True)
        Lt.reset_index(level=[1,2],inplace=True)

        # merge sensor data on time
        df = pd.merge_asof(Lt, Ed, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")
        df = pd.merge_asof(df, Lsky, left_index=True, right_index=True, tolerance=pd.Timedelta("2 seconds"),
                           direction="nearest")

        # add solar angle data and idpr
        # compute solar angle (mean between fisrt and last aqcuisition time
        df['sza', ''] = np.nan
        for index, row in df.iterrows():
            # print index
            sza = sunpos(index, lat, lon, alt)[1]
            df.at[index, 'sza'] = sza

        rho_h = awr.get_rho_values([df.sza.min()],[vza],[azi],wl=wl)
        rho15 = awr.get_rho_mobley(awr.rhoM2015,[df.sza.min()],[vza],[azi],[ws])
        rho99 = awr.get_rho_mobley(awr.rhoM1999,[df.sza.min()],[vza],[azi],[ws])

        Rrs_h =(df.loc[:,'Lt'] -rho_h*df.loc[:,'Lsky'])/ df.loc[:,'Ed']
        Rrs15 = (df.loc[:,'Lt'] -rho15*df.loc[:,'Lsky'])/ df.loc[:,'Ed']

        Rrs99 = (df.loc[:,'Lt'] -rho99*df.loc[:,'Lsky'])/ df.loc[:,'Ed']
        #plt.figure()


        def add_curve(ax,x,mean,std,c='red',label=''):
            ax.plot(x,mean, linestyle='solid', c=c, lw=2.5,
                alpha=0.8, label=label)
            ax.fill_between(x,
                    mean - std,
                    mean + std, alpha=0.35,color=c)
        add_curve(ax,wl,Rrs_swr.transpose().mean(axis=1),Rrs_swr.transpose().std(axis=1),label='swr',c='black')
        add_curve(ax,wl,Rrs15.transpose().mean(axis=1),Rrs15.transpose().std(axis=1),label='M2015')
        add_curve(ax,wl,Rrs99.transpose().mean(axis=1),Rrs99.transpose().std(axis=1),c='orange',label='M1999')
        add_curve(ax,wl,Rrs_h.transpose().mean(axis=1),Rrs_h.transpose().std(axis=1),c='grey',label='h')

        ax.set_title('azi='+str(azi)+', vza='+str(vza))


        ax.legend(loc='best', frameon=False)

        ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
        ax.set_xlabel(r'Wavelength (nm)')

Lt.index.names



