import os, sys
import pandas as pd
import numpy as np
import xarray as xr
import glob
import io
import matplotlib as mpl
jimport matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams.update({'font.size': 16})
opj = os.path.join
dir = '/DATA/projet/borges'
dirdata = opj(dir, 'data')
aerosols = ['fine', 'coarse']
aerosol = aerosols[1]
method = 'osoaa_' + aerosol

odir = opj(dirdata, 'L2', method)

files = glob.glob(opj(odir, 'awr_L2_*.csv'))
rho = pd.DataFrame()
for file in files:
    # test

    df = pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)
    rho = pd.concat([rho, df.loc[:, ['rho', 'rho_M99', 'sza', 'azi']]])

rho = rho.sort_index()
rho.columns = rho.columns.droplevel(1)
rho['doy'] = rho.index.dayofyear
rho['hour'] = rho.index.hour + rho.index.minute / 60

import seaborn as sns

colors = rho.doy.unique()
palette = dict(zip(colors,
                   sns.color_palette("rocket_r", len(colors))))
fig, axs = plt.subplots(ncols=2, figsize=(16, 6))
sns.scatterplot('hour', 'rho', hue='doy', data=rho, alpha=0.5, edgecolor="none", ax=axs[0], palette=palette)
sns.scatterplot('hour', 'rho_M99', hue='doy', data=rho, marker="+", alpha=0.5, ax=axs[0], legend=None)
axs[0].set_xlabel('Time of day (UTC)')
axs[0].set_ylabel('Reflection factor')
plt.legend(fontsize=11)
sns.scatterplot('hour', 'sza', hue='doy', data=rho, marker="+", alpha=0.5, ax=axs[1], legend=None)
sns.scatterplot('hour', 'azi', hue='doy', data=rho, marker="+", alpha=0.5, ax=axs[1], legend=None)

plt.savefig(opj(dir, 'fig', method + '_rho_values.png'))


# Rrs_basic = (df.Lt-0.028*df
# Rrs_basic.columns=pd.MultiInde.Lsky)/df.Ed
# wl = Rrs_basic.columnsx.from_product([['Rrs_basic'],wl])
#
# rel_diff = (Rrs_basic -df.Rrs)/df.Rrs
# rel_diff.columns=pd.MultiIndex.from_product([['rel_diff'],wl])
# rel_diff[np.abs(rel_diff)> 10]=np.nan
# df_tot = df.join(Rrs_basic).join(rel_diff)
# df_tot.to_csv(file.replace('.csv','_test.csv'))


# test on BRDF vs Rrs

def scat_angle(sza, vza, azi):
    '''
    self.azi: azimuth in rad for convention azi=0 when sun-sensenor in oppositioon
    :return: scattering angle in deg
    '''
    print(sza, vza, azi)
    sza = np.pi / 180 * sza
    vza = np.pi / 180 * vza
    azi = np.pi / 180 * azi
    ang = -np.cos(sza) * np.cos(vza) + np.sin(sza) * np.sin(vza) * np.cos(azi)
    # ang = np.cos(np.pi - sza) * np.cos(vza) - np.sin(np.pi - sza) * np.sin(vza) * np.cos(azi)
    ang = np.arccos(ang) * 180 / np.pi

    return ang

plt.ioff()
file = opj(odir, 'awr_L2_osoaa_coarse_2019-04-24_19:00:00.csv')
df = pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)
file = opj(odir, 'awr_L2_osoaa_coarse_2019-04-24_18:00:00.csv')
df = pd.concat([df, pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)])
file=opj(odir,'awr_L2_osoaa_coarse_2019-04-28_12:00:00.csv')
df= pd.concat([df,pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)])
file=opj(odir,'awr_L2_osoaa_coarse_2019-04-25_11:00:00.csv')
df= pd.concat([df,pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)])

cmap = plt.cm.get_cmap("Spectral").reversed()
norm = mpl.colors.Normalize(vmin=400, vmax=800)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
var = 'Rrs_osoaa_coarse'
scat = df.droplevel(1, 1).apply(lambda x: scat_angle(x.sza, 40, x.azi), axis=1)
sza, azi = df.sza, df.azi

fig, axs = plt.subplots(ncols=3, figsize=(18, 5))
axs = axs.ravel()
for i, vars in enumerate([(sza, 'SZA (deg)'), (azi, 'Azimuth (deg)'), (scat, 'Scattering angle (deg)')]):
    x = vars[0]
    xlabel = vars[1]
    for iwl, group in df[var].iloc[:, ::4].iteritems():
        wl = float(iwl)
        if wl > 800:
            continue

        y = group.values
        axs[i].plot(x, y, 'o', color=cmap(norm(wl)), lw=0.5, markersize=2, alpha=0.5)
    divider = make_axes_locatable(axs[i])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(sm, cax=cax, format=mpl.ticker.ScalarFormatter(),
                        shrink=1.0, fraction=0.1, pad=0)

    axs[i].set_ylabel('Rrs $(sr^{-1})$')
    axs[i].set_xlabel(xlabel)

plt.suptitle('Effect of in-water BRDF on Rrs')
plt.tight_layout(rect=[0.0, 0.0, 1, 0.92])

fig.savefig(opj(dir, 'fig', 'Rrs_inwater_BRDF_impact.png'), dpi=200)
