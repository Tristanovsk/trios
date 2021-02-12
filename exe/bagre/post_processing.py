import os, sys
import pandas as pd
import numpy as np
import xarray as xr
import glob
import io
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.ioff()
plt.rcParams.update({'font.size': 16})

plot = False

opj = os.path.join
dir = '/DATA/projet/bagre'
dirdata = opj(dir, 'data')
aerosols = ['fine', 'coarse']
aerosol = aerosols[1]
method = 'osoaa_' + aerosol
figdir = opj(dir, 'fig/trios/')

files = glob.glob(opj(dirdata, 'L2', method, 'trios*.csv'))


def scat_angle(sza, vza, azi):
    '''
    self.azi: azimuth in rad for convention azi=0 when sun-sensenor in opposition
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


rho = pd.DataFrame()
for file in files:
    # test

    df = pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)
    rho = pd.concat([rho, df.loc[:, ['rho', 'rho_M99', 'sza', 'azi']]])

rho = rho.sort_index()
rho.columns = rho.columns.droplevel(1)
rho['doy'] = rho.index.dayofyear
rho['hour'] = rho.index.hour + rho.index.minute / 60
if plot:
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
    sns.scatterplot('hour', 'sza', hue='doy', data=rho, edgecolor="none", palette=palette, alpha=0.5, ax=axs[1],
                    legend=None)
    # sns.scatterplot('hour', 'azi', hue='doy', data=rho, edgecolor="none", palette=palette, alpha=0.5, ax=axs[1], legend=None)
    axs[1].set_xlabel('Time of day (UTC)')
    axs[1].set_ylabel('Solar zenith angle (deg)')
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

    fig1, ax1 = plt.subplots(figsize=(8, 5))
    # test on BRDF vs Rrs

norm1 = mpl.colors.Normalize(vmin=0, vmax=90)
cmap1 = plt.cm.nipy_spectral
sm1 = mpl.cm.ScalarMappable(norm=norm1, cmap=cmap1)
var = 'Rrs_osoaa_coarse'
outputs = pd.DataFrame()
for file in files:
    # test
    basename = os.path.basename(file)
    ID = basename.split('_')[2]
    df = pd.read_csv(file, header=[0, 1], index_col=0, parse_dates=True)
    df.reset_index(inplace=True)
    var_median = df[var].median()
    qmax = df[var].quantile(.68)
    qmin = df[var].quantile(0.05)
    _mean = df[var].mean(axis=1)
    _qmax = df[var].mean(axis=1).quantile(.68)
    _qmin = df[var].mean(axis=1).quantile(0.05)
    #(df[var] > qmin) & (df[var] < qmax)
    dffiltered = df[(_mean > _qmin) & (_mean < _qmax)]
    var_std = dffiltered[var].std()
    var_mean = dffiltered[var].mean()
    _df =  dffiltered.mean()
    _df['date']=dffiltered.date.mean()
    _df['ID'] = ID
    var_std = pd.concat([var_std],keys=['Rrs_std'])
    outputs=outputs.append(_df.append(var_std).T,ignore_index=True)

    # ---------------------------------
    # Plotting section
    # ---------------------------------
    if plot:
        cmap = plt.cm.get_cmap("Spectral").reversed()
        norm = mpl.colors.Normalize(vmin=400, vmax=800)
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])

        scat = df.droplevel(1, 1).apply(lambda x: scat_angle(x.sza, 40, x.azi), axis=1)
        sza, azi = df.sza, df.azi

        fig, axs = plt.subplots(ncols=3, figsize=(18, 5))
        axs = axs.ravel()
        for i, vars in enumerate([(sza, 'SZA (deg)'), (scat, 'Scattering angle (deg)')]):
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

        ax = axs[-1]
        decimal = 3
        param = 'sza'
        cmin, cmax = (df[param].min() * 0.98).round(decimal), (df[param].max() * 1.02).round(decimal)
        wl = df[var].columns.astype('float')
        norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])

        for i, group in df.iterrows():
            y = group[var].values
            if isinstance(param, str):
                color = group[param].values[0]
                stitle = param
            else:
                color = group[param]
                stitle = '_'.join(param)

            ax.plot(wl, y, color=cmap(norm(color)), lw=1.5, alpha=0.75)
        ax.plot(wl, qmin, ':k', lw=1.5)
        ax.plot(wl, var_median, '--k', lw=1.5)
        ax.plot(wl, var_mean, '-r', lw=1.5)
        ax.plot(wl, qmax, ':k', lw=1.5)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(sm, cax=cax, format=mpl.ticker.ScalarFormatter(),
                            shrink=1.0, fraction=0.1, pad=0)

        ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
        ax.set_xlabel(r'Wavelength (nm)')
        ax.set_title('Color = ' + stitle)

        plt.suptitle('Effect of in-water BRDF on Rrs')
        plt.tight_layout(rect=[0.0, 0.0, 1, 0.92])

        fig.savefig(opj(figdir, basename.split('.')[0] + '_inwater_BRDF_impact.png'), dpi=200)
        plt.close(fig)


        ax1.plot(wl, var_mean, '-', color=cmap1(norm1(color)),lw=1.5,label=basename.split('_')[2])
#
# plt.legend()
# divider = make_axes_locatable(ax1)
# cax = divider.append_axes('right', size='5%', pad=0.05)
# cbar = fig.colorbar(sm1, cax=cax, format=mpl.ticker.ScalarFormatter(),
#                     shrink=1.0, fraction=0.1, pad=0)
# fig1.savefig(opj(figdir, 'bagre_trios_QC.png'), dpi=200)

# -------------------------
# reformat and save QC data
# -------------------------

outputs.columns=pd.MultiIndex.from_tuples(outputs.columns,names=['param','wl'])
for i, columns_old in enumerate(outputs.columns.levels):
    columns_new = np.where(columns_old.str.contains('Unnamed'), '', columns_old)
    outputs.rename(columns=dict(zip(columns_old, columns_new)), level=i, inplace=True)
outputs = outputs.set_index(['ID','date','lon','lat'])
outputs = outputs.drop_duplicates().dropna(axis=1,how='all')
outputs.to_csv(opj(dirdata,'L2','QC_trios_Rrs.csv'))