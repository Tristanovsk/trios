# !/usr/bin/env python
import os
import numpy as np
import cmocean
import matplotlib.pyplot as plt

from process.process import *

dirfig = os.path.abspath('/DATA/projet/gilerson/fig')
awr = awr_process()
ws = 2
aot = 0.1

szas = [10, 30, 50, 70]
Nx, Ny = 50, 50
delta = 0.8
azi_center = 90

vza = np.linspace(20, 60, Ny, endpoint=True)
azi = np.linspace(azi_center - 20, azi_center + 20, Nx, endpoint=True)

VZA, AZI = np.meshgrid(vza, azi)
VZAr, AZIr = np.radians(VZA), np.radians(AZI),

fig, axs = plt.subplots(2, 2, figsize=(12, 9))
fig.subplots_adjust(left=0.1, right=0.9, hspace=.35, wspace=0.25)
origin = 'lower'

for ax, sza in zip(axs.ravel(), szas):
    szar = np.radians(sza)
    scatt = np.arccos(-np.cos(VZAr) * np.cos(szar) + np.sin(VZAr) * np.sin(szar) * np.cos(AZIr))
    scatt = np.degrees(scatt)

    rho = np.empty((Nx, Ny))
    for ivza in range(Ny):
        rho[..., ivza] = awr.get_rho_values([sza], [vza[ivza]], azi, wl=[550], sunglint=True)

    rho = np.ma.array(rho)
    # mask pixels:
    # rho[...] = np.ma.masked

    CS = ax.contourf(AZI, VZA, rho, 8,  # [-1, -0.1, 0, 0.1],
                     # alpha=0.5,
                     cmap=cmocean.cm.thermal,
                     origin=origin)
    Cscatt = ax.contour(AZI, VZA, scatt, 10, cmap=plt.cm.gray)
    ax.clabel(Cscatt, Cscatt.levels[1::2],  # label every second level
              inline=1,
              fmt='%1.1f',
              fontsize=14)

    # CS2 = ax.contour(CS, levels=CS.levels[::2], colors = 'r', origin=origin,hold='on')

    ax.set_title('sza=' + str(sza) + ', wind=' + str(ws))
    ax.set_ylabel('Viewing angle (deg)')
    ax.set_xlabel('Azimuth (deg)')

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = fig.colorbar(CS, ax=ax, shrink=0.9)

    cbar.ax.set_ylabel(r'$\rho-factor$')
fig.savefig(os.path.join(dirfig, 'rho_g_snapshot_azi' + str(azi_center) + '_aot' + str(aot) + '_ws' + str(ws) + '.png'))
plt.close()
