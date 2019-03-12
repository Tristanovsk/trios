'''where you can set absolute and relative path used in the package'''
import os
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))

M2015_file = os.path.join(root, '../aux/rhoTable_Mobley2015.csv')
M1999_file = os.path.join(root, '../aux/rhoTable_Mobley1999.csv')
rhosoaa_fine_file = os.path.join(root, '../aux/surface_reflectance_factor_rho_fine_aerosol_rg0.06_sig0.46.csv')
rhosoaa_coarse_file = os.path.join(root, '../aux/surface_reflectance_factor_rho_coarse_aerosol_rg0.60_sig0.60.csv')
iopw_file =  os.path.join(root, '../aux/water_coef.txt')
F0_file =  os.path.join(root, '../aux/Thuillier_2003_0.3nm.dat')

# set common wavelengths on which spectra are reprojected
# default 320nm to 950 nm each wl_step nm
wl_step = 3
wl_common = np.arange(320,950+wl_step,wl_step)

