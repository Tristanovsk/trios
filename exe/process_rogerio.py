

import glob
import re

import matplotlib as mpl
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})

from trios.utils.utils import plot as up
from trios.utils import utils as u
from trios.process import *
from trios import __package__, __version__

type_list = {'awr': 'aw', 'iwr': 'uw', 'swr': 'Lu0+'}
type_description = {'awr': 'Above-Water Radiometry', 'iwr': 'In-Water Radiometry',
                    'swr': 'Surface-Water Radiometry'}

df,Lskyf,Ltf = [os.path.join(idir,f) for f in filenames.split(' ')]
    print('process files Ed: '+ Edf+', Lsky: '+Lskyf+', Lt: '+Ltf )
    uawr = u.awr_data(idpr, Edf=Edf, Lskyf=Lskyf, Ltf=Ltf)

if uawr.Edf:
    df, wl = uawr.reader(lat, lon, alt)
    date = df.index[0].date().__str__()
    if ofile:
        ofile = os.path.join(odir, ofile)
    else:
        ofile = os.path.join(odir, 'Rrs_awr_' + date + '_idpr' + idpr + name + '.csv')

    if noclobber and os.path.exists(ofile):
        print('Skip processing: data already processed with "--no_clobber" set')
        return

    if plot:
        figfile = os.path.join(figdir, 'trios_awr_' + date + '_idpr' + idpr + name + '.png')
    else:
        figfile = ""

    awr = awr_process(df, wl, name, idpr)
    Rrs = awr.call_process(method, ofile, vza=vza, azi=azi, ws=ws, aot=aot, plot_file=figfile)