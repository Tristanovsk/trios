''' Executable to process L1C images from Sentinel-2 and Landsat mission series

Usage:
  trios_processing <input_dir> <IDpr> <measurement_type> --lat <lat> --lon <lon> \
   [--altitude=alt] [--ofile <ofile>] [--odir <odir>] [--plot] [--figdir <figdir>] \
   [--name <name>] [--method <method>] [--no_clobber] \
   [--vza <vza>] [--azi <azi>] [--ws <ws>] [--aot <aot>]
  trios_processing -h | --help
  trios_processing -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <input_dir>  Input directory where the trios.csv files are located
  <IDpr>  IDpr label of the measurement sequence to process
  <measurement_type>   awr, swr, iwr
  --lat <lat>  Latitude in decimal degrees
  --lon <lon>  Longitude in decimal degrees

  --altitude=alt   altitude of the scene to be processed in meters
                   [default: 0]
  --odir odir   Path to the output directory [default: ./]
  --ofile ofile       basename of the output file.
  --plot  Plot output data and save figure in <figdir>
  --figdir figdir  Directory where figures are saved [default: ./]
  --name name  Keyword to append to file name and figures [default: ]
  --method method  Keyword for the method to apply for data processing.
                   For awr: M99, M15, osoaa, temp_opt [default: M99]
  --vza vza  Viewing zenith angle (in deg) of the awr radiometer (applicable for all awr methods)
            [default: 40]
  --azi azi  Relative azimuth angle (in deg) of the awr radiometer (applicable for all awr methods)
            [default: 135]
  --ws ws  Wind speed (in m/s) (applicable for awr methods: M99, M15, osoaa)
           [default: 2]
  --aot aot  Aerosol optical thickness at 550 nm (applicable for awr methods: osoaa)
             [default: 0.1]
  --no_clobber     Do not process  <input_dir> <IDpr> files if <output_file> already exists.
'''

from docopt import docopt
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


def main():
    args = docopt(__doc__, version=__package__ + ' ' + __version__)
    print(args)

    idir = os.path.abspath(args['<input_dir>'])
    idpr = args['<IDpr>']
    meas_type = args['<measurement_type>']
    method = args['--method']
    lat = float(args['--lat'])
    lon = float(args['--lon'])
    alt = float(args['--altitude'])
    odir = os.path.abspath(args['--odir'])
    ofile = args['--ofile']
    name = args['--name']
    plot = args['--plot']
    figdir = os.path.abspath(args['--figdir'])
    noclobber = args['--no_clobber']
    azi = float(args['--azi'])
    vza = float(args['--vza'])
    ws = float(args['--ws'])
    aot = float(args['--aot'])

    try:
        type_ = type_list[meas_type]
    except:
        raise SyntaxError('ERROR: bad request for <measurement_type>, should be either awr, iwr or swr')

    files = glob.glob(os.path.join(idir, type_ + '*' + idpr + '*.csv'))
    if not files:
        raise IOError('No file available for IDpr ' + idpr + ' in ' + idir + ' for ' + type_description[meas_type])

    if meas_type == 'swr':
        # -----------------------------------------------
        #   SWR processing
        # -----------------------------------------------
        uswr = u.swr_data(idpr, files)
        if uswr.file:
            df, wl = uswr.reader(lat, lon, alt)
            date = df.index[0].date().__str__()
            if ofile:
                ofile = os.path.join(odir, ofile)
            else:
                ofile = os.path.join(odir, 'Rrs_swr_' + date + '_idpr' + idpr + name + '.csv')

            if noclobber and os.path.exists(ofile):
                print('Skip processing: data already processed with "--no_clobber" set')
                return

            swr = swr_process(df, wl)
            Rrs = swr.call_process(ofile)

            if plot:
                mpl.rcParams.update({'font.size': 18})
                fig, ax = plt.subplots(figsize=(7, 6))
                up.add_curve(ax, wl, Rrs.transpose().mean(axis=1), Rrs.transpose().std(axis=1), label='swr',
                             c='black')
                ax.set_ylabel(r'$R_{rs}\  (sr^{-1})$')
                ax.set_xlabel(r'Wavelength (nm)')
                ax.set_title('ID: ' + idpr + ', ' + date + ', sza=' + str(round(df.sza.mean(), 2)))
                fig.savefig(os.path.join(figdir, 'trios_swr_' + date + '_idpr' + idpr + name + '.png'),
                            bbox_inches='tight')
                plt.close()

    elif meas_type == 'awr':
        # -----------------------------------------------
        #   AWR processing
        # -----------------------------------------------

        uawr = u.awr_data(idpr, files)
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

    elif meas_type == 'iwr':
        # -----------------------------------------------
        #   IWR processing
        # -----------------------------------------------

        uiwr = u.iwr_data(idpr, files)
        if uiwr.file:
            df, wl = uiwr.reader(lat, lon, alt)
            date = df.index[0].date().__str__()
            if ofile:
                ofile = os.path.join(odir, ofile)
            else:
                ofile = os.path.join(odir, 'Rrs_iwr_' + date + '_idpr' + idpr + name + '.csv')

            if noclobber and os.path.exists(ofile):
                print('Skip processing: data already processed with "--no_clobber" set')
                return

            if plot:
                figfile = os.path.join(figdir, 'trios_iwr_' + date + '_idpr' + idpr + name + '.pdf')
            else:
                figfile = ""

            iwr = iwr_process(df, wl, name, idpr)
            results = iwr.call_process(ofile, plot_file=figfile)


if __name__ == "__main__":
    main()
