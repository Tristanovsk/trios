# -*- coding: utf-8 -*-

"""
----------------------------------------------------------------
PYTHON SCRIPT
----------------------------------------------------------------
DBtrios.py
Created on: may 2018
Last update: 01/02/2019
Author: Nathalie REYNAUD

OBS2CO

Post-traitement des données TRIOS récoltées dans le cadre du
projet TELQUEL/OBS2CO - facultatif : ajout des valeurs de
profondeur tirées du capteur de pression HOBO

----------------------------------------------------------------
requirements warning: pandas
----------------------------------------------------------------
"""

""" Required libraries """

import os
import pyodbc
import pandas as pd
import glob

""" Environment """

""" Local variables """

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<I/O>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# répertoire projet
mydir = r'Z:\Hyax_Travail\LACS\TELQUEL'
# répertoire des données TRIOS
DBrep = r'Donnees\Donnees_Brutes\Terrain\TRIOS'

hobo_process = True
# dictionnaire du décalage (en mètres) entre les radiomètres et le capteur de pression
dec_m = {'Luz': 0.07, 'Edz': -0.28, 'reflectancez': -0.105}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<_>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# base de données TRIOS
DBfile = glob.glob(os.path.join(mydir, DBrep, 'dataTriOS_*.mdb'))[1]

""" Fonctions """


def odbc2lst(conn, query):
    cursor = conn.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    return result


""" Processing """

'''
STEP 1: trios TRIOS data
'''

# connect to bd TRIOS
odbc = pyodbc.connect('DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=' + DBfile)

# retrieve RAMSES data
query = 'SELECT * FROM tblData WHERE ((tblData.IDDataType LIKE \'SPECTRUM\') AND ' \
        '((tblData.IDDataTypeSub1 LIKE \'CALIBRATED\') OR (tblData.IDDataTypeSub1 LIKE \'CALCULATED\')))'
ramses_df = pd.read_sql(query, odbc)

# isolate measurement methods
dicmeth = {'Lu0+': ramses_df[ramses_df['CommentSub1'] == 'Lu0+'],
           'uw': ramses_df[ramses_df['CommentSub1'] == 'underwater'],
           'aw': ramses_df[ramses_df['CommentSub1'] == 'above water']}

# for each measurement method
for key in dicmeth.keys():
    lst_idpr = []
    # group by parameters and samples
    groups = dicmeth[key].groupby(['MethodName', 'CommentSub2', 'CommentSub3'])
    for name, group in groups:
        print(name)
        # retrieve and format data
        idpr = name[2].replace(' ', '').split('pr')[-1]
        # add idpr to list
        if idpr not in lst_idpr:
            lst_idpr.append(idpr)
        # check for retry
        if 'bis' in group['Comment'].values:
            group = group[group['Comment'] == 'bis']
        dt = group['DateTime']
        mes = group['Data'].transform(lambda x: x.split('\r\n '))

        # transpose measurement values per wavelength
        wls = [x.split(' ')[0] for x in mes.values[0] if x.split(' ')[0] not in (u'', u'0')]
        vals = []
        for row in mes:
            vals.append([x.split(' ')[1] for x in row if x.split(' ')[0] in wls])
        mes_df = pd.DataFrame(vals, columns=wls, index=mes.index)
        exp_df = pd.concat([dt, mes_df], axis=1, join='outer')
        # add depth information for underwater measurements
        if key == 'uw':
            z = name[1][2:]
            depth = pd.DataFrame([z] * len(exp_df), columns=['depth'], index=mes.index)
            exp_df = pd.concat([depth, exp_df], axis=1, join='outer')

        # save result to csv
        out_file = os.path.join(mydir, DBrep,
                                '{0}_{1}_idpr{2}.csv'.format(key, name[0], idpr))
        if 'SAM' in name[0]:
            out_file = os.path.join(mydir, DBrep,
                                    '{0}_{1}_{3}_idpr{2}.csv'.format(key, name[1], idpr, name[0].replace('_', ''), ))
        elif key == 'uw':
            out_file = os.path.join(mydir, DBrep,
                                    '{0}_{1}_idpr{2}.csv'.format(key, name[0] + z, idpr))
        exp_df.to_csv(out_file, sep=';', index=False)
    # compile data for underwater measurements at various depths in a single csv file
    if key == 'uw':
        for par in ('Ed', 'Lu', 'reflectance'):
            for idpr in lst_idpr:
                flst = glob.glob(os.path.join(mydir, DBrep,
                                              'uw_{0}[0,-][-,0-9]*{1}.csv'.format(par, idpr)))
                try:
                    foo = os.path.basename(flst[0]).split('_')[1]
                    comp_file = flst[0].replace(foo, foo.strip('-.0123456789') + 'z')
                    comp_df = pd.concat([pd.read_csv(f, sep=';') for f in flst])
                    comp_df.to_csv(comp_file, sep=';', index=False)
                    for f in flst:
                        os.remove(f)
                except IndexError:
                    pass

'''
STEP 2: trios HOBO data
'''

if hobo_process:

    # for each "underwater" RAMSES file
    uwfiles = glob.glob(os.path.join(mydir, DBrep, 'csv', 'uw_*z_*.csv'))
    for f in uwfiles:
        # find corresponding hobo file
        idpr = f[-7:-4]
        mes = os.path.basename(f).split('_')[1]
        hobo = glob.glob(os.path.join(mydir, r'Donnees\Donnees_Brutes\Terrain\*\TRIOS\HOBO_pression',
                                      '*_idpr{0}*.txt'.format(idpr)))
        # print message on error
        if len(hobo) == 0:
            print('ERREUR sur idpr{0} : fichier hobo manquant !'.format(idpr))
        elif len(hobo) > 1:
            print('ERREUR sur idpr{0} : plusieurs fichiers hobo trouvés.' \
                  'Vérifiez les idpr et/ou supprimer les fichiers superflus !'.format(idpr))
        else:
            # load data from RAMSES and HOBO
            dfh = pd.read_csv(hobo[0], sep=';', header=1, usecols=range(5), decimal=',',
                              names=['id', 'date', 'time', 'pres', 'temp'])
            dff = pd.read_csv(f, sep=';')
            # calculate depth from pressure (p0 = 1st pressure value after launching sensor)
            # formule: z = (p-p0)/(ro*g) avec ro = 0.999 (environ : dépend de la température) et g = 9.81 (standard)
            p0 = dfh['pres'][0]
            dfh = dfh.assign(prof=(dfh['pres'] - p0) / (0.999 * 9.81) + dec_m[mes])
            dt_df = pd.DataFrame(dfh['date'].str.replace('/', '-') + ' ' + dfh['time'],
                                 columns=['DateTime'], index=dfh.index)
            dfh = pd.concat([dfh, dt_df], axis=1, join='outer')
            # join data on datetime
            all_df = pd.merge(dff, dfh, how='left', on='DateTime')
            # export final results to csv
            fdf = pd.concat([all_df['prof'], all_df.iloc[:, 1:-7]], axis=1, join='outer')
            fdf.to_csv('_hobo.'.join(f.split('.')), sep=';', index=False)

""" ------------------------------END---------------------------------- """
