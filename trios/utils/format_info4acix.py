import os
import pandas as pd
import numpy as np
import datetime

infof = os.path.abspath("/DATA/OBS2CO/data/info/mesures_in_situ.csv")
acixf= os.path.abspath("/DATA/OBS2CO/data/info/acix_info.csv")
info = pd.read_csv(infof, sep=';')
info.sort_values('ID_prel',inplace=True)

privID=info.ID_prel
lat=info.Lat.round(5)
lon=info.Lon.round(5)
date=pd.to_datetime(info.Date_prel).dt.strftime('%m-%d-%Y')
time=pd.to_datetime(info.h_debut)+datetime.timedelta(minutes = 15)
time=time.dt.strftime('%H.%M')
time.name='start_time_plus15min'
acix_info=pd.concat([privID,lat,lon,date,time],axis=1)
acix_info.to_csv(acixf)