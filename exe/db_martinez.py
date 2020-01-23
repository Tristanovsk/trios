
import sys, os
import re
import datetime as dt
import pandas as pd

excel_jd = [42877.57571]# 42349.44914
date = pd.TimedeltaIndex(excel_jd, unit='d') + dt.datetime(1899,12,30)
print(date)