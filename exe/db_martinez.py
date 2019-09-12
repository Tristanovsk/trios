
import sys, os
import re
import datetime as dt
import julian
dublin_julian_day = 2415020


djd = 42349.44914
jd = djd + dublin_julian_day
dt = julian.from_jd(jd)
print(dt)