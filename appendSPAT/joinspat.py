import numpy as np
import pandas as pd
import math
import sys
import datetime
from pathlib import Path
import glob
import os
import csv
from hash import HashTable
import copy

df = pd.read_csv('outputclass.csv')
columns = ['track_id', 'x', 'y', 'class', 'timestamp']

spatdf=pd.read_csv('SPAT.csv',names=['timestamp', 'hexphase', 'binphase', 'redphase', 'cyclenumber'], index_col='timestamp')

ndf = df.copy(True)

spatH = HashTable(86400)
plist = []
rlist = []

prevt = -1
prevphase = 0
start = 0
timestamps = df['timestamp']
for t in timestamps:
    # find spat from SPAT
    if (t == prevt):
        phase = prevphase
    else:
        phase = spatdf.loc[t]
        prevt = t
        prevphase = phase
    plist.append(phase['hexphase'])
    rlist.append(phase['redphase'])


ndf['SPAT'] = plist
ndf['RedStatus'] = rlist
print (ndf)

ndf.to_csv('combo.csv', index=False)

