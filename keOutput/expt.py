import numpy as np
import pandas as pd
import math
import sys
import datetime
from pathlib import Path
import glob
import os
import csv
import copy

df = pd.read_csv('outputclass.csv')
columns = ['track_id', 'x', 'y', 'class', 'timestamp']

spatdf=pd.read_csv('SPAT.csv',names=['timestamp', 'hexphase', 'binphase', 'redphase', 'cyclenumber'], index_col='timestamp')
clusdf=pd.read_csv('cluster.csv',names=['track_id','cluster','phase','class','isAnomalous'], index_col='track_id')

ndf = df.copy(True)
plist = []
cyclelist = []
cluslist = []
alist = []

prevt = -1
prevclus = 0
tracks = df['track_id']
ctracks = clusdf.index.tolist()
for t in tracks:
    if (t == prevt):
        clus = prevclus
        anom = prevanom
    else:
        if (t in ctracks):
            clus = clusdf.loc[t]['cluster']
            anom = clusdf.loc[t]['isAnomalous']
        else:
            clus = -1
            anom = 0
        prevt = t
        prevclus = clus
        prevanom = anom
    cluslist.append(clus)
    alist.append(anom)

ndf['Cluster'] = cluslist
ndf['isAnomalous'] = alist

print (ndf)

ndf.to_csv('display.csv', index=False)

