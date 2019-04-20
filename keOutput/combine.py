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
    cyclelist.append(phase['cyclenumber'])

ndf['SPAT'] = plist
ndf['Cycle'] = cyclelist

tracks = df['track_id'].unique()

slistx=[]
slisty=[]

p=0
for tid in tracks:
    nndf = df.loc[df['track_id'] == tid]
    x = nndf.iloc[:,2].values
    y = nndf.iloc[:,3].values
    f = nndf.iloc[:,0].values
    num = len(x)
    for i in range (0, num):
        if (i == 0 or i == 1):
            slistx.append(0)
            slisty.append(0)
        else:
            speedx = 0.5 * ((x[i-2] - x[i-1])/(f[i-2] - f[i-1]) + (x[i-1] - x[i])/(f[i-1] - f[i]))
            speedy = 0.5 * ((y[i-2] - y[i-1])/(f[i-2] - f[i-1]) + (y[i-1] - y[i])/(f[i-1] - f[i]))
            slistx.append(speedx)
            slisty.append(speedy)

ndf['speed_x'] = slistx
ndf['speed_y'] = slisty

alist=[]

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

