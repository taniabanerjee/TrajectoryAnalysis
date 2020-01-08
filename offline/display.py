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
import pymysql

if (len(sys.argv) == 1):
        print ('Error: provide confile file')

configfile = sys.argv[1]

options = {}
def extractOptions(fp):
   line = fp.readline().rstrip('\n')
   while (line != ''):
       strlist = line.split()
       options[strlist[0]] = strlist[1]
       if (len(strlist) > 2):
           options[strlist[0]] = strlist[1] + ' ' + strlist[2]
       line = fp.readline().rstrip('\n')
   line = fp.readline().rstrip('\n')
   return line

filepath = configfile
with open(filepath) as fp:
   line = fp.readline().rstrip('\n')
   while (line):
       line = extractOptions(fp)

df = pd.read_csv(options['trackinfo'])
clusdf=pd.read_csv(options['clusterf'],names=['track_id','cluster','isAnomalous','startPoint','endPoint','phase','startphase','endphase','cycle'], index_col='track_id')
#rdf = pd.read_csv('../DataFiles/redjump.csv')
ctracks = clusdf.index.tolist()
#rtracks = rdf['track_id'].unique()
rtracks = []

plist = []
phlist = []
sphlist = []
ephlist = []
cyclelist = []
cluslist = []
alist = []

tracks = df['track_id'].unique()

nx=[]
ny=[]
tracklist=[]
slistx=[]
slisty=[]
framelist=[]
classlist=[]
timelist=[]
intersectionlist=[]
cluslist=[]
anomlist=[]
plist = []
cyclelist = []
rjlist = []
nmisslist = []
siglist = []

p=0
for tid in tracks:
    nndf = df.loc[df['track_id'] == tid]
    x = nndf.iloc[:,3].values
    y = nndf.iloc[:,4].values
    f = nndf.iloc[:,1].values
    c = nndf.iloc[:,5].values
    t = nndf.iloc[:,6].values
    ts = nndf.iloc[:,0].values
    intersection = nndf.iloc[:,7].values
    nm = nndf.iloc[:,10].values
    sig = nndf.iloc[:,11].values
    y = np.ma.compressed(np.ma.masked_where(x==0, y))
    f = np.ma.compressed(np.ma.masked_where(x==0, f))
    c = np.ma.compressed(np.ma.masked_where(x==0, c))
    t = np.ma.compressed(np.ma.masked_where(x==0, t))
    ts = np.ma.compressed(np.ma.masked_where(x==0, ts))
    intersection = np.ma.compressed(np.ma.masked_where(x==0, intersection))
    nm = np.ma.compressed(np.ma.masked_where(x==0, nm))
    sig = np.ma.compressed(np.ma.masked_where(x==0, sig))
    x = np.ma.compressed(np.ma.masked_where(x==0, x))
    if (x.size < 3):
        continue
    diffx = x[x.size-1]-x[0]
    diffy = y[y.size-1]-y[0]
    if (diffx > 10):
        x = np.maximum.accumulate(x)
    elif (diffx < -10):
        x = np.minimum.accumulate(x)
    if (diffy > 10):
        y = np.maximum.accumulate(y)
    elif (diffy < -10):
        y = np.minimum.accumulate(y)

    #speedx = nndf.iloc[:,8].values
    #speedy = nndf.iloc[:,9].values
    #sp = nndf.iloc[:,11].values
    #cycle = nndf.iloc[:,12].values
    num = len(x)
    startP = 0
    endP = num
    clus = 'anom'
    anom = 1
    ph = '4400bb'
    sph = '4400bb'
    eph = '4400bb'
    speedx = 0
    speedy = 0
    cycle = 0
    rj = 0
    if (tid in ctracks):
        #startP = clusdf.loc[tid]['startPoint']
        #endP = clusdf.loc[tid]['endPoint']
        clus = clusdf.loc[tid]['cluster']
        anom = clusdf.loc[tid]['isAnomalous']
        ph = clusdf.loc[tid]['phase']
        sph = clusdf.loc[tid]['startphase']
        eph = clusdf.loc[tid]['endphase']
        cycle = clusdf.loc[tid]['cycle']
    if (tid in rtracks):
        rj = 1
    for i in range (startP, endP):
        #frame_id,track_id,center_x,center_y,class,timestamp,intersection_id,SPAT,Cycle,speed_x,speed_y,Cluster,isAnomalous
        nx.append(x[i])
        ny.append(y[i])
        framelist.append(f[i])
        rjlist.append(rj)
        timelist.append(t[i]+'.'+str(f[i]%10))
        nmisslist.append(nm[i])
        siglist.append(sig[i])
        
        if (i == 0):
            slistx.append(speedx)
            slisty.append(speedy)
        else:
            if (np.absolute(ts[i-1] - ts[i]) > 0):
                speedx = (x[i-1] - x[i])/(ts[i-1] - ts[i])
                speedy = (y[i-1] - y[i])/(ts[i-1] - ts[i])
            slistx.append(speedx)
            slisty.append(speedy)

        #if (i == 0 or i == 1):
            #slistx.append(0)
            #slisty.append(0)
        #else:
            #speedx = 0.5 * ((x[i-2] - x[i-1])/(f[i-2] - f[i-1]) + (x[i-1] - x[i])/(f[i-1] - f[i]))
            #speedy = 0.5 * ((y[i-2] - y[i-1])/(f[i-2] - f[i-1]) + (y[i-1] - y[i])/(f[i-1] - f[i]))
            #slistx.append(speedx)
            #slisty.append(speedy)
    tracklist.extend([tid for k in range(endP - startP)])
    classlist.extend([c[0] for k in range(endP - startP)])
    intersectionlist.extend([intersection[0] for k in range(endP - startP)])
    cluslist.extend([clus for k in range(endP - startP)])
    anomlist.extend([anom for k in range(endP - startP)])
    phlist.extend([ph for k in range(endP - startP)])
    sphlist.extend([sph for k in range(endP - startP)])
    ephlist.extend([eph for k in range(endP - startP)])
    plist.extend([ph for k in range(endP - startP)])
    cyclelist.extend([cycle for k in range(endP - startP)])

ndf = pd.DataFrame()

ndf['frame_id'] = framelist
ndf['track_id'] = tracklist
ndf['center_x'] = nx
ndf['center_y'] = ny
ndf['class'] = classlist
ndf['timestamp'] = timelist
ndf['intersection_id'] = intersectionlist
ndf['SPAT'] = phlist
ndf['Cycle'] = cyclelist
ndf['speed_x'] = slistx
ndf['speed_y'] = slisty
ndf['Cluster'] = cluslist
ndf['isAnomalous'] = anomlist

#ndf['InstantPhase'] = plist
#ndf['startPhase'] = sphlist
#ndf['endPhase'] = ephlist
ndf['redJump'] = rjlist
ndf['nearmiss'] = nmisslist
ndf['signature'] = siglist

#print (ndf)

ndf.to_csv(options['display'], index=False, float_format='%g')
