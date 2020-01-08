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
from sklearn.neighbors import KDTree
from rdp import rdp
from shapely.geometry import Point, Polygon
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from ast import literal_eval

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

directions = {}
directions['NS'] = literal_eval(options['NS'])
directions['SN'] = literal_eval(options['SN'])
directions['EW'] = literal_eval(options['EW'])
directions['WE'] = literal_eval(options['WE'])
phase2 = options['phase2']
dir26 = 'NS'
dir62 = 'SN'
dir48 = 'EW'
dir84 = 'WE'
if (phase2 == 'WB'):
    dir26 = 'EW'
    dir62 = 'WE'
    dir48 = 'SN'
    dir84 = 'NS'
dir_2_x_26 = directions[dir26][0][0]  #point North in NS
dir_2_y_26 = directions[dir26][0][1]  #point North in NS
dir_6_x_26 = directions[dir26][1][0]  #point Sount in NS
dir_6_y_26 = directions[dir26][1][1]  #point Sount in NS
dir_2_x_62 = directions[dir62][1][0]  #point North in SN
dir_2_y_62 = directions[dir62][1][1]  #point North in SN
dir_6_x_62 = directions[dir62][0][0]  #point Sount in SN
dir_6_y_62 = directions[dir62][0][1]  #point Sount in SN
dir_4_x_48 = directions[dir48][0][0]  #point East in EW
dir_4_y_48 = directions[dir48][0][1]  #point East in EW
dir_8_x_48 = directions[dir48][1][0]  #point West in EW
dir_8_y_48 = directions[dir48][1][1]  #point West in EW
dir_8_x_84 = directions[dir84][0][0]  #point East in WE
dir_8_y_84 = directions[dir84][0][1]  #point East in WE
dir_4_x_84 = directions[dir84][1][0]  #point East in WE
dir_4_y_84 = directions[dir84][1][1]  #point East in WE

polylist = [literal_eval(options['polygon'])]
flat_list = [item for sublist in polylist for item in sublist]
polygon = Polygon(flat_list)

#polygon = Polygon([(263.784, 307.027), (432.432,659.459), (834.595,676.757), (1042.16,376.216), (955.676,192.432), (769.73,41.0811), (598.919,41.0811), (376.216,134.054), (263.784,307.027)])
df = pd.read_csv(options['trackinfo'])
spatdf=pd.read_csv(options['spat'],names=['timestamp', 'hexphase', 'binphase', 'redphase', 'cyclenumber'], index_col='timestamp')

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)

def checkDirection(Axi, Ayi, Bxi, Byi):
    [nAxi, nAyi] = 1/np.sqrt(Axi*Axi + Ayi*Ayi) * np.array([Axi, Ayi])
    [nBxi, nByi] = 1/np.sqrt(Bxi*Bxi + Byi*Byi) * np.array([Bxi, Byi])

    return nAxi*nBxi + nAyi*nByi

def polygon_area(x,y):
    correction = x[-1] * y[0] - y[-1]* x[0]
    main_area = np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:])
    return 0.5*np.abs(main_area + correction)

def getCosineOfAngle(Axi, Ayi, Bxi, Byi):
    a = 0
    if ((Axi != 0 or Ayi != 0) and (Bxi != 0 or Byi != 0)):
        [nAxi, nAyi] = 1/np.sqrt(Axi*Axi + Ayi*Ayi) * np.array([Axi, Ayi])
        [nBxi, nByi] = 1/np.sqrt(Bxi*Bxi + Byi*Byi) * np.array([Bxi, Byi])
        a = nAxi*nBxi + nAyi*nByi
    return a

def getPathLength(tid, x, y):
    length = 0
    if (x.size > 0):
        aa = np.asarray([[x, y] for x, y in zip(x, y)])
        a = aa[:x.size-1]
        b = np.copy(aa[1:])
        d = np.sqrt(np.einsum('ij,ij->i', (a-b), (a-b))) #https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
        length = np.sum(d)
    return length

def getShapeArea(rbit, x, y):
    area = 0
    if (rbit == 2 or rbit == 6 or rbit == 4 or rbit == 8):
        area = polygon_area(x, y)
    return area

def getGreenBit(x, y):
    bit = 0
    bit2 = 0
    throughcos = float(options['through'])
    curvedcos = float(options['curved'])
    a = []
    [Ax, Ay] = [x[x.size-1]-x[0], y[y.size-1]-y[0]]
    #Check for 2-> or 6->2 alignment
    a.append(getCosineOfAngle(Ax, Ay, dir_6_x_26 - dir_2_x_26, dir_6_y_26 - dir_2_y_26))
    if (a[0] > throughcos):
        bit = 2
    else:
        a.append(getCosineOfAngle(Ax, Ay, dir_2_x_62 - dir_6_x_62, dir_2_y_62 - dir_6_y_62))
        if (a[1] > throughcos):
            bit = 6
        else:
            a.append(getCosineOfAngle(Ax, Ay, dir_8_x_48 - dir_4_x_48, dir_8_y_48 - dir_4_y_48))
            if (a[2] > throughcos):
                bit = 4
            else:
                a.append(getCosineOfAngle(Ax, Ay, dir_4_x_84 - dir_8_x_84, dir_4_y_84 - dir_8_y_84))
                if (a[3] > throughcos):
                    bit = 8
                else: # check curved paths
                    #get the best match:
                    a.append(getCosineOfAngle(Ax, Ay, dir_4_x_84-dir_2_x_26, dir_4_y_84-dir_2_y_26))
                    a.append(getCosineOfAngle(Ax, Ay, dir_8_x_48-dir_6_x_62, dir_8_y_48-dir_6_y_62))
                    a.append(getCosineOfAngle(Ax, Ay, dir_6_x_26-dir_4_x_48, dir_6_y_26-dir_4_y_48))
                    a.append(getCosineOfAngle(Ax, Ay, dir_2_x_62-dir_8_x_84, dir_2_y_62-dir_8_y_84))
                    maxa = max(a[4:])
                    m = a.index(maxa)
                    if (maxa > curvedcos):
                        if (m==4):
                            #dist = getDistanceOfPointFromLine(dir_2_x_26, dir_2_y_26, dir_6_x_26, dir_6_y_26, x[0], y[0])
                            #if (np.absolute(y[0]-dir_2_y_26) < 50): #np.absolute needed to differentiate between left turn and the far right turn
                            if (a[3] > 0.3): #np.absolute needed to differentiate between left turn and the far right turn
                                bit = 5
                                bit2 = 2
                            else:
                                bit = 0
                                bit2 = 8
                        elif (m==5):
                            #dist = getDistanceOfPointFromLine(dir_6_x_62, dir_6_y_62, dir_2_x_62,  dir_2_y_62, x[0], y[0])
                            if (a[2] > 0.3):
                                bit = 1
                                bit2 = 6
                            else:
                                bit = 0
                                bit2 = 4
                        elif (m==6):
                            #dist = getDistanceOfPointFromLine(dir_4_x_48,dir_4_y_48,dir_8_x_48,dir_8_y_48,x[0],y[0])
                            if (a[0] > 0.3):
                                bit = 7
                                bit2 = 4
                            else:
                                bit = 0
                                bit2 = 2
                        else:
                            #dist = getDistanceOfPointFromLine(dir_8_x_84,dir_8_y_84,dir_4_x_84,dir_4_y_84,x[0],y[0])
                            if (a[1] > 0.3):
                                bit = 3
                                bit2 = 8
                            else:
                                bit = 0
                                bit2 = 6
    return bit, bit2


def getVector(rbit, rbit2):
    Ax = 0
    Ay = 0
    if (rbit == 1 or rbit == 6 or rbit2 == 6):
        Ax = dir_2_x_62 - dir_6_x_62
        Ay = dir_2_y_62 - dir_6_y_62
    elif (rbit == 5 or rbit == 2 or rbit2 == 2):
        Ax = dir_6_x_26 - dir_2_x_26
        Ay = dir_6_y_26 - dir_2_y_26
    elif (rbit == 3 or rbit == 8 or rbit2 == 8):
        Ax = dir_4_x_84 - dir_8_x_84
        Ay = dir_4_y_84 - dir_8_y_84
    elif (rbit == 7 or rbit == 4 or rbit2 == 4):
        Ax = dir_8_x_48 - dir_4_x_48
        Ay = dir_8_y_48 - dir_4_y_48

    return Ax, Ay

def getStart(x, y, rbit, rbit2, mask2, j):
    i = 0
    pt = [0, 0, 0]
    Ax, Ay = getVector(rbit, rbit2)
    throughcos = float(options['through'])
    while (i < 1 and j < x.size - 1):
        a = getCosineOfAngle(x[j+1]-x[j], y[j+1]-y[j], Ax, Ay)
        if (a > throughcos):
            pt[i] = j
            pt[i+1] = j+1
            i = i+1
            j = j+1
        else:
            mask2[pt[i]] = False
            j = j+1
    while (i < 2 and j < x.size - 1):
        pt[i+1] = j+1
        a = getCosineOfAngle(x[pt[i+1]]-x[pt[i]], y[pt[i+1]]-y[pt[i]], Ax, Ay)
        if (a > throughcos):
            i = i+1
        else:
            mask2[pt[i+1]] = False
        j = j+1
    if (pt[2] == j):
        j = j + 1
    return pt, j

def getThreshold(rbit, rbit2):
    threshold = 0.9
    if (rbit == 2 or rbit == 6 or rbit == 4 or rbit == 8):
        threshold = float(options['through'])
    return threshold

def smoothenTrack(x, y, tid):
    mask2 = [True for i in range(x.size)]
    rbit, rbit2 = getGreenBit(x, y)
    pathLength = getPathLength(tid, x, y)
    wsize = 3
    i = 0
    w = [0, 0, 0]
    pt, j = getStart(x, y, rbit, rbit2, mask2, 0)

    w[0] = pt[1]
    w[1] = pt[2]
    w[2] = j

    threshold = getThreshold(rbit, rbit2)
    while (j < x.size-1):
        a = getCosineOfAngle(x[w[2]] - x[w[1]], y[w[2]] - y[w[1]], x[w[1]] - x[w[0]], y[w[1]] - y[w[0]])
        j = j+1
        if (a > threshold):
            w[0] = w[1]
            w[1] = w[2]
            w[2] = j
        else:
            w[2] = j
            mask2[j-1] = False

    #Add the last point if its not too off
    if (mask2[x.size - 1] == False):
        a = getCosineOfAngle(x[x.size - 1] - x[pt[1]], y[x.size - 1] - y[pt[1]], x[pt[1]] - x[pt[0]], y[pt[1]] - y[pt[0]])
        if (a > 0.9):
            mask2[x.size - 1] = True
    return mask2

def smoothenTrackActual(x, y, tid):
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
    return


def isShapeWeird(tid, x, y):
    status = False
    [Ax, Ay] = [x[x.size-1]-x[0], y[y.size-1]-y[0]]
    Ax = x[1] - x[0]
    Ay = y[1] - y[0]
    snum = 0
    chdir = 0
    #two breaks each of more than 50 pixels in length
    startx = x[0]
    starty = y[0]
    for i in range (1, x.size-1):
        Bx = x[i+1] - x[i]
        By = y[i+1] - y[i]
        a = getCosineOfAngle(Ax, Ay, Bx, By)
        if (a < 0.35):
            status = True
            break
        Ax = Bx
        Ay = By

    return status

def isAnomalousShape(x, y, tid):
    print ('Processing track', tid, 'for anomaly')
    model = LinearRegression(fit_intercept=True, normalize=False, copy_X=True, n_jobs=None)
    transformer = PolynomialFeatures(degree=1, include_bias=False)
    pts = x.size-2
    status = False
    if (pts == 0):
        #status = True
        return status
    xran = x[x.size-1] - x[0]
    if (np.absolute(xran) < 20):
        yran = y[y.size - 1] - y[0]
        if (np.absolute(yran) > 20):
            step = int(yran/4)
            yrangearray = np.arange(y[0], y[y.size-1], step)
            index = []
            for value in yrangearray:
                idx = (np.abs(y - value)).argmin()
                index.append(idx)
            xselected = y[index]
            yselected = x[index]
            transformer.fit(xselected.reshape(-1,1))
            xsel_ = transformer.transform(xselected.reshape(-1,1))
            model.fit(xsel_, yselected)
            x_ = transformer.transform(y.reshape(-1,1))
            y_pred = model.predict(x_)
            diffy = np.absolute(np.subtract(x, y_pred))
            threshold = 50
            mask2 = np.where(diffy > threshold, False, True)
            if (np.any(mask2 == False)):
                if (isShapeWeird(tid, xselected, yselected)):
                    status = True
    else:
        step = int(xran/4)
        xrangearray = np.arange(x[0], x[x.size-1], step)
        index = []
        for value in xrangearray:
            idx = (np.abs(x - value)).argmin()
            index.append(idx)

        xselected = x[index]
        yselected = y[index]
        transformer.fit(xselected.reshape(-1,1))
        xsel_ = transformer.transform(xselected.reshape(-1,1))
        model.fit(xsel_, yselected)
        x_ = transformer.transform(x.reshape(-1,1))
        y_pred = model.predict(x_)
        diffy = np.absolute(np.subtract(y, y_pred))
        threshold = 50 #width of a lane
        mask2 = np.where(diffy > threshold, False, True)
        if (np.any(mask2 == False)):
            if (isShapeWeird(tid, xselected, yselected)):
                status = True

    return status

def getRelevantCoordinates(pol, x, y, isPed, tid):
    s = len(x)
    excl = 0
    p = Point(x[0], y[0])
    i = 1
    while (i<s and pol.contains(p) == False):
        p = Point(x[i], y[i])
        i = i + 1

    if (i == s):
        excl = 1
        return excl, s, 0

    si = i - 1
    p = Point(x[s-1], y[s-1])
    i = s-2
    while (pol.contains(p) == False):
        p = Point(x[i], y[i])
        i = i - 1
    se = i + 1

    if (se <= si):
        excl = 1
        return excl, s, 0

    return excl, si, se

ndf = pd.DataFrame()
ldf = pd.DataFrame()

cluslist = []
alist = []

tracks = df['track_id'].unique()

nnx=[]
nny=[]
tracklist=[]
slistx=[]
slisty=[]
framelist=[]
classlist=[]
timelist=[]
dtimelist=[]
intersectionlist=[]
llist=[]
wlist=[]
hlist=[]
lenlist=[]
scalelist=[]

plist = []
rlist = []
cyclelist = []
maxpathlength = float(options['pathlength'])

p=0
for tid in tracks:
    nndf = df.loc[df['track_id'] == tid]
    tlen = 0
    nx = nndf.iloc[:,3].values
    ny = nndf.iloc[:,4].values
    f = nndf.iloc[:,0].values
    fid = nndf.iloc[:,1].values
    c = nndf.iloc[:,5].values
    ti = nndf.iloc[:,6].values
    intersection = nndf.iloc[:,7].values
    w = nndf.iloc[:,8].values
    h = nndf.iloc[:,9].values
    ny = np.ma.compressed(np.ma.masked_where(nx==0, ny))
    f = np.ma.compressed(np.ma.masked_where(nx==0, f))
    fid = np.ma.compressed(np.ma.masked_where(nx==0, fid))
    c = np.ma.compressed(np.ma.masked_where(nx==0, c))
    ti = np.ma.compressed(np.ma.masked_where(nx==0, ti))
    intersection = np.ma.compressed(np.ma.masked_where(nx==0, intersection))
    w = np.ma.compressed(np.ma.masked_where(nx==0, w))
    h = np.ma.compressed(np.ma.masked_where(nx==0, h))
    nx = np.ma.compressed(np.ma.masked_where(nx==0, nx))
    if (nx.size < 3):
        continue
    diffx = nx[nx.size-1]-nx[0]
    diffy = ny[ny.size-1]-ny[0]
    if (diffx > 10):
        nx = np.maximum.accumulate(nx)
    elif (diffx < -10):
        nx = np.minimum.accumulate(nx)
    if (diffy > 10):
        ny = np.maximum.accumulate(ny)
    elif (diffy < -10):
        ny = np.minimum.accumulate(ny)

    nsx = []
    nsy = []
    length = []
    num = len(nx)
    for i in range (0, num):
        if (i == 0):
            nsx.append(0.0)
            nsy.append(0.0)
            length.append(0.0)
        else:
            speedx = 0
            speedy = 0
            if (np.absolute(f[i-1] - f[i]) > 0):
                speedx = (nx[i-1] - nx[i])/(f[i-1] - f[i])
                speedy = (ny[i-1] - ny[i])/(f[i-1] - f[i])
            tlen = math.sqrt((nx[i-1]-nx[i])*(nx[i-1]-nx[i]) + (ny[i-1]-ny[i])*(ny[i-1]-ny[i]))
            nsx.append(speedx)
            nsy.append(speedy)
            length.append(tlen)
    if (c[0] != 'p'):
        excl, si, se = getRelevantCoordinates(polygon, nx, ny, 0, tid)
        if (excl):
            continue
        x = nx[si:se+1]
        y = ny[si:se+1]
        sx = nsx[si:se+1]
        sy = nsy[si:se+1]
        length = length[si:se+1]
        f = f[si:se+1]
        fid = fid[si:se+1]
        c = c[si:se+1]
        ti = ti[si:se+1]
        intersection = intersection[si:se+1]
        w = w[si:se+1]
        h = h[si:se+1]
        if (x.size < 5):
            continue
        if (isAnomalousShape(x, y, tid)):
            print ("Anomalous", tid)
            #continue
        pathLength = getPathLength(tid, x, y)
        if (pathLength < maxpathlength):
            continue

    X = np.array([[x, y] for x, y in zip(x, y)])
    #mask = rdp(pairs)
    e = 1
    #if (c[0] == 'pedestrian'):
        #e = 5
    mask = rdp(X, epsilon=e, return_mask=True)
    #if (c[0] != 'pedestrian'):
        #smoothen_kinks(x, y, mask)
    num = len(x)
    previndex = 0
    currindex = 0
    scale = 1
    k = 0
    lx=[]
    ly=[]
    ftime=[]
    mask2 = smoothenTrack(x, y, tid)
    for i in range (0, num):
        if (mask[i] == False):
            continue
        if (mask2[i] == False):
            continue
        if (x[i] == 0.0 and y[i] == 0.0):
            continue
        tindex = ti[i]#.replace('-24', '-04')
        phase = spatdf.loc[tindex]
        plist.append(phase['hexphase'])
        rlist.append(phase['redphase'])
        cyclelist.append(phase['cyclenumber'])
        dtimelist.append(f[i])
        nnx.append(x[i])
        nny.append(y[i])
        lx.append(x[i])
        ly.append(y[i])
        ftime.append(f[i])
        timelist.append(ti[i])
        framelist.append(fid[i])
        if (k == 0):
            slistx.append(sx[i])
            slisty.append(sy[i])
            llist.append(length[i])
        else:
            speedx = sx[i]
            speedy = sy[i]
            if (np.absolute(ftime[k-1] - ftime[k]) > 0):
                speedx = (lx[k-1] - lx[k])/(ftime[k-1] - ftime[k])
                speedy = (ly[k-1] - ly[k])/(ftime[k-1] - ftime[k])
            l = math.sqrt((lx[k-1]-lx[k])*(lx[k-1]-lx[k]) + (ly[k-1]-ly[k])*(ly[k-1]-ly[k])) * scale
            slistx.append(speedx)
            slisty.append(speedy)
            llist.append(l)
        scalelist.append(scale)
        wlist.append(w[i] * scale)
        hlist.append(h[i] * scale)
        tracklist.append(tid)
        classlist.append(c[0])
        intersectionlist.append(intersection[0])
        k = k + 1

ndf['time'] = dtimelist
ndf['frame_id'] = framelist
ndf['track_id'] = tracklist
ndf['center_x'] = nnx
ndf['center_y'] = nny
ndf['class'] = classlist
ndf['timestamp'] = timelist
ndf['intersection_id'] = intersectionlist
ndf['speed_x'] = slistx
ndf['speed_y'] = slisty
ndf['dist'] = llist
ndf['SPAT'] = plist
ndf['Cycle'] = cyclelist


ndf.to_csv(options['combo'], index=False)
