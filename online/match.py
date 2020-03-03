import numpy as np
import pandas as pd
import math
import datetime
import time
import pymysql
import sys
import multiprocessing as mp
import scipy as sp
from scipy import signal
from ast import literal_eval
from hash import HashTable
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
from shapely.geometry import Point, Polygon
from sqlalchemy import create_engine
from rdp import rdp

if (len(sys.argv) == 1):
        print ('Error: provide config file')

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

def movingaverage(values, window):
    weights = np.repeat(1.0, window)/window
    smas = np.convolve(values, weights, 'valid')
    return smas

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
# NS ((1122.16, 456.216), (125.405, 451.872))
# SN ((125.405, 451.872), (1122.16, 456.216))
# EW ((624.865, 644.324), (674.595, 276.757))
# WE ((540.541, 136.216), (557.838, 856.216))
phase2stopbar = literal_eval(options['phase2stopbar'])
phase4stopbar = literal_eval(options['phase4stopbar'])
phase6stopbar = literal_eval(options['phase6stopbar'])
phase8stopbar = literal_eval(options['phase8stopbar'])
polylist = [literal_eval(options['polygon'])]
flat_list = [item for sublist in polylist for item in sublist]
polygon = Polygon(flat_list)
pedpolygon = polygon

finalCentroids = []
finalCentroidName = []
objects=options['objects'].split(',')

#@dataclass
class Track:
    def __init__(self, tid, uid, x, y, origx, origy, t, ts, f, iid, classtype, rrange, nm, sig):
        self.tid = tid        #track id
        self.rrange = rrange
        self.x = x
        self.arr = y
        self.t = t
        self.ts = ts
        self.f = f
        self.sx = np.zeros(shape=x.size)
        self.sy = np.zeros(shape=x.size)
        self.origx = origx
        self.origy = origy
        self.cid = -1
        self.phasetype = '' #major, minor, leftmajor, leftminor
        self.classtype = classtype #other types are bus, truck, ped
        self.cycles = []      # cycle numbers for which the vehicle is present in the intersection
        self.slope = 0
        self.isAnomalous = 0
        self.isrj = 0
        self.start = 0
        self.end = 0
        self.scale = 1
        self.len = 0
        self.iid = str(iid)
        self.clus = 0
        self.iscent = 0
        self.isthrough = True
        self.gbit = 0
        self.gbit2 = 0
        self.phase = 0
        self.uid = uid
        self.issmall = 0
        self.isoutside = 0
        self.hispeed = 0
        self.oddshape = 0
        self.nm = nm
        self.sig = sig

class PhaseCycle:
    def __init__(self, sTS, eTS, rTS, r, g, m, isG, isY, isR, ped):
        self.eTimestamp = eTS  #end timestamp
        self.rTimestamp = rTS  #start of red timestamp
        self.aor = r           #number of arrivals on red for cycle
        self.aog = g           #number of arrivals on green for cycle
        self.isMaxed = m       #1 of the phase maxed out or was forced off
        self.isGreen = isG     #1 if phase is green, 0 otherwise
        self.isYellow = isY     #1 if phase is green, 0 otherwise
        self.isRed = isR     #1 if phase is green, 0 otherwise
        self.isPedGreen = isG
        self.isPedYellow = isY
        self.isPedRed = isR
        self.pedstatus = ped   #1 if ped call recorded, 0 otherwise
        self.skipDetector = 0
        self.finish = 0
        self.sIndex = 0
        self.eIndex = 0
        self.cumaor = 0
        self.cumaog = 0

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
np.seterr(all='raise')
tracksH = HashTable(50000)
debugH = HashTable(100000)

def swap (x1, x2):
    tmp = x1
    x1 = x2
    x2 = tmp
    return x1, x2

#vlines = [[(0, 0), (0, 960)], [(332.973, 646.486), (315.676, 315.676)], [(436.757, 261.622), (490.811, 741.622)], [(518.919, 263.784), (577.297, 745.946)], [(605.405, 261.622), (640.743, 743.784)], [(737.297, 274.595), (750.27, 683.243)], [(880.339, 459), (882.162, 611.892)], [(1280, 0), (1280, 960)]]
#hlines = [[(322.162, 0), (877.838, 0)], [(322.162, 352.432), (877.838, 352.432)], [(322.162, 406.486), (910.27, 406.486)], [(332.541, 454.054), (966.486, 454.054)], [(322.162, 516.757), (877.838, 508.108)], [(321.162, 562.162), (877.838, 560)], [(322.162, 607.568), (877.838, 611.892)], [(402.162, 735.135), (715.676, 735.135)], [(332.160, 960), (332.160, 960)]]
#https://math.stackexchange.com/questions/274712/calculate-on-which-side-of-a-straight-line-is-a-given-point-located
vlinetuples=literal_eval(options['verticalgrid'])
hlinetuples=literal_eval(options['horizontalgrid'])
vlines=[list(ele) for ele in vlinetuples]
hlines=[list(ele) for ele in hlinetuples]

lenhlines = len(hlines)
lenvlines = len(vlines)
for i in range(0,lenhlines,2):
    x1 = hlines[i][0]
    y1 = hlines[i][1]
    x2 = hlines[i+1][0]
    y2 = hlines[i+1][1]
    if (x1 > x2):
        #must be arranged in increasing order of x
        tmpx = hlines[i][0]
        tmpy = hlines[i][1]
        hlines[i][0] = hlines[i+1][0]
        hlines[i][1] = hlines[i+1][1]
        hlines[i+1][0] = tmpx
        hlines[i+1][1] = tmpy

for i in range(0,lenvlines,2):
    x1 = vlines[i][0]
    y1 = vlines[i][1]
    x2 = vlines[i+1][0]
    y2 = vlines[i+1][1]
    if (y1 > y2):
        #must be arranged in increasing order of y
        tmpx = vlines[i][0]
        tmpy = vlines[i][1]
        vlines[i][0] = vlines[i+1][0]
        vlines[i][1] = vlines[i+1][1]
        vlines[i+1][0] = tmpx
        vlines[i+1][1] = tmpy

lowerpvlines = vlines[0::2]
upperpvlines = vlines[1::2]
leftphlines = hlines[0::2]
rightphlines = hlines[1::2]
pminusupperlower = np.subtract(upperpvlines,lowerpvlines)
pminusrightleft = np.subtract(rightphlines, leftphlines)
pminusupperlower[:,[0,1]] = pminusupperlower[:,[1,0]]
pminusrightleft[:,[0,1]] = pminusrightleft[:,[1,0]]

def locatePointOnGridOld(x, y):
    lenhlines = len(hlines)
    lenvlines = len(vlines)
    d1 = -1
    for i in range(0,lenhlines,2):
        x1 = hlines[i][0]
        y1 = hlines[i][1]
        x2 = hlines[i+1][0]
        y2 = hlines[i+1][1]
        #must be arranged in increasing order of x
        if (x1 > x2):
            x1, x2 = swap (x1, x2)
            y1, y2 = swap (y1, y2)
        d = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
        sameSign = (d * d1) > 0
        if (sameSign == False):
            break
        d1 = d

    d1 = -1
    for j in range(0, lenvlines,2):
        x1 = vlines[j][0]
        y1 = vlines[j][1]
        x2 = vlines[j+1][0]
        y2 = vlines[j+1][1]
        #must be arranged in increasing order of y
        if (y1 > y2):
            x1, x2 = swap (x1, x2)
            y1, y2 = swap (y1, y2)
        d = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
        sameSign = (d * d1) > 0
        if (sameSign == False):
            break
        d1 = d
    return j, i #vertical, horizontal

def locatePointOnGrid(x, y):
    i=0
    j=0
    length = int(lenvlines/2)
    pxy = np.array([(x,y)] * length)
    pminuspointlower = np.subtract(pxy,lowerpvlines)
    multres = np.multiply(pminuspointlower, pminusupperlower)
    d = np.subtract(multres[:,0], multres[:,1])
    dsign = np.sign(d)
    signchange = ((np.roll(dsign, 1) - dsign) != 0).astype(int)
    wheresignchange = np.where(signchange==1)
    if (wheresignchange[0].size):
        i = wheresignchange[0][-1]

    length = int(lenhlines/2)
    pxy = np.array([(x,y)] * length)
    pminuspointleft = np.subtract(pxy,leftphlines)
    multres = np.multiply(pminuspointleft, pminusrightleft)
    d = np.subtract(multres[:,0], multres[:,1])
    dsign = np.sign(d)
    signchange = ((np.roll(dsign, 1) - dsign) != 0).astype(int)
    wheresignchange = np.where(signchange==1)
    if (wheresignchange[0].size):
        j = wheresignchange[0][-1]

    return i, j #vertical, horizontal

def getEuclideanDistance(Ax, Ay, Bx, By):
    return math.sqrt((Ax-Bx)**2 + (Ay-By)**2)

def isCurved(nt):
    ntsize = nt.x.size-1
    a = 0 #orthogonal
    status = True

    if (ntsize < 2):
        return a, status

    #midpoint = int(ntsize/2)
    [Axi, Ayi] = [nt.x[ntsize] - nt.x[0], nt.arr[ntsize] - nt.arr[0]]
    dirval = directions.values()
    for i in dirval:
        [Bxi, Byi] = [i[1][0] - i[0][0], i[1][1] - i[0][1]]
        a = np.absolute(getCosineOfAngle(Axi, Ayi, Bxi, Byi))
        if (a > 0.95):
            status = False
            break

    return status

def distanceOld(nt, nt1, boxDim):
    tdist = 100
    distance = 100000
    [Ax, Ay] = [nt.x[0], nt.arr[0]]
    [Bx, By] = [nt1.x[0], nt1.arr[0]]
    startdim = getEuclideanDistance(Ax, Ay, Bx, By)
    if (startdim < boxDim):
        [Ax, Ay] = [nt.x[nt.x.size-1], nt.arr[nt.x.size-1]]
        [Bx, By] = [nt1.x[nt1.x.size-1], nt1.arr[nt1.x.size-1]]
        enddim = getEuclideanDistance(Ax, Ay, Bx, By)
        if (enddim < boxDim):
            ntpair = [[x, y] for x, y in zip(nt.x, nt.arr)]
            nt1pair = [[x, y] for x, y in zip(nt1.x, nt1.arr)]
            distance, path = fastdtw (ntpair, nt1pair, dist=euclidean)
    lenratio = 0
    maxdiff = 100
    angled=180
    if (distance < 35000):
        first = 1
        prevpair = (0,0)
        longpath = 0
        area = 0
        z1=0
        z=0
        first = 1
        distlist = []
        el = 0
        el1 = 0
        for pair in path:
            if (first):
                first = 0
            elif ((pair [0] != prevpair[0]) and (pair[1] != prevpair[1])):
                [Ax, Ay] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Bx, By] = [nt1.x[pair[1]], nt1.arr[pair[1]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea1 = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen1 = math.sqrt((nt1.x[pair[1]] - nt1.x[prevpair[1]])**2 + (nt1.arr[pair[1]]-nt1.arr[prevpair[1]])**2)
                area = area + secarea1
                z1 = z1 + seclen1
                [Ax, Ay] = [nt.x[prevpair[0]], nt.arr[prevpair[0]]]
                [Bx, By] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen = math.sqrt((nt.x[prevpair[0]]-nt.x[pair[0]])**2 + (nt.arr[prevpair[0]]-nt.arr[pair[0]])**2)
                area = area + secarea
                z = z + seclen
                if ((seclen + seclen1) > 0):
                    distlist.append((secarea + secarea1)/(seclen + seclen1))
            elif ((pair[0] != prevpair[0])):
                [Ax, Ay] = [nt.x[prevpair[0]], nt.arr[prevpair[0]]]
                [Bx, By] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen = math.sqrt((nt.x[prevpair[0]]-nt.x[pair[0]])**2 + (nt.arr[prevpair[0]]-nt.arr[pair[0]])**2)
                area = area + secarea
                z = z + seclen
                if (seclen > 0):
                    distlist.append(secarea/seclen)
            elif ((pair[1] != prevpair[1])):
                [Ax, Ay] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Bx, By] = [nt1.x[pair[1]], nt1.arr[pair[1]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen = math.sqrt((nt1.x[pair[1]] - nt1.x[prevpair[1]])**2 + (nt1.arr[pair[1]]-nt1.arr[prevpair[1]])**2)
                area = area + secarea
                z1 = z1 + seclen
                if (seclen > 0):
                    distlist.append(secarea/seclen)
            prevpair = pair

        distarray = np.asarray(distlist)
        distarray = distarray - np.mean(distarray)
        maxdiff = np.amax(np.absolute(distarray))
        if ((z + z1) > 0):
            tdist = area/(z + z1)
            lenratio = (z + z1)/(z + z1 + el + el1)
        longlength = z
        if (longpath == 1):
            longlength = z1
    return tdist+maxdiff, maxdiff, distance, lenratio

def distanceR1(nt, nt1, boxDim):
    tdist = 100
    distance = 100000
    prevpair = (0,0)
    area = 0
    z1=0
    z=0
    first = 1
    ntpair = [[x, y] for x, y in zip(nt.x, nt.arr)]
    nt1pair = [[x, y] for x, y in zip(nt1.x, nt1.arr)]
    distance, path = fastdtw (ntpair, nt1pair, dist=euclidean)
    if (distance < 35000):
        for pair in path:
            if (first):
                first = 0
            elif ((pair [0] != prevpair[0]) and (pair[1] != prevpair[1])):
                [Ax, Ay] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Bx, By] = [nt1.x[pair[1]], nt1.arr[pair[1]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea1 = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen1 = math.sqrt((nt1.x[pair[1]] - nt1.x[prevpair[1]])**2 + (nt1.arr[pair[1]]-nt1.arr[prevpair[1]])**2)
                area = area + secarea1
                z1 = z1 + seclen1
                [Ax, Ay] = [nt.x[prevpair[0]], nt.arr[prevpair[0]]]
                [Bx, By] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen = math.sqrt((nt.x[prevpair[0]]-nt.x[pair[0]])**2 + (nt.arr[prevpair[0]]-nt.arr[pair[0]])**2)
                area = area + secarea
                z = z + seclen
            elif ((pair[0] != prevpair[0])):
                [Ax, Ay] = [nt.x[prevpair[0]], nt.arr[prevpair[0]]]
                [Bx, By] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen = math.sqrt((nt.x[prevpair[0]]-nt.x[pair[0]])**2 + (nt.arr[prevpair[0]]-nt.arr[pair[0]])**2)
                area = area + secarea
                z = z + seclen
            elif ((pair[1] != prevpair[1])):
                [Ax, Ay] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Bx, By] = [nt1.x[pair[1]], nt1.arr[pair[1]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen = math.sqrt((nt1.x[pair[1]] - nt1.x[prevpair[1]])**2 + (nt1.arr[pair[1]]-nt1.arr[prevpair[1]])**2)
                area = area + secarea
                z1 = z1 + seclen
            prevpair = pair

    if ((z + z1) > 0):
        tdist = area/(z + z1)
    longlength = z
    return tdist, distance

def distance(nt, nt1, boxDim):
    tdist = 100
    distance = 100000
    ntpair = [[x, y] for x, y in zip(nt.x, nt.arr)]
    nt1pair = [[x, y] for x, y in zip(nt1.x, nt1.arr)]
    distance, path = fastdtw (ntpair, nt1pair, dist=euclidean)
    lenratio = 0
    maxdiff = 100
    angled=180
    firstdiff = 0
    lastdiff = 0
    if (distance > 0):
    #if (distance < 55000):
        first = 1
        prevpair = (0,0)
        longpath = 0
        area = 0
        z1=0
        z=0
        first = 1
        distlist = []
        alist = []
        secalist = []
        secalist1 = []
        zlist = []
        z1list = []
        el = 0
        el1 = 0
        for pair in path:
            if (first):
                first = 0
            elif ((pair [0] != prevpair[0]) and (pair[1] != prevpair[1])):
                [Ax, Ay] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Bx, By] = [nt1.x[pair[1]], nt1.arr[pair[1]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea1 = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen1 = math.sqrt((nt1.x[pair[1]] - nt1.x[prevpair[1]])**2 + (nt1.arr[pair[1]]-nt1.arr[prevpair[1]])**2)
                area = area + secarea1
                z1 = z1 + seclen1
                [Ax, Ay] = [nt.x[prevpair[0]], nt.arr[prevpair[0]]]
                [Bx, By] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                seclen = math.sqrt((nt.x[prevpair[0]]-nt.x[pair[0]])**2 + (nt.arr[prevpair[0]]-nt.arr[pair[0]])**2)
                area = area + secarea
                z = z + seclen
                alist.append(secarea1+secarea)
                zlist.append(seclen)
                z1list.append(seclen1)
                secalist.append(secarea)
                secalist1.append(secarea1)
                if ((seclen + seclen1) > 0):
                    distlist.append((secarea + secarea1)/(seclen + seclen1))
#            elif ((pair[0] != prevpair[0])):
#                [Ax, Ay] = [nt.x[prevpair[0]], nt.arr[prevpair[0]]]
#                [Bx, By] = [nt.x[pair[0]], nt.arr[pair[0]]]
#                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
#                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
#                seclen = math.sqrt((nt.x[prevpair[0]]-nt.x[pair[0]])**2 + (nt.arr[prevpair[0]]-nt.arr[pair[0]])**2)
#                area = area + secarea
#                z = z + seclen
#                if (seclen > 0):
#                    distlist.append(secarea/seclen)
#            elif ((pair[1] != prevpair[1])):
#                [Ax, Ay] = [nt.x[pair[0]], nt.arr[pair[0]]]
#                [Bx, By] = [nt1.x[pair[1]], nt1.arr[pair[1]]]
#                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
#                secarea = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
#                seclen = math.sqrt((nt1.x[pair[1]] - nt1.x[prevpair[1]])**2 + (nt1.arr[pair[1]]-nt1.arr[prevpair[1]])**2)
#                area = area + secarea
#                z1 = z1 + seclen
#                if (seclen > 0):
#                    distlist.append(secarea/seclen)

            prevpair = pair

        distarray = np.asarray(distlist)
        distarray = distarray - np.mean(distarray)
        maxdiff = np.amax(np.absolute(distarray))
        zlistlen = len(zlist)
        n = zlistlen-2
        if (zlistlen > 1):
            if (zlist[1] > 0 and z1list[1] > 0):
                firstdiff = max(secalist[1]/zlist[1], secalist1[1]/z1list[1])
            elif (zlist[1] + z1list[1] > 0):
                firstdiff = (alist[1]/(zlist[1] + z1list[1]))
        if (n > 0):
            if (zlist[n] > 0 and z1list[n] > 0):
                lastdiff = max(secalist[n]/zlist[n], secalist1[n]/z1list[n])
            elif (zlist[n] + z1list[n] > 0):
                lastdiff = alist[n]/(zlist[n] + z1list[n])
        if ((z + z1) > 0):
            tdist = area/(z + z1)
            lenratio = (z + z1)/(z + z1 + el + el1)
        longlength = z
        if (longpath == 1):
            longlength = z1

    return tdist, maxdiff, distance, lenratio, firstdiff, lastdiff

def getMxmVal(x, y, mxm, ovl):
    if (x == y):
        return 0.0
    l = (int(x), int(y))
    return (mxm.get(l, -1), ovl.get(l, -1))

def computeCentroid(trackseries, centroids, currentIndex, mxm, ovl, boxDim):
    mindist = 10000
    mintrack = -1
    maxlen = 0
    tracklist = []
    lengthlist = []
    slen = trackseries.size
    if (slen == 0):
        return centroids[currentIndex]

    num = max(1, int(0.5*slen))
    #Compute centroid among 10% of members chosen randomly
    if (num == 1):
        return tracksH[trackseries.iloc[0]]

    rseries = trackseries.sample(n=num, random_state=1)
    if (currentIndex >= 0):
        cseries = pd.Series(centroids[currentIndex].tid)
        rseries = rseries.append(cseries)
    for i, track in rseries.iteritems():
        nt1 = tracksH[track]
        dist = 0
        for j, otrack in trackseries.iteritems():
            nt2 = tracksH[otrack]
            if (nt1.tid != nt2.tid):
                val, ovlp = getMxmVal(nt1.tid, nt2.tid, mxm, ovl)
                if (val != -1):
                    dist = dist + val
                else:
                    val, maxdiff, dtwdist, ovlp = distance(nt1, nt2, boxDim)
                    mxm[nt1.tid, nt2.tid] = val
                    mxm[nt2.tid, nt1.tid] = val
                    ovl[nt1.tid, nt2.tid] = ovlp
                    ovl[nt2.tid, nt1.tid] = ovlp
                    dist = dist + val
        if (dist < mindist):
            mindist = dist
            mintrack = track
    return tracksH[mintrack]

def getCosineOfAngle(Axi, Ayi, Bxi, Byi):
    a = 0
    if ((Axi != 0 or Ayi != 0) and (Bxi != 0 or Byi != 0)):
        [nAxi, nAyi] = 1/np.sqrt(Axi*Axi + Ayi*Ayi) * np.array([Axi, Ayi])
        [nBxi, nByi] = 1/np.sqrt(Bxi*Bxi + Byi*Byi) * np.array([Bxi, Byi])
        a = nAxi*nBxi + nAyi*nByi
    return a

def computeDistance(x1, y1, x2, y2):
    return np.sqrt(np.add(np.square(x1-x2), np.square(y1-y2)))

def isLengthDifferent(nt1, nt2):
    len1 = nt1.len
    len2 = nt2.len

    status = False
    biggertrack = nt1
    smallertrack = nt2
    if (len1 < 0.75 * len2 or len2 < 0.75 * len1):
        status = True
        if (len1 < len2):
            biggertrack = nt2
            smallertrack = nt1
    else:
        starttostart = computeDistance(nt1.x[0], nt1.arr[0], nt2.x[0], nt2.arr[0])
        endtoend = computeDistance(nt1.x[nt1.x.size-1], nt1.arr[nt1.x.size-1], nt2.x[nt2.x.size-1], nt2.arr[nt2.x.size-1])
        if (starttostart > 0.3 * len1 or starttostart > 0.3 * len2 or endtoend > 0.3 * len1 or endtoend > 0.3 * len2):
            status = True
            if (len1 < len2):
                biggertrack = nt2
                smallertrack = nt1

    return biggertrack, smallertrack, status

def isCollinear(nt1, nt2):
    #https://www.gamedev.net/forums/topic/556821-check-if-vectors-are-parallel-and-pointing-in-the-same-direction-with-tolerance/
    nt1size = nt1.x.size-1
    nt2size = nt2.x.size-1

    #create vectors
    [Axi, Ayi] = [nt1.x[0] - nt1.x[nt1size], nt1.arr[0] - nt1.arr[nt1size]]
    [Bxi, Byi] = [nt1.x[0] - nt2.x[0], nt1.arr[0] - nt2.arr[0]]
    [Cxi, Cyi] = [nt1.x[0] - nt2.x[nt2size], nt1.arr[0] - nt2.arr[nt2size]]

    biggerTrack, smallertrack, status = isLengthDifferent(nt1, nt2)

    if (status):
        bsize = biggerTrack.x.size - 1
        ssize = smallertrack.x.size - 1
        dist1 = computeDistance(biggerTrack.x[0], biggerTrack.arr[0], smallertrack.x[0], smallertrack.arr[0])
        dist2 = computeDistance(biggerTrack.x[0], biggerTrack.arr[0], smallertrack.x[ssize], smallertrack.arr[ssize])
        dist3 = computeDistance(biggerTrack.x[bsize], biggerTrack.arr[bsize], smallertrack.x[0], smallertrack.arr[0])
        dist4 = computeDistance(biggerTrack.x[bsize], biggerTrack.arr[bsize], smallertrack.x[ssize], smallertrack.arr[ssize])
        fpdist = dist1 + dist2
        epdist = dist3 + dist4
        [Axi, Ayi] = [biggerTrack.x[0] - biggerTrack.x[bsize], biggerTrack.arr[0] - biggerTrack.arr[bsize]]
        if (fpdist < epdist):
            [Bxi, Byi] = [biggerTrack.x[bsize] - smallertrack.x[0], biggerTrack.arr[bsize] - smallertrack.arr[0]]
            [Cxi, Cyi] = [biggerTrack.x[bsize] - smallertrack.x[ssize], biggerTrack.arr[bsize] - smallertrack.arr[ssize]]
        else:
            [Bxi, Byi] = [biggerTrack.x[0] - smallertrack.x[0], biggerTrack.arr[0] - smallertrack.arr[0]]
            [Cxi, Cyi] = [biggerTrack.x[0] - smallertrack.x[ssize], biggerTrack.arr[0] - smallertrack.arr[ssize]]
    #normalize

    a1 = np.absolute(getCosineOfAngle(Axi, Ayi, Bxi, Byi))
    a2 = np.absolute(getCosineOfAngle(Axi, Ayi, Cxi, Cyi))

    status = False
    if (a1 > 0.99 and a2 > 0.99):
        status = True

    return status

def checkDirection(nt1, nt2):
    #https://www.gamedev.net/forums/topic/556821-check-if-vectors-are-parallel-and-pointing-in-the-same-direction-with-tolerance/
    nt1size = nt1.x.size-1
    nt2size = nt2.x.size-1
    #create vectors
    [Axi, Ayi] = [nt1.x[0] - nt1.x[nt1size], nt1.arr[0] - nt1.arr[nt1size]]
    [Bxi, Byi] = [nt2.x[0] - nt2.x[nt2size], nt2.arr[0] - nt2.arr[nt2size]]

    return getCosineOfAngle(Axi, Ayi, Bxi, Byi)

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

def strictly_decreasing(L):
    return all(x>y for x, y in zip(L, L[1:]))

def findSegmentEndPoint(x, y, x0, y0):
    idx = (np.abs(x - x0)).argmin()
    if (strictly_increasing(x)):
        if (x[idx] > x0):
            idx = idx + 1
    elif (strictly_decreasing(x)):
        if (x[idx] < x0):
            idx = idx + 1
    else: #non-monotonic
        idx = (np.abs(y - y0)).argmin()
        if (strictly_increasing(y)):
            if (y[idx] > y0):
                idx = idx + 1
        elif (strictly_decreasing(y)):
            if (y[idx] < y0):
                idx = idx + 1

    if (idx == x.size):
        idx = idx - 2
    elif (idx == x.size-1):
        idx = idx - 1

    return idx

def isOnSameLane(nt1, nt2):
    dist = 100
    #compare with respect to the bigger track
    big = nt1
    small = nt2
    if (nt1.len < nt2.len):
        big = nt2
        small = nt1
    x0 = small.x[0]
    y0 = small.arr[0]
    pts = small.x.size
    ptsidx = 1
    if (pts > 4):
        ptsidx = int(pts/2)
    xn = small.x[ptsidx]
    yn = small.arr[ptsidx]
    xlast = small.x[small.x.size-1]
    ylast = small.arr[small.x.size-1]
    n = big.x.size
    status = False

    if (big.isthrough == False):
        distlist = []
        xdiff = np.subtract(big.x, x0)
        ydiff = np.subtract(big.arr, y0)
        distdiff = np.sqrt(np.add(np.square(xdiff), np.square(ydiff)))
        mindiststart = np.min(distdiff)
        xdiff2 = np.subtract(big.x, xlast)
        ydiff2 = np.subtract(big.arr, ylast)
        distdiff2 = np.sqrt(np.add(np.square(xdiff2), np.square(ydiff2)))
        mindistlast = np.min(distdiff2)
        if (mindiststart < 32):
            idx = distdiff.argmin()
            if (idx == 0):
                x1 = big.x[0]
                y1 = big.arr[0]
                x2 = big.x[1]
                y2 = big.arr[1]
            else:
                x1 = big.x[idx-1]
                y1 = big.arr[idx-1]
                x2 = big.x[idx]
                y2 = big.arr[idx]
            a = getCosineOfAngle(x1-x2, y1-y2, x0-xn, y0-yn)
            if (a > 0.85):
                #dist = np.absolute((y2-y1) * x0 - (x2 - x1) * y0 + x2*y1 - y2*x1)/(np.sqrt(np.add(np.square(y2-y1), np.square(x2-x1))))
                dist = getDistanceOfPointFromLine(x1, y1, x2, y2, x0, y0)
                status = dist < 32
                #print ('isOnSameLane', mindiststart, mindistlast, dist)
            else:
                if (idx < big.x.size-2):
                    x2 = big.x[idx+1]
                    y2 = big.arr[idx+1]
                    a = getCosineOfAngle(x1-x2, y1-y2, x0-xn, y0-yn)
                    if (a > 0.85):
                        #dist = np.absolute((y2-y1) * x0 - (x2 - x1) * y0 + x2*y1 - y2*x1)/(np.sqrt(np.add(np.square(y2-y1), np.square(x2-x1))))
                        dist = getDistanceOfPointFromLine(x1, y1, x2, y2, x0, y0)
                        status = dist < 32
                        #print ('isOnSameLane', mindiststart, mindistlast, dist)

    else:
        x1 = big.x[0]
        y1 = big.arr[0]
        x2 = big.x[n-1]
        y2 = big.arr[n-1]

        dist1 = np.absolute((y2-y1) * x0 - (x2 - x1) * y0 + x2*y1 - y2*x1)/(np.sqrt(np.add(np.square(y2-y1), np.square(x2-x1))))
        m = (y2 - y1)/(x2 - x1)
        x = (m * (m*x1 - y1) + (x0 + m * y0))/(1 + m*m)
        y = (m * (x0 + m * y0) - (m*x1-y1))/(1 + m*m)
        dist = np.sqrt(np.add(np.square(x0 - x), np.square(y0 - y)))
        status = dist < 20
    return status

#point is (x0, y0), line is given by two points (x1, y1) and (x2, y2)
def getDistanceOfPointFromLine(x1, y1, x2, y2, x0, y0):
    dist = np.absolute((y2-y1) * x0 - (x2 - x1) * y0 + x2*y1 - y2*x1)/(np.sqrt(np.add(np.square(y2-y1), np.square(x2-x1))))
    return dist

def printClusterCenter(centroids, cid, rrange, cid_list, x_list, y_list):
    center = centroids[cid]
    y1 = center.arr
    for i in range (0, rrange):
        if (y1[i] > 0):
            cid_list.append(cid)
            x_list.append(i)
            y_list.append(y1[i])

def printClustersToFile(cdf, filename, mode):
    df_tracks_with_clusters = pd.DataFrame()
    df_tracks_with_clusters['track_id'] = cdf['track_id'].astype(int)
    df_tracks_with_clusters['closest'] = cdf['closest']
    df_tracks_with_clusters['isAnomalous'] = cdf['isAnomalous']
    df_tracks_with_clusters['startPoint'] = cdf['startPoint']
    df_tracks_with_clusters['endPoint'] = cdf['endPoint']
    df_tracks_with_clusters['phase'] = cdf['phase']
    df_tracks_with_clusters['startphase'] = cdf['startphase']
    df_tracks_with_clusters['endphase'] = cdf['endphase']

    df_tracks_with_clusters.dropna(inplace=True)
    df_tracks_with_clusters.to_csv(filename,  float_format='%.f', mode=mode, header=False, index=False)

def getTrackMatchingThresholdForMembers(nt1, nt2):
    #compare lengths of nt1 and nt2
    threshold = 1.0
    tdistthresh = 25
    c1 = nt1.isthrough
    c2 = nt2.isthrough
    if (c1 == False and c2 == False):
        threshold = 0.98
    elif (c1 == True and c2 == True):
        biggerTrack, smallertrack, status = isLengthDifferent(nt1, nt2)
        lb = biggerTrack.len
        ls = smallertrack.len
        if (biggerTrack.isthrough and lb > 250 and ls > 0):
            threshold = 0.70
            tdistthresh = 60 * lb/ls
        else:
            threshold = 0.75 #getSweep(biggerTrack)
            tdistthresh = 60
    else:
        biggerTrack, smallertrack, status = isLengthDifferent(nt1, nt2)
        lb = biggerTrack.len
        ls = smallertrack.len
        if (biggerTrack.isthrough and lb > 250 and ls > 0):
            threshold = 0.70
            tdistthresh = 60 * lb/ls
        else:
            threshold = 0.97
            tdistthresh = 50
    return threshold, tdistthresh

def closeToMembers(centi, centj, isPed):
    #find the longest track
    status = True
    direction = checkDirection(centi, centj)
    scale = 1
    tdist = 100

    if (isPed):
        if (direction < 0.82):
            tdist, dtwdist = distance(centi, centj, 1000)
            tdist = tdist * scale
            if (tdist > 75 or dtwdist > 12000):
                status = False
            #print(centi.tid, centj.tid, tdist, maxdiff, dtwdist, ovlp, direction, status)
    else:
        threshold, tdistthresh = getTrackMatchingThresholdForMembers(centi, centj)
        if (direction < threshold):
            status = False
        else:
            isSameLane = True #isOnSameLane(centi, centj)
            if (isSameLane == False):
                status = False
            else:
                tdist, dtwdist = distance(centi, centj, 1000)
                tdist = tdist * scale
                #if (tdist > 15 or maxdiff > 200 or (dtwdist > 10000 and maxdiff < 20)):  # second half is an indication that the path diverges out
                    #status = False
                if (tdist > tdistthresh):
                    if (tdist < 50 and isCollinear(centi, centj)): 
                        nt1 = centi
                        nt2 = centj
                        n1 = nt1.x.size
                        n2 = nt2.x.size
                        status = True
                        #if (tdist < 25):
                            #status = True
                        #elif (nt1.x[n1-1] < nt2.x[0] or nt2.x[n2-1] < nt1.x[0] and strictly_increasing(nt1.x)):
                            #status = True
                        #elif (nt1.x[n1-1] > nt2.x[0] or nt2.x[n2-1] > nt1.x[0] and strictly_decreasing(nt1.x)):
                            #status = True
                        #else:
                            #status = False
                    #elif (tdist < 20 and maxdiff < 25):
                        #status = True
                    elif (tdist < 45 and centj.isthrough):
                        status = True
                    else:
                        status = False
    return status, tdist

def getPathLength(tid, x, y):
    length = 0
    if (x.size > 0):
        aa = np.asarray([[x, y] for x, y in zip(x, y)])
        a = aa[:x.size-1]
        b = np.copy(aa[1:])
        d = np.sqrt(np.einsum('ij,ij->i', (a-b), (a-b))) #https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
        length = np.sum(d)
    return length

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
                            bit = 0
                            bit2 = 8
                            dist = min(np.sqrt(np.add(np.square(polygon.centroid.x - x), np.square(polygon.centroid.y - y))))
                            if (dist < 130): #np.absolute needed to differentiate between left turn and the far right turn
                                bit = 5
                                bit2 = 2
                            else:
                                dist = getDistanceOfPointFromLine(phase6stopbar[0][0],phase6stopbar[0][1],phase6stopbar[1][0],phase6stopbar[1][1],x[0],y[0])
                                if (dist > 300):
                                    bit = 5
                                    bit2 = 2
                                #[Bx, By] = [x[1] - x[0], y[1] - y[0]]
                                #[Ex, Ey] = [x[x.size-1]-x[x.size-2], y[y.size-1]-y[y.size-2]]
                                #begangle = getCosineOfAngle(Bx, By, dir_4_x_84 - dir_8_x_84, dir_4_y_84 - dir_8_y_84)
                                #endangle = getCosineOfAngle(Ex, Ey, dir_4_x_84 - dir_8_x_84, dir_4_y_84 - dir_8_y_84)
                                #if (begangle < 0.1 or endangle > 0.75):
                                    #bit = 5
                                    #bit2 = 2
                        elif (m==5):
                            #dist = getDistanceOfPointFromLine(dir_6_x_62, dir_6_y_62, dir_2_x_62,  dir_2_y_62, x[0], y[0])
                            bit = 0
                            bit2 = 4
                            dist = min(np.sqrt(np.add(np.square(polygon.centroid.x - x), np.square(polygon.centroid.y - y))))
                            if (dist < 130):
                                bit = 1
                                bit2 = 6
                            else:
                                dist = getDistanceOfPointFromLine(phase2stopbar[0][0],phase2stopbar[0][1],phase2stopbar[1][0],phase2stopbar[1][1],x[0],y[0])
                                if (dist > 300):
                                    bit = 1
                                    bit2 = 6
                        elif (m==6):
                            #dist = getDistanceOfPointFromLine(dir_4_x_48,dir_4_y_48,dir_8_x_48,dir_8_y_48,x[0],y[0])
                            bit = 0
                            bit2 = 2
                            dist = min(np.sqrt(np.add(np.square(polygon.centroid.x - x), np.square(polygon.centroid.y - y))))
                            if (dist < 130):
                                bit = 7
                                bit2 = 4
                            else:
                                dist = getDistanceOfPointFromLine(phase8stopbar[0][0],phase8stopbar[0][1],phase8stopbar[1][0],phase8stopbar[1][1],x[0],y[0])
                                if (dist > 300):
                                    bit = 7
                                    bit2 = 4
                        else:
                            #dist = getDistanceOfPointFromLine(dir_8_x_84,dir_8_y_84,dir_4_x_84,dir_4_y_84,x[0],y[0])
                            bit = 0
                            bit2 = 6
                            dist = min(np.sqrt(np.add(np.square(polygon.centroid.x - x), np.square(polygon.centroid.y - y))))
                            if (dist < 130):
                                bit = 3
                                bit2 = 8
                            else:
                                dist = getDistanceOfPointFromLine(phase4stopbar[0][0],phase4stopbar[0][1],phase4stopbar[1][0],phase4stopbar[1][1],x[0],y[0])
                                if (dist > 300):
                                    bit = 3
                                    bit2 = 8
    return bit, bit2

def convertToPixelsPerSec(spkmph):
    return 50 * spkmph / (3 * 2.23694)

def setSpeed(x, y, ts, nt, cent):
    slistx = []
    slisty = []
    index = 0
    check = False
    slistx.append(0)
    slisty.append(0)
    if (x.size > 1):
        aa = [(x, y) for x, y in zip(x, y)]
        a = aa[:x.size-1]
        b = np.copy(aa[1:])
        t1 = ts[:x.size-1]
        t2 = ts[1:]
        dsub = np.subtract(b, a)
        tsub = np.subtract(t2, t1)
        t_2d = np.asarray([[tx, ty] for tx, ty in zip(tsub, tsub)])
        sp = np.divide(dsub, t_2d)
        sp2 = np.sqrt(np.einsum('ij,ij->i', sp, sp))
        speedlimit = float(options['speedlimit'])
        speedlimitwithbuffer = speedlimit + 15
        maxspeedlimit = speedlimit*2
        splpixels = convertToPixelsPerSec(speedlimitwithbuffer)
        maxpixels = convertToPixelsPerSec(maxspeedlimit)
        if (np.amax(sp2) > maxpixels):
            check = True
        else:
            highsp = np.where(sp2>splpixels)
            spthresh = convertToPixelsPerSec(20)
            highspentries = sp[highsp]
            if (nt.isAnomalous == 0 and np.logical_and(abs(highspentries[:, 0]) > spthresh, abs(highspentries[:, 1]) > spthresh).any()):
                spq = highspentries[np.logical_and(abs(highspentries[:, 0]) > spthresh , abs(highspentries[:, 1]) > spthresh)]
                Ax = cent.x[cent.x.size-1] - cent.x[0]
                Ay = cent.arr[cent.x.size-1] - cent.arr[0]
                Bx = spq[0][0] * 0.1 #speed * time = distance
                By = spq[0][1] * 0.1
                a = np.absolute(getCosineOfAngle(Ax, Ay, Bx, By))
                throughthresh = float(options['through'])
                curvedthresh = float(options['curved'])
                if (cent.isthrough and a < throughthresh or cent.isthrough == False and a < curvedthresh):
                    check = True
        index = np.argmax(sp2 > 100)
        sx = sp[:, [0]].reshape(x.size-1,)
        sy = sp[:, [1]].reshape(x.size-1,)
        slistx.extend(np.round(sx, 2).tolist())
        slisty.extend(np.round(sy, 2).tolist())
    return slistx, slisty, index, check

isBaseSet = 0
i = 0
cycle_no = 0
cycle_mark = -1
j = -1
timeIndex = 1
eventCIndex = 2
eventPIndex = 3

step = datetime.timedelta(seconds=1)
sigstatus = bytearray(b'000000000000000011111111')
sigPedstatus = bytearray(b'000000000000000011111111')
previousTime = -1
phaseInfo = []
numPhases = 16
currentPhase = ''
currentPedPhase = ''
k = 0
while k <= numPhases:      #assign and initialize a phase cycle list for each phase
    cycleList = []
    cycleList.insert(0, PhaseCycle(0, 0, 0, 0, 0, 0, 0, 0, 1, 0)) #Wont the first location have the last phase information?
    phaseInfo.append(cycleList)
    k = k + 1

def processEventList(eventList, currentTime, reportTime, spatrows):
    global isBaseSet
    global i
    global cycle_no
    global cycle_mark
    global j
    global currentPhase
    global currentPedPhase
    if (i != 0 and currentTime > reportTime): #time to dump data, one second is over
        while (reportTime < currentTime):
            k = 1
            while k <= 8: #Iterate over phases to gather AOR/AOG/MaxOut/ForceOff
                #Remember: an intersection cycle may have multiple phase cycles
                #Sum up values from all cycles of phase k
                #The last cycle in the list is the oldest one, and the first is the newest
                if (len(phaseInfo[k]) > 0):
                    sigstatus[k-1] = ord(str(phaseInfo[k][0].isGreen))
                    sigstatus[k+7] = ord(str(phaseInfo[k][0].isYellow))
                    sigstatus[k+15] = ord(str(phaseInfo[k][0].isRed))
                    sigPedstatus[k-1] = ord(str(phaseInfo[k][0].isPedGreen))
                    sigPedstatus[k+7] = ord(str(phaseInfo[k][0].isPedYellow))
                    sigPedstatus[k+15] = ord(str(phaseInfo[k][0].isPedRed))
                    #print (k, sigstatus.decode())
                k = k+1     #Next phase
            if ('-' not in sigstatus.decode()):
                hexstr = "%x" % int(sigstatus.decode(), 2)
                l = len(hexstr)
                paddedzeroes = 6-l
                fhexstr = '0' * paddedzeroes + hexstr
                redstr = hexstr[0:1]
                #print (fhexstr)
                j = j + 1
                #convert reportTime
                #roundedtime = reportTime.round('1s')
                phasestr = '0' * paddedzeroes + hexstr
                #if (debugH[str(reportTime.strftime("%Y-%m-%d %H:%M:%S.%f"))] != None):
                    #debug = 1
                if (currentPhase != phasestr):
                    spatrow = (options['cityintid'], options['cameraid'], reportTime.strftime("%Y-%m-%d %H:%M:%S.%f"), phasestr, str(cycle_no), 'Gainesville','FL',1)
                    #debugH[str(reportTime.strftime("%Y-%m-%d %H:%M:%S.%f"))] = 1
                    spatrows.append(spatrow)
                    currentPhase = phasestr
            if ('-' not in sigPedstatus.decode()):
                hexstr = "%x" % int(sigPedstatus.decode(), 2)
                l = len(hexstr)
                paddedzeroes = 6-l
                fhexstr = '0' * paddedzeroes + hexstr
                redstr = hexstr[0:1]
                #print (fhexstr)
                j = j + 1
                #convert reportTime
                #roundedtime = reportTime.round('1s')
                phasestr = '0' * paddedzeroes + hexstr
                #if (debugH[str(reportTime.strftime("%Y-%m-%d %H:%M:%S.%f"))] != None):
                    #debug = 1
                if (currentPedPhase != phasestr):
                    spatrow = (options['cityintid'], options['cameraid'], reportTime.strftime("%Y-%m-%d %H:%M:%S.%f"), phasestr, str(cycle_no), 'Gainesville','FL',0)
                    #debugH[str(reportTime.strftime("%Y-%m-%d %H:%M:%S.%f"))] = 1
                    spatrows.append(spatrow)
                    currentPedPhase = phasestr
            #reportTime = reportTime + step
            reportTime = currentTime
    for row in eventList:
        if (row[eventCIndex] == 1):     #Phase begin green
            if (isBaseSet == 0):
                baseTime = currentTime
                isBaseSet = 1
            if (i == 0):       #The first phase turns green in the CSV
                i = i+1
                cycle_mark = row[eventPIndex]
                #reportTime = currentTime.round('1s') + step
                reportTime = currentTime
            if (i > 0):
                currentCycle = phaseInfo[row[eventPIndex]][0]
                currentCycle.isGreen = 1            #Set the isGreen flag for the phase
                currentCycle.isYellow = 0            #Set the isGreen flag for the phase
                currentCycle.isRed = 0            #Set the isGreen flag for the phase
                if (row[eventPIndex] == cycle_mark):
                    cycle_no = cycle_no + 1
        elif (i==0):           #Have not yet encountered the first green in the CSV, so continue to next CSV row
            continue
        elif (row[eventCIndex] == 9):                   #Phase end yellow clearence
            #length = len(phaseInfo[row[2]])
            currentCycle = phaseInfo[row[eventPIndex]][0]
            currentCycle.isGreen = 0
            currentCycle.isYellow = 0
            currentCycle.isRed = 1
        elif (row[eventCIndex] == 8):                   # Phase begin yellow clearence
            currentCycle = phaseInfo[row[eventPIndex]][0]
            if (currentCycle.isGreen == 1):
                currentCycle.isGreen = 0
                currentCycle.isYellow = 1
                currentCycle.isRed = 0
        elif (row[eventCIndex] == 21):                   #Ped begin walk
            #length = len(phaseInfo[row[2]])
            currentCycle = phaseInfo[row[eventPIndex]][0]
            currentCycle.isPedGreen = 1
            currentCycle.isPedYellow = 0
            currentCycle.isPedRed = 0
        elif (row[eventCIndex] == 23):                   #Ped dont walk
            #length = len(phaseInfo[row[2]])
            currentCycle = phaseInfo[row[eventPIndex]][0]
            currentCycle.isPedGreen = 0
            currentCycle.isPedYellow = 0
            currentCycle.isPedRed = 1
        elif (row[eventCIndex] == 22):                   # Ped begin clearence (countdown)
            currentCycle = phaseInfo[row[eventPIndex]][0]
            if (currentCycle.isPedGreen == 1):
                currentCycle.isPedGreen = 0
                currentCycle.isPedYellow = 1
                currentCycle.isPedRed = 0
    return reportTime

def convertToSPaT(sdf, currentTime):
    global previousTime
    spatrows = []
    #reportTime = pd.to_datetime('today')
    reportTime = -1
    eventList = []
    for row_id, row in enumerate(sdf.values):
        #print (lineno)
        currentTime = pd.to_datetime(row[timeIndex], format = '%m/%d/%y %H:%M:%S.%f')
        if (reportTime == -1):
            reportTime = currentTime
        elif (currentTime != previousTime):
            reportTime = processEventList(eventList, previousTime, reportTime, spatrows)
            eventList.clear()
            previousTime = currentTime
        if (row[eventCIndex] == 9 or row[eventCIndex] == 8 or row[eventCIndex] == 1):
            eventList.append(row)
        if (row[eventCIndex] == 21 or row[eventCIndex] == 22 or row[eventCIndex] == 23):
            eventList.append(row)
    if (reportTime != currentTime):
        reportTime = processEventList(eventList, previousTime, reportTime, spatrows)
        eventList.clear()
        previousTime = currentTime
    return spatrows

def getSPaTIndex(spatdf, ts):
    index = np.searchsorted(spatdf.timestamp.values, ts) - 1
    return index

def adjustTrackBeginning(begIndex, track):
    track.x = track.x[begIndex:]
    track.arr = track.arr[begIndex:]
    track.t = track.t[begIndex:]
    track.ts = track.ts[begIndex:]
    track.f = track.f[begIndex:]
    track.sx = track.sx[begIndex:]
    track.sy = track.sy[begIndex:]
    track.cycles = track.cycles[begIndex:]

def adjustTrackEnding(endIndex, track):
    track.x = track.x[0:endIndex]
    track.arr = track.arr[0:endIndex]
    track.t = track.t[0:endIndex]
    track.ts = track.ts[0:endIndex]
    track.f = track.f[0:endIndex]
    track.sx = track.sx[0:endIndex]
    track.sy = track.sy[0:endIndex]
    track.cycles = track.cycles[0:endIndex]

def trackStartsAndEndsBeforeStopBar(pol, x, y, isPed):
    status = False
    if (isPed):
        return status

    s = len(x)
    p = Point(x[0], y[0])
    i = 1
    while (i<s and pol.contains(p) == False):
        p = Point(x[i], y[i])
        i = i + 1

    if (i == s):
        status = True

    return status

def processTrack(tid, videodf, spatdf, geoclusters, centroid_list, centroids, tH):
    updaterows = []
    insertrows = []
    inserttrackrows = []
    updateindex = -1
    oldclus = 'anom'
    s1 = time.time()
    ndf = videodf.loc[videodf['track_id'] == tid]
    t = ndf['time']
    f = ndf['frame_id']
    origx = ndf['center_x'].copy()
    origy = ndf['center_y'].copy()
    x = ndf['center_x'].replace(to_replace=0, method='ffill').replace(to_replace=0, method='bfill')
    y = ndf['center_y'].replace(to_replace=0, method='ffill').replace(to_replace=0, method='bfill')
    c = ndf['class']
    ts = ndf['timestamp']
    iid = ndf['intersection_id']
    nm = ndf['nearmiss']
    sig = ndf['signature']
    uid = ndf['unique_ID']
    debug = 1
    s2 = time.time()
    x.reset_index(inplace=True, drop=True)
    y.reset_index(inplace=True, drop=True)
    if (x.size > 0 and x[0] == 0 and y[0] == 0):
        print ("Ignoring 0 coordinate track or track section")
        return updaterows, '', inserttrackrows, '', 1
        #return createDummyTrack(tH, tid, updateindex, uid.iloc[0], x, y, origx, origy, 1, 0, c.iloc[0], iid.iloc[0], t, ts, f, 0, 0, 1260, updaterows, inserttrackrows, nm[0], sig[0], 'small', 1)
    if (x.size > 2):
        x = sp.signal.medfilt(x, 3)
        y = sp.signal.medfilt(y, 3)
    else:
        x = np.asarray(x)
        y = np.asarray(y)
    t = np.asarray(t)
    f = np.asarray(f)
    ts = np.asarray(ts)
    nm = np.asarray(nm)
    sig = np.asarray(sig)
    origx = np.asarray(origx)
    origy = np.asarray(origy)
    s5 = time.time()
    nt = tH[tid]
    rrange = 1260
    if (nt != None and nt.iscent != 1):
        updateindex = nt.x.size
        nt.x = np.concatenate([nt.x, x])
        nt.arr = np.concatenate([nt.arr, y])
        nt.t = np.concatenate([nt.t, t])
        nt.ts = np.concatenate([nt.ts, ts])
        nt.f = np.concatenate([nt.f, f])
        nt.origx = origx
        nt.origy = origy
        nt.sx = np.zeros(shape=nt.x.size)
        nt.sy = np.zeros(shape=nt.x.size)
        oldclus = nt.clus
        #have to make new assessments for the following
        nt.issmall = 0
        nt.isoutside = 0
        nt.oddshape = 0
        nt.hispeed = 0
        nt.isAnomalous = 0
    else:
        nt = Track(tid, uid.iloc[0], x, y, origx, origy, t, ts, f, iid.iloc[0], c.iloc[0], rrange, nm[0], sig[0])
    s6 = time.time()

    isPed = False
    if (c.iloc[0]=='pedestrian'):
        isPed = True
    #check for tracks with very small lengrh with respect to the centroid
    if (nt.x.size <= 2 ):
        phaseindex = 0
    else:
        phaseindex = ts.size - 2
    if (isPed == False):
        vehspat = spatdf[spatdf['type']==1]
        exactmatch = vehspat.iloc[np.searchsorted(vehspat.timestamp.values, ts[phaseindex]+(f[phaseindex]%10)*np.timedelta64(100,'ms')) - 1]
    else:
        pedspat = spatdf[spatdf['type']==0]
        if (pedspat.empty):
            return createTrackData(nt, updateindex, updaterows, inserttrackrows, 'pednotallowed', 0)
        exactmatch = pedspat.iloc[np.searchsorted(pedspat.timestamp.values, ts[phaseindex]+(f[phaseindex]%10)*np.timedelta64(100,'ms')) - 1]
    #exactmatch = getSPaTIndex(spatdf, ts[phaseindex]+(f[phaseindex]%10)*np.timedelta64(100,'ms'))
    rs = exactmatch['hexphase']
    cycle = exactmatch['cycle']
    nt.phasetype = rs
    nt.cycles.append(cycle)
    if (nt.x.size <= 3 ):
        nt.issmall = 1
        return createTrackData(nt, updateindex, updaterows, inserttrackrows, 'small', 0)
    if (trackStartsAndEndsBeforeStopBar(polygon, nt.x, nt.arr, isPed)):
        nt.isoutside = 1
        return createTrackData(nt, updateindex, updaterows, inserttrackrows, 'outside', 0)
    nt.len = getPathLength(tid, x, y)
    cent = findMatchingCentroid(nt, geoclusters, centroid_list, centroids)
    sx, sy, spindex, check = setSpeed(nt.x, nt.arr, nt.t, nt, cent)
    nt.sx = np.asarray(sx)
    nt.sy = np.asarray(sy)
    if (check):
        nt.hispeed = 1
        #return createTrackData(nt, updateindex, updaterows, inserttrackrows, 'hispeed', 0)
    if (nt.oddshape == 1):
        nt.oddshape = 1
        return createTrackData(nt, updateindex, updaterows, inserttrackrows, 'shape', 0)
    #check for the length of the track relative to the centroid; ed = euclidean diatance
    edtrack = np.sqrt(np.add(np.square(nt.x[0]-nt.x[nt.x.size-1]), np.square(nt.arr[0]-nt.arr[nt.arr.size-1])))
    edcentroid = np.sqrt(np.add(np.square(cent.x[0]-cent.x[cent.x.size-1]), np.square(cent.arr[0]-cent.arr[cent.arr.size-1])))
    threshold = 0.1
    if (edtrack/edcentroid < 0.1):
        nt.issmall = 1
        return createTrackData(nt, updateindex, updaterows, inserttrackrows, 'small', 0)

    skipbeg, skipend, skipmid, skipangular = trimTrack(nt, cent)
    s7mid = time.time()
    isrj = 0
    speedphase = rs
    s8 = time.time()

    si = 0
    se = np.size(x)
    r = np.size(x)
    phasetype=rs
    pathLength = getPathLength(tid, x, y)
    cindex = objects.index(nt.classtype)
    pedindex = objects.index('pedestrian')
    s9 = time.time()
    if (cindex != pedindex):
        nt.isrj = checkPhaseAndMovement(nt, vehspat)
    s10 = time.time()
    si = 0
    if (updateindex > -1):
        si = updateindex
        updaterow = createUpdateRow(nt, nt.clus)
        updaterows.append(updaterow)
        print ('update', str(nt.phase), nt.clus, str(nt.isAnomalous), str(nt.isrj), str(tid), str(nt.uid))
    s11 = time.time()
    if (nt.x.size != skipbeg.size):
        print (nt.tid, nt.x.size, skipbeg.size)
    outputdf = createTrackFrame(nt, si, skipbeg, skipend, skipmid, skipangular)
    #for index in range(si, ts.size-1):
        #pdts = pd.to_datetime(ts[index])
        #tstr = pdts.strftime("%Y-%m-%d %H:%M:%S") + '.' + str(f[index]%10)
        #insertrow = (str(f[index]), str(tid), str(x[index]), str(y[index]), str(c[index-si]), tstr, options['cityintid'], str(iid[index-si]), rs, str(cycle), str(sx[index]), str(sy[index]))
        #insertrows.append(insertrow)
    s12 = time.time()
    if (updateindex == -1):
        pdts = pd.to_datetime(ts[0])
        tstr = pdts.strftime("%Y-%m-%d %H:%M:%S")
        inserttrackrow = (options['cityintid'], str(iid.iloc[0]), tstr, str(nt.tid), str(nt.phase), nt.clus, str(nt.isAnomalous), str(nt.isrj), str(nm[0]), str(sig[0]), str(nt.uid), str(nt.issmall), str(nt.isoutside), str(nt.hispeed), str(nt.oddshape))
        inserttrackrows.append(inserttrackrow)
        print (options['cityintid'], str(iid.iloc[0]), tstr, str(nt.tid), str(nt.phase), cent.phase, nt.clus, str(nt.isAnomalous), str(nt.isrj), str(nm[0]), str(sig[0]), str(nt.uid), str(nt.issmall), str(nt.isoutside), str(nt.hispeed), str(nt.oddshape))
        #mysqlwb = time.time()
        #mycursor.execute(writequery, insertrow)
        #myoutputdb.commit()
    s13 = time.time()
    #print (tid, ts.size, s2-s1, s3-s2, s4-s3, s5-s4, s6-s5, s7-s6, s7mid-s7, s8-s7mid, s9-s8, s10-s9, s11-s10, s12-s11, s13-s12)
    return updaterows, outputdf, inserttrackrows, nt, 0

def createUpdateRow(nt, clus):
    pdts = pd.to_datetime(nt.ts[0])
    tstr = pdts.strftime("%Y-%m-%d %H:%M:%S")
    updaterow = (str(nt.phase), clus, str(nt.isAnomalous), str(nt.isrj), str(nt.issmall), str(nt.isoutside), str(nt.hispeed), str(nt.oddshape), tstr, str(nt.tid), str(nt.uid))
    return updaterow

def createTrackFrame(nt, si, skipbeg, skipend, skipmid, skipangular):
    outputdf = pd.DataFrame()
    outputdf['frame_id'] = nt.f[si:nt.ts.size]
    outputdf['track_id'] = [nt.tid] * (nt.ts.size-si)
    outputdf['center_x'] = nt.x[si:nt.ts.size]
    outputdf['center_y'] = nt.arr[si:nt.ts.size]
    outputdf['class'] = [nt.classtype] * (nt.ts.size-si)
    outputdf['timestamp'] = pd.to_datetime(nt.ts[si:nt.ts.size]) + nt.f[si:nt.ts.size]%10*np.timedelta64(100,'ms')
    outputdf['intersection_id'] = [options['cityintid']] * (nt.ts.size-si)
    outputdf['camera_id'] = [nt.iid] * (nt.ts.size-si)
    outputdf['SPAT'] = [nt.phasetype] * (nt.ts.size-si)
    if (len(nt.cycles) > 0):
        outputdf['cycle'] = [nt.cycles[0]] * (nt.ts.size-si)
    else:
        outputdf['cycle'] = [0] * (nt.ts.size-si)
    outputdf['speed_x'] = nt.sx[si:nt.ts.size]
    outputdf['speed_y'] = nt.sy[si:nt.ts.size]
    outputdf['unique_ID'] = [nt.uid] * (nt.ts.size-si)
    outputdf['skip_begin'] = skipbeg[si:nt.ts.size]
    outputdf['skip_end'] = skipend[si:nt.ts.size]
    outputdf['skip_mid'] = skipmid[si:nt.ts.size]
    outputdf['skip_angle'] = skipangular[si:nt.ts.size]
    outputdf['raw_x'] = nt.origx
    outputdf['raw_y'] = nt.origy
    print (outputdf)
    #outputdf = pd.DataFrame()
    #outputdf['frame_id'] = nt.f[si:nt.ts.size-1]
    #outputdf['track_id'] = [tid] * (nt.ts.size-1-si)
    #outputdf['center_x'] = nt.x[si:nt.ts.size-1]
    #outputdf['center_y'] = nt.arr[si:nt.ts.size-1]
    #outputdf['class'] = [c] * (nt.ts.size-1-si)
    #outputdf['timestamp'] = pd.to_datetime(nt.ts[si:nt.ts.size-1]) + nt.f[si:nt.ts.size-1]%10*np.timedelta64(100,'ms')
    #outputdf['intersection_id'] = [options['cityintid']] * (nt.ts.size-1-si)
    #outputdf['camera_id'] = [iid] * (nt.ts.size-1-si)
    #outputdf['SPAT'] = [rs] * (nt.ts.size-1-si)
    #outputdf['cycle'] = [cycle] * (nt.ts.size-1-si)
    #outputdf['speed_x'] = nt.sx[si:nt.ts.size-1]
    #outputdf['speed_y'] = nt.sy[si:nt.ts.size-1]
    #outputdf['unique_ID'] = [nt.uid] * (nt.ts.size-1-si)
    #outputdf['skipBegin'] = [0] * (nt.ts.size-1-si)
    #outputdf['skipEnd'] = [0] * (nt.ts.size-1-si)
    #outputdf['skipMiddle'] = [0] * (nt.ts.size-1-si)
    #outputdf['medFilter'] = [0] * (nt.ts.size-1-si)
    return outputdf

#def createDummyTrack(tH, tid, si, uid, x, y, origx, origy, issmall, isoutside, classtype, iid, t, ts, f, rs, cycle, ranget, updaterows, inserttrackrows, nm, sig, clusname, ignore):
#    skipbeg = np.zeros((x.size,), dtype=int)
#    skipend = np.zeros((x.size,), dtype=int)
#    skipmid = np.zeros((x.size,), dtype=int)
#    skipangular = np.zeros((x.size,), dtype=int)
#    if (si != -1):
#        nt = tH[tid]
#        outputdf = createTrackFrame(nt, si, skipbeg, skipend, skipmid, skipangular)
#        updaterow = createUpdateRow(nt, clusname)
#        updaterows.append(updaterow)
#        return updaterows, outputdf, inserttrackrows, nt, 0
#    nt = Track(tid, uid, x, y, origx, origy, t, ts, f, iid, classtype, ranget)
#    nt.clus = clusname
#    nt.isAnomalous = 1
#    nt.sx = [0] * x.size
#    nt.sy = [0] * x.size
#    nt.origx = origx
#    nt.origy = origy
#    nt.issmall = 1
#    nt.isoutside = 1
#    pdts = pd.to_datetime(nt.ts[0])
#    tstr = pdts.strftime("%Y-%m-%d %H:%M:%S")
#    outputdf = createTrackFrame(nt, 0, skipbeg, skipend, skipmid, skipangular)
#    inserttrackrow = (options['cityintid'], str(iid), tstr, str(nt.tid), str(nt.phase), nt.clus, str(nt.isAnomalous), str(nt.isrj), str(nm), str(sig), str(nt.uid), str(nt.issmall), str(nt.isoutside), str(nt.hispeed), str(nt.oddshape))
#    inserttrackrows.append(inserttrackrow)
#    return updaterows, outputdf, inserttrackrows, nt, ignore

def createTrackData(nt, si, updaterows, inserttrackrows, clusname, ignore):
    skipbeg = np.zeros((nt.x.size,), dtype=int)
    skipend = np.zeros((nt.x.size,), dtype=int)
    skipmid = np.zeros((nt.x.size,), dtype=int)
    skipangular = np.zeros((nt.x.size,), dtype=int)
    nt.isAnomalous = 1
    nt.clus = clusname
    if (si != -1):
        outputdf = createTrackFrame(nt, si, skipbeg, skipend, skipmid, skipangular)
        updaterow = createUpdateRow(nt, clusname)
        updaterows.append(updaterow)
        return updaterows, outputdf, inserttrackrows, nt, 0
    pdts = pd.to_datetime(nt.ts[0])
    tstr = pdts.strftime("%Y-%m-%d %H:%M:%S")
    outputdf = createTrackFrame(nt, 0, skipbeg, skipend, skipmid, skipangular)
    inserttrackrow = (options['cityintid'], str(nt.iid), tstr, str(nt.tid), str(nt.phase), nt.clus, str(nt.isAnomalous), str(nt.isrj), str(nt.nm), str(nt.sig), str(nt.uid), str(nt.issmall), str(nt.isoutside), str(nt.hispeed), str(nt.oddshape))
    inserttrackrows.append(inserttrackrow)
    return updaterows, outputdf, inserttrackrows, nt, 0

def checkPhaseAndMovement(nt, spatdf):
    x = nt.x
    y = nt.arr
    isrj = 0
    rbit, rbit2 = getGreenBit(x, y)
    nt.gbit = rbit
    nt.gbit2 = rbit2
    if (rbit == 0):
        #nt.phase = rbit2
        permittedPhases = options['permittedright']
        if str(rbit2) not in permittedPhases:
            rbit = rbit2
    if (rbit == 2 or rbit == 5):
        stopbar = phase2stopbar
    elif (rbit == 4 or rbit == 7):
        stopbar = phase4stopbar
    elif (rbit == 6 or rbit == 1):
        stopbar = phase6stopbar
    elif (rbit == 8 or rbit == 3):
        stopbar = phase8stopbar
    if (rbit != 0):
        crossts = timeTrackCrossesStopBar(nt, stopbar)
        if (crossts != 0):
            exactmatch = spatdf.iloc[getSPaTIndex(spatdf, crossts)]
            speedphase = exactmatch['hexphase']
            bin_string = bin(int('1'+speedphase, 16))[3:]
            if (bin_string[rbit+16-1] == '1'):
                if (rbit2 > 0):
                    if (bin_string[rbit2+16-1] == '1'):
                        isrj = 1
                else:
                    isrj = 1
            #if (debug == 1):
                #if (nt.tid==2359):
                    #print ('debug', rs, crossts, exactmatch, speedphase, bin_string, rbit, rbit2, bin_string[rbit2+16-1])
    return isrj

def angle(dir):
    """
    Returns the angles between vectors.

    Parameters:
    dir is a 2D-array of shape (N,M) representing N vectors in M-dimensional space.

    The return value is a 1D-array of values of shape (N-1,), with each value
    between 0 and pi.

    0 implies the vectors point in the same direction
    pi/2 implies the vectors are orthogonal
    pi implies the vectors point in opposite directions
    """
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return np.arccos((dir1*dir2).sum(axis=1)/(
        np.sqrt((dir1**2).sum(axis=1)*(dir2**2).sum(axis=1))))

def trimTrack(track, centroid):
    mask = np.ones((track.x.size,), dtype=int)
    skipbeg = np.zeros((track.x.size,), dtype=int)
    skipmid = np.zeros((track.x.size,), dtype=int)
    skipend = np.zeros((track.x.size,), dtype=int)
    skipangular = np.zeros((track.x.size,), dtype=int)
    ntpair = [[x, y] for x, y in zip(track.x, track.arr)]
    nt1pair = [[x, y] for x, y in zip(centroid.x, centroid.arr)]
    dist, path = fastdtw(ntpair, nt1pair, dist=euclidean)
    first = 1
    count = 0
    found = 0
    prevpair = path[0]
    for pair in path[1:]:
        if (pair[0] != prevpair[0] and (pair[1] == prevpair[1]) and pair[1] != centroid.x.size-1):
            mask[prevpair[0]] = 0
        #if centroid track coordinate is 0 and actual track coordinate is > 0, mark skipbeg index
        if (prevpair[1] == 0 and pair[0] > 0):
            skipbeg[prevpair[0]] = 1
        #else if centroid track coordinate is the last one for two consecutive entries, mark the skipend index
        elif (prevpair[1] == pair[1]):
            skipend[pair[0]] = 1
        elif (pair[1] == prevpair[1]):
            skipmed[prevpair[0]] = 1
        prevpair = pair

    if (track.tid==1484):
        print (1484, path)
        print (mask)
        print (track.arr)
    modx = track.x[np.where(mask)]
    mody = track.arr[np.where(mask)]
    if (track.tid==1484):
        print (track.arr)
        print (track.x)
    tolerance = 20
    origxy = [[x, y] for x, y in zip(track.x, track.arr)]
    xy = [[x, y] for x, y in zip(modx, mody)]
    simplified = np.array(rdp(xy, tolerance))
    min_angle = np.pi*0.35
    directions = np.diff(simplified, axis=0)
    theta = angle(directions)
    idx = np.where(theta>min_angle)[0]+1
    if (track.tid == 1484):
        print (simplified)
        print (min_angle, idx, theta)
    if (idx.any()):
        coord = simplified[idx[[idx.size-1]]]
        print (coord)
        midx = np.argwhere(np.all(origxy == coord, axis=1))[0]
        skipangular = np.zeros((track.x.size,), dtype=int)
        skipangular[:midx[0]] = 1
        if (track.tid==1484):
            print (mask)
            print (track.arr)
            print (track.x)
        
    return skipbeg, skipend, skipmid, skipangular

def trimTrack2(track, centroid):
    mask = np.ones((track.x.size,), dtype=int)
    skipbeg = np.zeros((track.x.size,), dtype=int)
    skipmid = np.zeros((track.x.size,), dtype=int)
    skipend = np.zeros((track.x.size,), dtype=int)
    skipangular = np.ones((track.x.size,), dtype=int)
    ntpair = [[x, y] for x, y in zip(track.x, track.arr)]
    nt1pair = [[x, y] for x, y in zip(centroid.x, centroid.arr)]
    dist, path = fastdtw(ntpair, nt1pair, dist=euclidean)
    first = 1
    count = 0
    found = 0
    for pair in path:
        if (first):
            first = 0
        elif (pair[0] != prevpair[0] and (pair[1] == prevpair[1]) and pair[1] != centroid.x.size-1):
            mask[prevpair[0]] = 0
        #if centroid track coordinate is 0 and actual track coordinate is > 0, mark skipbeg index
        if (prevpair[1] == 0 and pair[0] > 0):
            skipbeg[prevpair[0]] = 1
        #else if centroid track coordinate is the last one for two consecutive entries, mark the skipend index
        elif (prevpair[1] == pair[1]):
            skipend[pair[0]] = 1
        elif (pair[1] == prevpair[1]):
            skipmed[prevpair[0]] = 1
        prevpair = pair

    if (track.tid==1484):
        print (1484, path)
        print (mask)
        print (track.arr)
    track.x = track.x[np.where(mask)]
    track.arr = track.arr[np.where(mask)]
    track.t = track.t[np.where(mask)]
    track.ts = track.ts[np.where(mask)]
    track.f = track.f[np.where(mask)]
    #track.sx = track.sx[np.where(mask)] #commented out because speed is set later
    #track.sy = track.sy[np.where(mask)]
    track.x = sp.signal.medfilt(track.x, 3)
    track.arr = sp.signal.medfilt(track.arr, 3)
    #track.cycles = track.cycles[np.where(mask)]
    if (track.tid==1484):
        print (track.arr)
        print (track.x)
    tolerance = 20
    xy = [[x, y] for x, y in zip(track.x, track.arr)]
    simplified = np.array(rdp(xy, tolerance))
    min_angle = np.pi*0.35
    directions = np.diff(simplified, axis=0)
    theta = angle(directions)
    idx = np.where(theta>min_angle)[0]+1
    if (track.tid == 1484):
        print (simplified)
        print (min_angle, idx, theta)
    if (idx.any()):
        coord = simplified[idx[[idx.size-1]]]
        print (coord)
        midx = np.argwhere(np.all(xy == coord, axis=1))[0]
        mask = np.ones((track.x.size,), dtype=int)
        mask[:midx[0]] = 0
        track.x = track.x[np.where(mask)]
        track.arr = track.arr[np.where(mask)]
        track.t = track.t[np.where(mask)]
        track.ts = track.ts[np.where(mask)]
        track.f = track.f[np.where(mask)]
        track.sx = track.sx[np.where(mask)]
        track.sy = track.sy[np.where(mask)]
        if (track.tid==1484):
            print (mask)
            print (track.arr)
            print (track.x)
        
    return skipbeg, skipend, skipmid, skipangular

def readCentroidData(myoutputdb):
    geoClusters = [[[] for j in range(len(hlines))] for i in range(len(vlines))]
    centquery='select * from CentroidsTable where intersection_id=%(name1)s and camera_id=%(name2)s;'
    centdf = pd.read_sql(centquery, myoutputdb, params={'name1': options['cityintid'], 'name2' : options['cameraid']})
    centroid_tracks = centdf['track_id'].unique()
    for tid in centroid_tracks:
        ndf = centdf.loc[centdf['track_id'] == tid]
        x = ndf.iloc[:,1].values
        y = ndf.iloc[:,2].values
        c = ndf.iloc[:,3].values
        iid = ndf.iloc[:,4].values
        ps = ndf.iloc[:,6].values #phases
        name = ndf.iloc[:,7].values
        r = np.size(x)
        cindex = objects.index(c[0])
        r = np.size(x)
        nt = Track(tid, 0, x, y, x, y, 0, 0, 0, iid[0], c[0], 0, 0, -1)
        nt.start = 0
        nt.end = x.size
        nt.phasetype = ps[r-2]
        nt.len = getPathLength(tid, x, y)
        nt.iscent = 1
        if (c[0] != 'pedestrian'):
            #rbit, rbit2 = getGreenBit(x, y)
            #tracksH[tid] = nt
            #nt.gbit = rbit
            #nt.gbit2 = rbit2
            nt.phase = int(ps[0])
        finalCentroids.append(nt)
        finalCentroidName.append(name[0])
        if (isCurved(nt)):
            nt.isthrough = False
    allocateClustersToGrid(geoClusters)
    return geoClusters

def main():
    engine_string = 'mysql+pymysql://username:password@host/dbname'
    engine_string = engine_string.replace('username', options['user'])
    engine_string = engine_string.replace('password', options['passwd'])
    engine_string = engine_string.replace('host', options['host'])
    engine_string = engine_string.replace('dbname', options['testdb'])
    engine = create_engine(engine_string)

    myinputdb = pymysql.connect(host=options['host'], user=options['user'], passwd=options['passwd'], db=options['db'], port=int(options['port']))
    myoutputdb = pymysql.connect(host=options['host'], user=options['user'], passwd=options['passwd'], db=options['testdb'], port=int(options['port']))
    mycursor = myoutputdb.cursor()

    writespat='insert into OnlineSPaT(intersection_id, camera_id, timestamp, hexphase, cycle, city, state, type) values (%s, %s, %s, %s, %s, %s, %s, %s);'
    writequery='insert into RealDisplayInfo (frame_id, track_id, center_x, center_y, class, timestamp, intersection_id, camera_id, SPAT, Cycle, speed_x, speed_y, unique_ID, skip_begin, skip_end, skip_mid, skip_angle, raw_x, raw_y) values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);'
    writetrackquery='insert into RealTrackProperties (intersection_id, camera_id, timestamp, track_id, phase, cluster, isAnomalous, redJump, nearmiss, signature, unique_ID, small, outside, hispeed, oddshape) values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);'
    updatetrackquery='update RealTrackProperties set phase=%s, cluster=%s, isAnomalous=%s, redJump=%s, small=%s, outside=%s, hispeed=%s, oddshape=%s where timestamp=%s and track_id=%s and unique_ID=%s;'
    query = 'select * from testATSPM where timestamp between %(name)s and %(name2)s and intersection_id=%(name3)s;'
    spatquery = 'select timestamp,hexphase,cycle,type from OnlineSPaT where timestamp between %(name2)s and %(name3)s and intersection_id=%(name)s;'

    cindex = 4
    phases = []
    signal_violate = []

    minstep = datetime.timedelta(minutes=1)
    maxphases = 80
    start = time.time()
    stime = pd.to_datetime(options['start'])
    ftime = stime + minstep
    etime = pd.to_datetime(options['end'])
    intersection_num = int(options['cameraid'])
    #query = 'select time,frame_id,track_id,center_x,center_y,class,timestamp,intersection_id,nearmiss,signature from OnlineTrackInfo where ID <= %(name)s and intersection_id=6 and timestamp > %(name2)s;'
    querylimit = 'select time,frame_id,track_id,center_x,center_y,class,timestamp,intersection_id,nearmiss,signature,unique_ID from OnlineTrackInfo where timestamp between %(name2)s and %(name3)s and intersection_id=%(name)s;'
    videodf = pd.read_sql(querylimit, myinputdb, params={'name' : str(int(options['videocameraid'])), 'name2' : stime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
    tracks = videodf['track_id'].unique()
    ts = videodf.loc[0][6] #change later to now = datetime.now()f
    #ets = videodf.tail(1)[6]
    currentTime = pd.to_datetime(ts) #, format = '%m/%d/%y %H:%M:%S.%f')
    startTime = pd.to_datetime(stime) - datetime.timedelta(minutes=15)
    sdf = pd.read_sql(query, myoutputdb, params={'name': startTime.strftime("%Y-%m-%d %H:%M:%S"), 'name2' : ftime.strftime("%Y-%m-%d %H:%M:%S"), 'name3' : options['cityintid']})
    spatrows = convertToSPaT(sdf, startTime)
    mycursor.executemany(writespat, spatrows)
    myoutputdb.commit()
    videodf = pd.read_sql(querylimit, myinputdb, params={'name' : str(int(options['videocameraid'])), 'name2' : stime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
    spatdf = pd.DataFrame(spatrows, columns =['intersection_id', 'camera_id', 'timestamp', 'hexphase', 'cycle', 'city', 'state','type'])
    spatdf['timestamp'] = pd.to_datetime(spatdf['timestamp'])
    spatdf['timestamp'] = spatdf['timestamp'].dt.round('100ms')
    tracks = videodf['track_id'].unique()
    recs = 0
    updaterows = []
    insertrows = []
    outputdf = pd.DataFrame()
    inserttrackrows = []
    phase2vehiclets = []
    ci = int(options['clusterinterval'])
    clusterinterval = datetime.timedelta(minutes=ci)
    wait = datetime.timedelta(minutes=ci)
    t2r = 0
    t2g = 0
    fftime = 0
    stopbar = phase2stopbar
    offset = datetime.timedelta(seconds=int(options['offset']))
    threshtime = datetime.timedelta(seconds=30)
    geoClusters = readCentroidData(myoutputdb)
    clusternow = ftime + datetime.timedelta(minutes=int(options['collecttime']))
    while (ftime < etime):
        if (ftime > clusternow):
            geoClusters = readCentroidData(myoutputdb)
            clusternow = ftime + datetime.timedelta(minutes=int(options['clusterinterval'])-int(options['collecttime']))
        results = []
        myimpl = time.time()
        with mp.Pool(10) as pool:
            for tid in tracks:
                results.append(pool.apply_async(processTrack, args=(tid, videodf, spatdf, geoClusters, finalCentroidName, finalCentroids, tracksH)))
                #updaterows, outputdf, inserttrackrows, nt, ignore = processTrack(tid, videodf, spatdf, geoClusters, finalCentroidName, finalCentroids, tracksH)
                #if (ignore == 0):
                    #tracksH[tid] = nt
            prevstime = stime
            prevftime = ftime
            prevspatdf = spatdf
            stime = ftime + step
            ftime = stime + minstep
            mysqldfread = time.time()
            videodf = pd.read_sql(querylimit, myinputdb, params={'name' : str(int(options['videocameraid'])), 'name2' : stime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
            mysqldfreaddone = time.time()
            mysqlsdfread = time.time()
            sdf = pd.read_sql(query, myoutputdb, params={'name': stime.strftime("%Y-%m-%d %H:%M:%S"), 'name2' : ftime.strftime("%Y-%m-%d %H:%M:%S"), 'name3' : options['cityintid']})
            mysqlsdfreaddone = time.time()
            spatrows = convertToSPaT(sdf, stime)
            mytimespatrows =  time.time()
            if (spatrows == []):
                spatdf = prevspatdf
            else:
                mycursor.executemany(writespat, spatrows)
                mytimespatwrite =  time.time()
                startTime = ftime - datetime.timedelta(minutes=10)
                #spatdf = pd.DataFrame(spatrows, columns =['intersection_id', 'camera_id', 'timestamp', 'hexphase', 'cycle', 'city', 'state', 'type'])
                #spatdf['timestamp'] = pd.to_datetime(spatdf['timestamp'])
            mytimespatread = time.time()
            spatdf = pd.read_sql(spatquery, myoutputdb, params={'name' : options['cityintid'], 'name2' : startTime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
            mytimespatreaddone = time.time()
            spatdf['timestamp'] = spatdf['timestamp'].dt.round('100ms')
            mytimespatround = time.time()
            #print ("Finished processing 1 minute worth of records", (mysqldfreaddone-mysqldfread), 'df read')
            #print ("Finished processing 1 minute worth of records", (mysqlsdfreaddone-mysqlsdfread), 'sdf read', (mytimespatrows-mysqlsdfreaddone), 'convertToSPaT')
            #print ("Finished processing 1 minute worth of records", (mytimespatwrite-mytimespatrows), 'db write spat', (mytimespatreaddone-mytimespatread), 'spatdf read', (mytimespatround-mytimespatreaddone), 'spat round deci')
            myimpldone = time.time()
            pool.close()
            pool.join()
        mypropcopy = time.time()
        for i in range(0, len(tracks)):
            gotresults = results[i].get()
            ignore = gotresults[4]
            if (ignore == 0):
                updaterows.extend(gotresults[0])
                outputdf = outputdf.append(gotresults[1])
                inserttrackrows.extend(gotresults[2])
                tracksH[tracks[i]] = gotresults[3]
        #print (outputdf)
        outputdf = outputdf.astype({"unique_ID":'int64'})
        mypropcopydone = time.time()
        mysqlws = time.time()
        outputdf.sort_values(by='frame_id')
        outputdf.to_sql('RealDisplayInfo', con=engine, if_exists='append', index = False)
        outputdf = outputdf.iloc[0:0]
        mycursor.executemany(updatetrackquery, updaterows)
        #print ('Here are the insert rows')
        print (inserttrackrows)
        mycursor.executemany(writetrackquery, inserttrackrows)
        mysqlwb = time.time()
        myoutputdb.commit()
        mysqlwe = time.time()
        #print ("Finished processing 1 minute worth of records", mysqlwb-mysqlws, 'track write', mysqlwe-mysqlwb, 'displayinfo write')
        #print ("Finished processing 1 minute worth of records",mysqldfread-myimpl, 'impl', len(tracks), (mypropcopy-myimpldone), 'pool', (mypropcopydone-mypropcopy), 'copy res')
        updaterows.clear()
        insertrows.clear()
        inserttrackrows.clear()
        spatrows.clear()
        wait = wait + step + minstep
        if (wait > clusterinterval):
            sigtime, phasechangestatus = getPhase2Signal(prevstime, prevftime, prevspatdf)
            if (phasechangestatus):
                tenseconds = datetime.timedelta(seconds=10)
                lowthresh = np.datetime64(sigtime - threshtime)
                highthresh = np.datetime64(sigtime + tenseconds)
                threshtime = tenseconds
                for tid in tracks:
                    nt = tracksH[tid]
                    if (nt != None):
                        ts = staticTrackCrossesStopBar(nt, stopbar, sigtime)
                        if (ts != 0 and ts > lowthresh and ts < highthresh):
                            phase2vehiclets.append(ts)
                if (phase2vehiclets):
                    phase2vehiclets.sort()
                    lowestts = pd.to_datetime(phase2vehiclets[0])
                    offset = offset + sigtime-(lowestts-datetime.timedelta(seconds=2))
                    print ('offset is', offset)
                    wait = datetime.timedelta(seconds=0)
                phase2vehiclets.clear()
        videodf.timestamp = videodf.timestamp + offset
        tracks = videodf['track_id'].unique()
        end = time.time()
        print (end-myimpl)
#
    mycursor.close()
    myoutputdb.close()

def getPhase2Signal(stime, ftime, spatdf):
    status = False
    ts = 0
    vehspat = spatdf[spatdf['type'] == 1]
    df = vehspat[vehspat['timestamp'] >= stime]
    mask = df.hexphase.str[0] == '4'
    pos = np.flatnonzero(mask)
    if (pos.size != 0):
        ts = df.iloc[pos[0]]['timestamp'] 
        prevphase = df.iloc[pos[0]-1]['hexphase']
        bin_string = bin(int('1'+prevphase, 16))[3:]
        if (bin_string[1] == '0' and bin_string[4] == '0'):
            status = True
    return ts, status

def getfirstvehphase2(t2g, t2r, phase2vehiclets):
    llen = len(phase2vehiclets)
    firstveh = 0
    for i in range(0, llen-1):
        ts1 = phase2vehiclets[i]
        ts2 = phase2vehiclets[i+1]
        diff = ts2 - ts1
        difflight = t2g - t2r
        if (diff >= difflight):
            firstveh = ts2
            break
    return firstveh

def timeTrackCrossesStopBar(nt, stopbar): #the car must have arrived before the signal turned green which is checked by speed and a comparison with stime
    ts = 0
    xx = nt.x
    yy = nt.arr
    tss = nt.ts
    fss = nt.f
    xlen = xx.size
    for i in range(0,xlen):
        x = xx[i]
        y = yy[i]
        x1 = stopbar[0][0]
        y1 = stopbar[0][1]
        x2 = stopbar[1][0]
        y2 = stopbar[1][1]
        #must be arranged in increasing order of x
        d = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
        if (i > 0):
            sameSign = (d * d1) > 0
            if (sameSign == False):
                ts = tss[i] + fss[i]%10*np.timedelta64(100,'ms')
                break
        d1 = d
    return ts

def staticTrackCrossesStopBar(nt, stopbar, stime): #the car must have arrived before the signal turned green which is checked by speed and a comparison with stime
    ts = 0
    if ((nt.phase == 2 or nt.phase == 5) and pd.to_datetime(nt.ts[0]) < stime):
        if (nt.x.size < 10):
            return ts
        sp = np.asarray([[x, y] for x, y in zip(nt.sx, nt.sy)])
        sp2 = np.sqrt(np.einsum('ij,ij->i', sp, sp))
        spzero = np.argwhere(sp2==0)
        if (spzero.size <= 1):
            return ts
        xx = nt.x
        yy = nt.arr
        tss = nt.ts
        xlen = xx.size
        for i in range(0,xlen):
            x = xx[i]
            y = yy[i]
            x1 = stopbar[0][0]
            y1 = stopbar[0][1]
            x2 = stopbar[1][0]
            y2 = stopbar[1][1]
            #must be arranged in increasing order of x
            d = (x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)
            if (i > 0):
                sameSign = (d * d1) > 0
                if (sameSign == False):
                    ts = tss[i]
                    break
            d1 = d
    return ts

def allocateClustersToGrid(geoClusters):
    clen = len(finalCentroids)
    for j in range(0, clen):
        cent = finalCentroids[j]
        xlen = cent.x.size
        i = 0
        while (i < xlen-1):
            v, h = locatePointOnGrid(cent.x[i], cent.arr[i])
            geoClusters[v][h].append(finalCentroidName[j])
            #advance i
            #i = advanceForPointsOnSameCell(i, cent.x, cent.arr, h, v)
            i = i+1
    lenhlines = len(hlines)
    lenvlines = len(vlines)

    for i in range(0, lenvlines):
        for j in range (0, lenhlines):
            glen = len(geoClusters[i][j])
            if glen > 1:
                indexlist = []
                for d in range (0, glen-1):
                    if (d in indexlist):
                        continue
                    cname0 = geoClusters[i][j][d]
                    for k in range (d+1, glen):
                        if (cname0 ==  geoClusters[i][j][k]):
                            #collapse k
                            indexlist.append(k)
                indexlist.sort(reverse=True)
                for d in indexlist:
                    del geoClusters[i][j][d]

def isInDifferentCell(x, y, v, h):
    status = True
    if (v > 0 and h > 0):
        vv1 = [lowerpvlines[v-1][0], lowerpvlines[v-1][1]]
        vv2 = [upperpvlines[v-1][0], upperpvlines[v-1][1]]
        v1 = [lowerpvlines[v][0], lowerpvlines[v][1]]
        v2 = [upperpvlines[v][1], upperpvlines[v][1]]
        p = [x, y]
        dvv = (x-vv1[0])*(vv2[1] - vv1[1]) - (y-vv1[1])*(vv2[0] - vv1[0])
        dv = (x - v1[0])*(v2[1] - v1[1]) - (y - v1[1]) * (v2[0] - v1[0])
        #print ('isInDifferentCell', x, y, dvv, dv)
        if (dv * dvv < 0):
            hh1 = [leftphlines[h-1][0], leftphlines[h-1][1]]
            hh2 = [rightphlines[h-1][0], rightphlines[h-1][1]]
            h1 = [leftphlines[h][0], leftphlines[h][1]]
            h2 = [rightphlines[h][0], rightphlines[h][1]]
            dhh = (x-hh1[0])*(hh2[1] - hh1[1]) - (y-hh1[1])*(hh2[0] - hh1[0])
            dh = (x-h1[0])*(h2[1] - h1[1]) - (y-h1[1])*(h2[0] - h1[0])
            #print ('isInDifferentCell', x, y, dhh, dh)
            if (dhh * dh < 0):
                status = False
    #print ('isInDifferentCell', status)
    return status

def findMatchingCentroid(track, geoclusters, centroid_list, centroids):
    if (len(centroids) == 0):
        return track
    x = track.x
    y = track.arr
    xlen = x.size
    cset = set(centroid_list)
    clist = []
    i = 0
    prevx = 0
    prevy = 0
    prevv = -1
    prevh = -1
#    while (i < xlen-1):
#        if ((prevx != int(x[i]) or prevy != int(y[i])) and np.sqrt(np.square(x[i]-prevx) + np.square(y[i]-prevy)) > 20):
#            #if (prevv == -1 or prevh == -1 or isInDifferentCell(x[i], y[i], prevv, prevh)):
#            v, h = locatePointOnGrid(x[i], y[i])
#            #print (track.tid, x[i], y[i], v, h, len(vlines)/2, len(hlines)/2, vlines[v], hlines[h], prevx, prevy, np.sqrt(np.square(x[i]-prevx) + np.square(y[i]-prevy)))
#            prevv = v
#            prevh = h
#            clist.extend(geoclusters[v][h])
#            prevx = int(x[i])
#            prevy = int(y[i])
#        i = i+1
#    cset = set(clist)
    found = 0
    mindist = 500
    mincent = ''
    min_centroid_track = centroids[0]
    distace_thresh = 30
    first_thresh = 19
    last_threah = 19
    if (isCurved(track)):
        track.isthrough = False
        distace_thresh = 39
        first_thresh = 20
        last_threah = 22
    if (track.classtype == 'pedestrian'):
        distace_thresh = 50
        first_thresh = 30
        last_thresh = 30
    minimax = 0
    for cent in cset:
        isPed = 0
        idx = centroid_list.index(cent)
        centroid_track = centroids[idx]
        #if (track.gbit == centroid_track.gbit and track.gbit2 == centroid_track.gbit2):
        if (centroid_track.classtype == track.classtype): # and centroid_track.phase == track.phase):
            #isRelated, dist = closeToMembers(track, centroid_track, isPed)
            tdist, maxdiff, dist, lenratio, firstdiff, lastdiff = distance(track, centroid_track, 200)
            tdist1, maxdiff1, dist, lenratio, firstdiff1, lastdiff1 = distance(centroid_track, track, 200)
            (d, f, l) = ((tdist+tdist1)/2, min(firstdiff, firstdiff1), min(lastdiff, lastdiff1))
            print ("Centroid determination:", track.tid, cent, d, f, l, min(maxdiff,maxdiff1))
            if (d < distace_thresh and f < first_thresh and l < last_threah):
                found = 1
                if (d < mindist):
                    mindist = d
                    mincent = cent
                    idx = centroid_list.index(cent)
                    centroid_track = centroids[idx]
                    min_centroid_track = centroid_track
                    minimax = min(maxdiff,maxdiff1)
    if (found == 0):
        #print ('match_result', track.tid, 'anomalous')
        track.isAnomalous = 1
        track.clus ="shape"
        track.oddshape = 1
    else:
        #print ('match_result', track.tid, mincent[6:], track.phasetype, mincent)
        track.clus = mincent
        track.phase = min_centroid_track.phase
        if (track.phase > 8):
            track.phase = 2*(track.phase-8)

    text = track.clus
    if (minimax > 70):
        text = "convert to shape"
    print ("Final Centroid determination:", track.tid, track.clus, track.isthrough, text)
    return min_centroid_track
if __name__=="__main__":
        main()
