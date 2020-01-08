import numpy as np
import pandas as pd
import math
import datetime
import time
import pymysql
import sys
import multiprocessing as mp
from ast import literal_eval
from hash import HashTable
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
from shapely.geometry import Point, Polygon

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

directions = {}
directions['NS'] = literal_eval(options['NS'])
directions['SN'] = literal_eval(options['SN'])
directions['EW'] = literal_eval(options['EW'])
directions['WE'] = literal_eval(options['WE'])
dir_2_x_2 = directions['NS'][0][0]  #point North in NS
dir_2_y_2 = directions['NS'][0][1]  #point North in NS
dir_6_x_2 = directions['NS'][1][0]  #point Sount in NS
dir_6_y_2 = directions['NS'][1][1]  #point Sount in NS
dir_2_x_6 = directions['SN'][0][0]  #point North in SN
dir_2_y_6 = directions['SN'][0][1]  #point North in SN
dir_6_x_6 = directions['SN'][1][0]  #point Sount in SN
dir_6_y_6 = directions['SN'][1][1]  #point Sount in SN
dir_4_x_8 = directions['EW'][0][0]  #point East in EW
dir_4_y_8 = directions['EW'][0][1]  #point East in EW
dir_8_x_8 = directions['EW'][1][0]  #point West in EW
dir_8_y_8 = directions['EW'][1][1]  #point West in EW
dir_8_x_4 = directions['WE'][0][0]  #point East in WE
dir_8_y_4 = directions['WE'][0][1]  #point East in WE
dir_4_x_4 = directions['WE'][1][0]  #point East in WE
dir_4_y_4 = directions['WE'][1][1]  #point East in WE

polylist = [literal_eval(options['polygon'])]
flat_list = [item for sublist in polylist for item in sublist]
polygon = Polygon(flat_list)
pedpolygon = polygon

finalCentroids = []
finalCentroidName = []
objects=options['objects'].split(',')

#@dataclass
class Track:
    def __init__(self, tid, x, y, rrange):
        self.tid = tid        #track id
        self.rrange = rrange
        self.x = x
        self.arr = y
        self.t = np.empty(shape=rrange)
        self.ts = np.empty(shape=rrange)
        self.f = np.empty(shape=rrange)
        self.cid = -1
        self.phasetype = '' #major, minor, leftmajor, leftminor
        self.classtype = '' #other types are bus, truck, ped
        self.cycles = []      # cycle numbers for which the vehicle is present in the intersection
        self.slope = 0
        self.isAnomalous = 0
        self.isrj = 0
        self.start = 0
        self.end = 0
        self.scale = 1
        self.len = 0
        self.iid = 0
        self.clus = 0
        self.iscent = 0
        self.isthrough = True

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
        self.pedstatus = ped   #1 if ped call recorded, 0 otherwise
        self.skipDetector = 0
        self.finish = 0
        self.sIndex = 0
        self.eIndex = 0
        self.cumaor = 0
        self.cumaog = 0

centdf = pd.read_csv(options['centroid'])
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
np.seterr(all='raise')
tracksH = HashTable(50000)

def swap (x1, x2):
    tmp = x1
    x1 = x2
    x2 = tmp
    return x1, x2

#vlines = [[(0, 0), (0, 960)], [(332.973, 646.486), (315.676, 315.676)], [(436.757, 261.622), (490.811, 741.622)], [(518.919, 263.784), (577.297, 745.946)], [(605.405, 261.622), (640.743, 743.784)], [(737.297, 274.595), (750.27, 683.243)], [(880.339, 459), (882.162, 611.892)], [(1280, 0), (1280, 960)]]
#hlines = [[(322.162, 0), (877.838, 0)], [(322.162, 352.432), (877.838, 352.432)], [(322.162, 406.486), (910.27, 406.486)], [(332.541, 454.054), (966.486, 454.054)], [(322.162, 516.757), (877.838, 508.108)], [(321.162, 562.162), (877.838, 560)], [(322.162, 607.568), (877.838, 611.892)], [(402.162, 735.135), (715.676, 735.135)], [(332.160, 960), (332.160, 960)]]
#https://math.stackexchange.com/questions/274712/calculate-on-which-side-of-a-straight-line-is-a-given-point-located
vlines=literal_eval(options['verticalgrid'])
hlines=literal_eval(options['horizontalgrid'])
geoClusters = [[[] for j in range(len(hlines))] for i in range(len(vlines))]
def locatePointOnGrid(x, y):
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

    d1 = 1
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

def distance(nt, nt1, boxDim):
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

    if (isCurved(big)):
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
                dist = np.absolute((y2-y1) * x0 - (x2 - x1) * y0 + x2*y1 - y2*x1)/(np.sqrt(np.add(np.square(y2-y1), np.square(x2-x1))))
                status = dist < 32
                print ('isOnSameLane', mindiststart, mindistlast, dist)
            else:
                if (idx < big.x.size-2):
                    x2 = big.x[idx+1]
                    y2 = big.arr[idx+1]
                    a = getCosineOfAngle(x1-x2, y1-y2, x0-xn, y0-yn)
                    if (a > 0.85):
                        dist = np.absolute((y2-y1) * x0 - (x2 - x1) * y0 + x2*y1 - y2*x1)/(np.sqrt(np.add(np.square(y2-y1), np.square(x2-x1))))
                        status = dist < 32
                        print ('isOnSameLane', mindiststart, mindistlast, dist)

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
    c1 = isCurved(nt1)
    c2 = isCurved(nt2)
    if (c1 == False and c2 == False):
        threshold = 0.98
    elif (c1 == True and c2 == True):
        biggerTrack, smallertrack, status = isLengthDifferent(nt1, nt2)
        lb = biggerTrack.len
        ls = smallertrack.len
        if (isCurved(biggerTrack) and lb > 250 and ls > 0):
            threshold = 0.70
            tdistthresh = 60 * lb/ls
        else:
            threshold = 0.75 #getSweep(biggerTrack)
            tdistthresh = 60
    else:
        biggerTrack, smallertrack, status = isLengthDifferent(nt1, nt2)
        lb = biggerTrack.len
        ls = smallertrack.len
        if (isCurved(biggerTrack) and lb > 250 and ls > 0):
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
            tdist, maxdiff, dtwdist, ovlp = distance(centi, centj, 1000)
            tdist = tdist * scale
            if (tdist > 75 or dtwdist > 12000):
                status = False
            print(centi.tid, centj.tid, tdist, maxdiff, dtwdist, ovlp, direction, status)
    else:
        threshold, tdistthresh = getTrackMatchingThresholdForMembers(centi, centj)
        if (direction < threshold):
            status = False
        else:
            isSameLane = isOnSameLane(centi, centj)
            if (isSameLane == False):
                status = False
            else:
                tdist, maxdiff, dtwdist, ovlp = distance(centi, centj, 1000)
                tdist = tdist * scale
                #if (tdist > 15 or maxdiff > 200 or (dtwdist > 10000 and maxdiff < 20)):  # second half is an indication that the path diverges out
                    #status = False
                if (tdist > tdistthresh):
                    if (tdist < 50 and isCollinear(centi, centj)): 
                        nt1 = centi
                        nt2 = centj
                        n1 = nt1.x.size
                        n2 = nt2.x.size
                        if (tdist < 25):
                            status = True
                        elif (nt1.x[n1-1] < nt2.x[0] or nt2.x[n2-1] < nt1.x[0] and strictly_increasing(nt1.x)):
                            status = True
                        elif (nt1.x[n1-1] > nt2.x[0] or nt2.x[n2-1] > nt1.x[0] and strictly_decreasing(nt1.x)):
                            status = True
                        else:
                            status = False
                    #elif (tdist < 20 and maxdiff < 25):
                        #status = True
                    elif (tdist < 45 and isCurved(centj)):
                        status = True
                    else:
                        status = False
                print(centi.tid, centj.tid, tdist, maxdiff, dtwdist, ovlp, direction, status)
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
    [Ax, Ay] = [x[x.size-1]-x[0], y[y.size-1]-y[0]]
    #Check for 2-> or 6->2 alignment
    a = getCosineOfAngle(Ax, Ay, dir_6_x_6 - dir_2_x_6, dir_6_y_6 - dir_2_y_6)
    throughcos = float(options['through'])
    curvedcos = float(options['curved'])
    if (a > throughcos):
        bit = 2
    else:
        a = getCosineOfAngle(Ax, Ay, dir_2_x_2 - dir_6_x_2, dir_2_y_2 - dir_6_y_2)
        if (a > throughcos):
            bit = 6
        else:
            a = getCosineOfAngle(Ax, Ay, dir_8_x_8 - dir_4_x_8, dir_8_y_8 - dir_4_y_8)
            if (a > throughcos):
                bit = 4
            else:
                a = getCosineOfAngle(Ax, Ay, dir_4_x_4 - dir_8_x_4, dir_4_y_4 - dir_8_y_4)
                if (a > throughcos):
                    bit = 8
                else: # check curved paths
                    a = getCosineOfAngle(Ax, Ay, dir_4_x_8-dir_2_x_6, dir_4_y_8-dir_2_y_6)
                    if (a > curvedcos):
                        if (np.absolute(y[0]-dir_2_y_6) < 50): #np.absolute needed to differentiate between left turn and the far right turn
                            bit = 5
                            bit2 = 2
                        else:
                            bit = 0
                            bit2 = 8
                    else:
                        a = getCosineOfAngle(Ax, Ay, dir_8_x_8-dir_6_x_2, dir_8_y_8-dir_6_y_2)
                        if (a > curvedcos):
                            if (np.absolute(y[0]-dir_6_y_6) < 50):
                                bit = 1
                                bit2 = 6
                            else:
                                bit = 0
                                bit2 = 4
                        else:
                            a = getCosineOfAngle(Ax, Ay, dir_6_x_2-dir_4_x_8, dir_6_y_2-dir_4_y_8)
                            if (a > curvedcos):
                                if (np.absolute(x[0]-dir_4_x_8) < 50):
                                    bit = 7
                                    bit2 = 4
                                else:
                                    bit = 0
                                    bit2 = 2
                            else:
                                a = getCosineOfAngle(Ax, Ay, dir_2_x_6-dir_8_x_8, dir_2_y_6-dir_8_y_8)
                                if (a > curvedcos):
                                    if (np.absolute(x[0]-dir_8_x_8) < 50):
                                        bit = 3
                                        bit2 = 8
                                    else:
                                        bit = 0
                                        bit2 = 6
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
        splpixels = convertToPixelsPerSec(speedlimitwithbuffer)
        highsp = np.where(sp2>np.sqrt(splpixels))
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
                print ('Signs of track jump', nt.tid, a)
            else:
                maxspeed = convertToPixelsPerSec(speedlimit*2)
                if (np.amax(sp2) > maxspeed):
                    check = True
                    print ('Signs of track jump', nt.tid, a)
        index = np.argmax(sp2 > 50)
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
previousTime = -1
phaseInfo = []
numPhases = 16
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
                roundedtime = reportTime.round('1s')
                phasestr = '0' * paddedzeroes + hexstr
                spatrow = (options['cityintid'], roundedtime.strftime("%Y-%m-%d %H:%M:%S"), phasestr, str(cycle_no), 'Gainesville','FL')
                spatrows.append(spatrow)
            reportTime = reportTime + step
    for row in eventList:
        if (row[eventCIndex] == 1):     #Phase begin green
            if (isBaseSet == 0):
                baseTime = currentTime
                isBaseSet = 1
            if (i == 0):       #The first phase turns green in the CSV
                i = i+1
                cycle_mark = row[eventPIndex]
                reportTime = currentTime.round('1s') + step
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
    if (reportTime != currentTime):
        reportTime = processEventList(eventList, previousTime, reportTime, spatrows)
        eventList.clear()
        previousTime = currentTime
    return spatrows

def processTrack(tid, df, spatdf, geoclusters, centroid_list, centroids, tH):#, mycursor, myoutputdb):
    updaterows = []
    insertrows = []
    updateindex = -1
    oldclus = 'anom'
    ndf = df.loc[df['track_id'] == tid]
    t = ndf.iloc[:,0].values
    f = ndf.iloc[:,1].values
    x = ndf.iloc[:,3].values
    y = ndf.iloc[:,4].values
    c = ndf.iloc[:,5].values
    ts = ndf.iloc[:,6].values
    iid = ndf.iloc[:,7].values
    nm = ndf.iloc[:,8].values
    sig = ndf.iloc[:,9].values
    y = np.ma.compressed(np.ma.masked_where(x==0, y))
    t = np.ma.compressed(np.ma.masked_where(x==0, t))
    f = np.ma.compressed(np.ma.masked_where(x==0, f))
    ts = np.ma.compressed(np.ma.masked_where(x==0, ts))
    nm = np.ma.compressed(np.ma.masked_where(x==0, nm))
    sig = np.ma.compressed(np.ma.masked_where(x==0, sig))
    x = np.ma.compressed(np.ma.masked_where(x==0, x))
    if (x.size == 0):
        nt = Track(tid, x, y, 1260)
        nt.classtype = c[0]
        nt.iid = iid[0]
        nt.t = t
        nt.ts = ts
        nt.f = f
        return updaterows, insertrows, nt
    phaseindex = ts.size - 2
    if (phaseindex < 0):
        phaseindex = ts.size - 1
    exactmatch=spatdf[spatdf['timestamp']==ts[phaseindex]]
    if (exactmatch.empty == False):
        rs = exactmatch['hexphase'].values[0]
        cycle = exactmatch['cycle'].values[0]
    else:
        temp = ts[phaseindex]
        while (exactmatch.empty):
            exactmatch=spatdf[spatdf['timestamp']==(temp - np.timedelta64(1,'s'))]
            temp = temp - np.timedelta64(1,'s')
        if (exactmatch.empty == False):
            rs = exactmatch['hexphase'].values[0]
            cycle = exactmatch['cycle'].values[0]
    nt = tH[tid]
    if (nt != None and nt.iscent != 1):
        updateindex = nt.x.size
        x = np.concatenate([nt.x, x])
        y = np.concatenate([nt.arr, y])
        t = np.concatenate([nt.t, t])
        ts = np.concatenate([nt.ts, ts])
        f = np.concatenate([nt.f, f])
        oldclus = nt.clus
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

    rstart = 0
    rend = 0
    rrange = 1260
    nt = Track(tid, x, y, rrange)
    nt.classtype = c[0]
    nt.phasetype = rs
    nt.len = getPathLength(tid, x, y)
    nt.iid = iid[0]
    nt.t = t
    nt.ts = ts
    nt.f = f
    cent = findMatchingCentroid(nt, geoclusters, centroid_list, centroids)
    sx, sy, spindex, check = setSpeed(x, y, t, nt, cent)
    if (check):
        nt.isAnomalous = 1
    isrj = 0
    speedphase = rs

    if (spindex > -1):
        exactmatch=spatdf[spatdf['timestamp']==ts[spindex]]
        if (exactmatch.empty == False):
            speedphase = exactmatch['hexphase'].values[0]
        else:
            temp = ts[spindex]
            while (exactmatch.empty):
                exactmatch=spatdf[spatdf['timestamp']==(temp - np.timedelta64(1,'s'))]
                temp = temp - np.timedelta64(1,'s')
            if (exactmatch.empty == False):
                speedphase = exactmatch['hexphase'].values[0]
    si = 0
    se = np.size(x)
    r = np.size(x)
    phasetype=rs
    pathLength = getPathLength(tid, x, y)
    cindex = objects.index(c[0])
    pedindex = objects.index('pedestrian')
    if (cindex != pedindex):
        rbit, rbit2 = getGreenBit(x, y)
        if (rbit >= 1):
            bin_string = bin(int('1'+speedphase, 16))[3:]
            if (bin_string[rbit+16-1] == '1'):
                if (rbit2 > 0):
                    if (bin_string[rbit2+16-1] == '1'):
                        isrj = 1
                else:
                    isrj = 1
    nt.isrj = isrj
    si = 0
    writequery='insert into OnlineDisplayInfo2(frame_id, track_id, center_x, center_y, class, timestamp, intersection_id, SPAT, Cycle, speed_x, speed_y, cluster, isAnomalous, redJump, nearmiss, signature) values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);'
    updatequery='update OnlineDisplayInfo2 set cluster=%s, isAnomalous=%s, redJump=%s where frame_id=%s and track_id=%s;'
    if (updateindex > -1):
        si = updateindex
        if (oldclus != nt.clus):
            for index in range(0, si):
                pdts = pd.to_datetime(ts[index])
                updaterow = (nt.clus, str(nt.isAnomalous), str(nt.isrj), str(f[index]), str(tid))
                updaterows.append(updaterow)
                #mycursor.execute(updatequery, updaterow)
                #myoutputdb.commit()
    for index in range(si, ts.size-1):
        pdts = pd.to_datetime(ts[index])
        tstr = pdts.strftime("%Y-%m-%d %H:%M:%S") + '.' + str(f[index]%10)
        insertrow = (str(f[index]), str(tid), str(x[index]), str(y[index]), str(c[index-si]), tstr, str(iid[index-si]), rs, str(cycle), str(sx[index-si]), str(sy[index-si]), nt.clus, str(nt.isAnomalous), str(nt.isrj), str(nm[index-si]), str(sig[index-si]))
        insertrows.append(insertrow)
        #mysqlwb = time.time()
        #mycursor.execute(writequery, insertrow)
        #myoutputdb.commit()
    return updaterows, insertrows, nt

def main():
    cindex = 4
    phases = []
    centroid_tracks = centdf['track_id'].unique()
    signal_violate = []
    for tid in centroid_tracks:
        ndf = centdf.loc[centdf['track_id'] == tid]
        x = ndf.iloc[:,1].values
        y = ndf.iloc[:,2].values
        c = ndf.iloc[:,3].values
        iid = ndf.iloc[:,4].values
        ps = ndf.iloc[:,5].values #phases
        name = ndf.iloc[:,6].values
        r = np.size(x)
        cindex = objects.index(c[0])
        r = np.size(x)
        nt = Track(tid, x, y, 0)
        nt.start = 0
        nt.end = x.size
        nt.classtype = c[0]
        nt.phasetype = ps[r-2]
        nt.len = getPathLength(tid, x, y)
        nt.iid = iid[0]
        nt.iscent = 1
        tracksH[tid] = nt
        finalCentroids.append(nt)
        finalCentroidName.append(name[0])
        if (isCurved(nt)):
            nt.isthrough = False
    allocateClustersToGrid(geoClusters)

    myinputdb = pymysql.connect(host=options['host'], user=options['user'], passwd=options['passwd'], db=options['db'], port=int(options['port']))
    myoutputdb = pymysql.connect(host=options['host'], user=options['user'], passwd=options['passwd'], db=options['testdb'], port=int(options['port']))
    mycursor = myoutputdb.cursor()
    writespat='insert into OnlineSPaT(intersection_id, timestamp, hexphase, cycle, city, state) values (%s, %s, %s, %s, %s, %s);'
    writequery='insert into OnlineDisplayInfo2 (frame_id, track_id, center_x, center_y, class, timestamp, intersection_id, SPAT, Cycle, speed_x, speed_y, cluster, isAnomalous, redJump, nearmiss, signature) values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);'
    updatequery='update OnlineDisplayInfo2 set cluster=%s, isAnomalous=%s, redJump=%s where frame_id=%s and track_id=%s;'
    query = 'select * from testATSPM where timestamp between %(name)s and %(name2)s and intersection_id=%(name3)s;'
    spatquery = 'select timestamp,hexphase,cycle from OnlineSPaT where timestamp between %(name2)s and %(name3)s and intersection_id=%(name)s;'

    minstep = datetime.timedelta(minutes=1)
    maxphases = 80
    start = time.time()
    stime = pd.to_datetime(options['start'])
    ftime = stime + minstep
    etime = pd.to_datetime(options['end'])
    intersection_num = int(options['id'])
    #query = 'select time,frame_id,track_id,center_x,center_y,class,timestamp,intersection_id,nearmiss,signature from OnlineTrackInfo where ID <= %(name)s and intersection_id=6 and timestamp > %(name2)s;'
    querylimit = 'select time,frame_id,track_id,center_x,center_y,class,timestamp,intersection_id,nearmiss,signature from OnlineTrackInfo where timestamp between %(name2)s and %(name3)s and intersection_id=%(name)s;'
    df = pd.read_sql(querylimit, myinputdb, params={'name' : str(intersection_num), 'name2' : stime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
    tracks = df['track_id'].unique()
    ts = df.loc[0][6] #change later to now = datetime.now()f
    #ets = df.tail(1)[6]
    currentTime = pd.to_datetime(ts) #, format = '%m/%d/%y %H:%M:%S.%f')
    startTime = pd.to_datetime(stime) - datetime.timedelta(minutes=15)
    sdf = pd.read_sql(query, myoutputdb, params={'name': startTime.strftime("%Y-%m-%d %H:%M:%S"), 'name2' : ftime.strftime("%Y-%m-%d %H:%M:%S"), 'name3' : options['cityintid']})
    spatrows = convertToSPaT(sdf, startTime)
    mycursor.executemany(writespat, spatrows)
    myoutputdb.commit()
    df = pd.read_sql(querylimit, myinputdb, params={'name' : str(intersection_num), 'name2' : stime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
    spatdf = pd.read_sql(spatquery, myoutputdb, params={'name' : options['cityintid'], 'name2' : startTime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
    tracks = df['track_id'].unique()
    recs = 0
    updaterows = []
    insertrows = []
    while (ftime < etime):
        results = []
        pool = mp.Pool(40)
        for tid in tracks:
            results.append(pool.apply_async(processTrack, args=(tid, df, spatdf, geoClusters, finalCentroidName, finalCentroids, tracksH)))
            #updaterows, insertrows, nt = processTrack(tid, df, spatdf, geoClusters, finalCentroidName, finalCentroids, tracksH, mycursor, myoutputdb)
            #tracksH[tid] = nt
        pool.close()
        pool.join()
        for i in range(0, len(tracks)):
            if (results[i].ready()):
                updaterows.extend(results[i].get()[0])
                insertrows.extend(results[i].get()[1])
                tracksH[tracks[i]] = results[i].get()[2]
        insertrows.sort(key=lambda x: x[0])
        mysqlws = time.time()
        mycursor.executemany(updatequery, updaterows)
        mysqlwb = time.time()
        mycursor.executemany(writequery, insertrows)
        mysqlwe = time.time()
        end = time.time()
        print ("Finished processing 1 minute worth of records",end-start, mysqlwb-mysqlws, len(updaterows), mysqlwe-mysqlwb, len(insertrows))
        updaterows.clear()
        insertrows.clear()
        start = end
        stime = ftime + step
        ftime = stime + minstep
        df = pd.read_sql(querylimit, myinputdb, params={'name' : str(intersection_num), 'name2' : stime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
        tracks = df['track_id'].unique()
        sdf = pd.read_sql(query, myoutputdb, params={'name': stime.strftime("%Y-%m-%d %H:%M:%S"), 'name2' : ftime.strftime("%Y-%m-%d %H:%M:%S"), 'name3' : options['cityintid']})
        spatrows = convertToSPaT(sdf, stime)
        mycursor.executemany(writespat, spatrows)
        myoutputdb.commit()
        startTime = ftime - datetime.timedelta(minutes=10)
        spatdf = pd.read_sql(spatquery, myoutputdb, params={'name' : options['cityintid'], 'name2' : startTime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
    mycursor.close()
    myoutputdb.close()

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

def findMatchingCentroid(track, geoclusters, centroid_list, centroids):
    x = track.x
    y = track.arr
    xlen = x.size
    cset = set(centroid_list)
    clist = []
    i = 0
    while (i < xlen-1):
        v, h = locatePointOnGrid(x[i], y[i])
        clist.extend(geoclusters[v][h])
        i = i+1
    cset = set(clist)
    found = 0
    mindist = 500
    mincent = ''
    centroid_track = centroids[0]
    for cent in cset:
        idx = centroid_list.index(cent)
        centroid_track = centroids[idx]
        isPed = 0
        if (track.phasetype == centroid_track.phasetype):
            isRelated, dist = closeToMembers(track, centroid_track, isPed)
            if (isRelated):
                found = 1
                if (dist < mindist):
                    mindist = dist
                    mincent = cent
    if (found == 0):
        print ('match_result', track.tid, 'anomalous')
        track.isAnomalous = 1
        track.clus ="anom"
    else:
        print ('match_result', track.tid, mincent[6:], track.phasetype, mincent)
        track.clus = mincent[6:]

    return centroid_track
if __name__=="__main__":
        main()
