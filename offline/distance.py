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
from sklearn.cluster import KMeans
from ast import literal_eval
from hash import HashTable
import time
from scipy.sparse import dok_matrix
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
from sqlalchemy import create_engine

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

distancedict = {}
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
        self.phasetype = 0 #major, minor, leftmajor, leftminor
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
        self.clus = ''
        self.iscent = 0

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
gIndex = 0

np.seterr(all='raise')
tracksH = HashTable(50000)
lengthH = HashTable(50000)

def swap (x1, x2):
    tmp = x1
    x1 = x2
    x2 = tmp
    return x1, x2

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

    deletequery='delete from CentroidsTable where intersection_id=name1 and camera_id=name2;'
    deletequery=deletequery.replace("name1", options['cityintid'])
    deletequery=deletequery.replace("name2", options['cameraid'])
    writequery='insert into CentroidsTable (track_id, x, y, class, camera_id, intersection_id, phase, name) values (%s, %s, %s, %s, %s, %s, %s, %s)'
    querylimit = 'select time,frame_id,track_id,center_x,center_y,class,timestamp,intersection_id,nearmiss,signature from OnlineTrackInfo where timestamp between %(name2)s and %(name3)s and intersection_id=%(name)s;'
    spatquery = 'select timestamp,hexphase,cycle from OnlineSPaT where timestamp between %(name2)s and %(name3)s and intersection_id=%(name)s;'

    start = time.time()
    cindex = 4
    rangex = 1280
    rangey = 960
    maxtrackid = -10000
    stime = pd.to_datetime(options['start'])
    ftime = stime + datetime.timedelta(minutes=int(options['collecttime']))
    etime = pd.to_datetime(options['end'])-datetime.timedelta(minutes=int(options['collecttime']))
    ctime = stime + datetime.timedelta(minutes=1) #replace with now
    interval = datetime.timedelta(minutes=int(options['clusterinterval']))
    endctime = ctime + interval
    while (ctime < etime):
        #sleep for (ftime - ctime)
        while (ctime < ftime):
            #wait for a minute
            ctime = ctime + datetime.timedelta(minutes=1)
        df = pd.read_sql(querylimit, myinputdb, params={'name' : options['cameraid'], 'name2' : stime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
        spatdf = pd.read_sql(spatquery, myoutputdb, params={'name' : options['cityintid'], 'name2' : stime.strftime("%Y-%m-%d %H:%M:%S.%f"), 'name3' : ftime.strftime("%Y-%m-%d %H:%M:%S.%f")})
        #df = pd.read_csv(options['combo'])
        tracks = df['track_id'].unique()
        objects = options['objects'].split(',')
        phases = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        gtracks = [[[] for j in range(len(phases))] for i in range(len(objects))]
        ntrackobj = np.zeros(shape=(len(objects),len(phases)))
    
        p=0
        excluded_car_tracks = []
        excluded_ped_tracks = []
        track_with_two_points = []
        blobs = []
        track_with_length_greater_than_1350 = []
        track_with_length_greater_than_2000 = []
        high_speed_ped_track_200 = []
        low_speed_veh_track_150 = []
        track_with_length_smaller_than_200 = []
        track_with_length_smaller_than_150_after_zoning = []
        high_speed_ped_track_130 = []
        low_speed_veh_track_100 = []
        tweezer_track = []
        tracklist = []
        lenlist = []
        speedlist = []
        isvehlist = []
        car_bb=[]
        bus_bb=[]
        motorbike_bb=[]
        ped_bb=[]
        car_ratio=[]
        bus_ratio=[]
        motorbike_ratio=[]
        ped_ratio=[]
        new_ped=[]
        new_motorbike=[]
        signal_violate=[]
        track_shape=[]
        track_lane_change=[]
        weird_shape=[]
    
        nummixedphases = 0
        pedindex = -1
        if ('pedestrian' in objects):
            pedindex = objects.index('pedestrian')
        maxpathlength = float(options['pathlength'])
        for tid in tracks:
            ndf = df.loc[df['track_id'] == tid]
            t = ndf.iloc[:,0].values
            f = ndf.iloc[:,1].values
            nx = ndf.iloc[:,3].values
            ny = ndf.iloc[:,4].values
            c = ndf.iloc[:,5].values
            ts = ndf.iloc[:,6].values
            iid = ndf.iloc[:,7].values
            #nsx = ndf.iloc[:,8].values
            #nsy = ndf.iloc[:,9].values
            #rs = ndf.iloc[:,11].values #phases
            #ps = ndf.iloc[:,11].values #phases
            #cycle = ndf.iloc[:,12].values #cycle
            si = 0
            se = np.size(nx)
            r = np.size(nx)
            exactmatch = spatdf.iloc[getSPaTIndex(spatdf, ts[r-2]+(f[r-2]%10)*np.timedelta64(100,'ms'))]
            rs = exactmatch['hexphase']
            cyclenum=exactmatch['cycle']
            phasetype = rs

            pathLength = getPathLength(tid, nx, ny)
            if (pathLength < maxpathlength):
                #print ("Ignored track of length < 159:", tid)
                track_with_length_smaller_than_200.append(tid)
                continue
            cindex = objects.index(c[0])
            #if (c[0] == 'pedestrian'):
            if (cindex == pedindex):
                excl, si, se = getRelevantCoordinates(pedpolygon, nx, ny, cindex==pedindex, tid)
                if (excl):
                    excluded_ped_tracks.append(tid)
                    continue
                x = nx[si:se]
                y = ny[si:se]
                t = t[si:se]
                f = f[si:se]
                ts = ts[si:se]
            else:
                excl, si, se = getRelevantCoordinates(polygon, nx, ny, cindex==pedindex, tid)
                if (excl):
                    excluded_car_tracks.append(tid)
                    continue
                x = nx[si:se]
                y = ny[si:se]
                t = t[si:se]
                f = f[si:se]
                ts = ts[si:se]
            xlen = x.size
            if (xlen <= 2):
                continue
            pathLength = getPathLength(tid, x, y)
            X = np.array([[x, y] for x, y in zip(x, y)])
            e = 1
            mask = rdp(X, epsilon=e, return_mask=True)
            t = np.ma.compressed(np.ma.masked_where(mask==False, t))
            f = np.ma.compressed(np.ma.masked_where(mask==False, f))
            x = np.ma.compressed(np.ma.masked_where(mask==False, x))
            y = np.ma.compressed(np.ma.masked_where(mask==False, y))
            ts = np.ma.compressed(np.ma.masked_where(mask==False, ts))
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

            pathLength = getPathLength(tid, x, y)
            lengthH[tid] = pathLength
            maxpathlength = float(options['pathlength'])
            if (pathLength < maxpathlength):
                continue
    
            rbit, rbit2 = getGreenBit(x,y)
            pindex = rbit - 1
            if (rbit == 0):
                pindex = int(8 + rbit2/2 - 1)
            if (pindex < 0):
                continue
            rrange = rangex
            nt = Track(tid, x, y, rrange)
            nt.x = x
            nt.arr = y
            nt.start = si
            nt.end = se
            nt.classtype = c[0]
            nt.phasetype = phasetype
            nt.cycle = cyclenum
            nt.len = pathLength
            nt.iid = iid[0]
            sx, sy, spindex, check = setSpeed(x, y, t, nt)
            maxSpeed = getMaxSpeed(tid, sx, sy)
            if (cindex == pedindex and maxSpeed > 130):
                continue
            elif (cindex != pedindex and maxSpeed < 100):
                continue
            gtracks[cindex][pindex].append(tid)
            ntrackobj[cindex][pindex] = ntrackobj[cindex][pindex] + 1
            #Create interpolated versions of each track
    
            tracksH[tid] = nt
            if (tid > maxtrackid):
                maxtrackid = tid
    
        print ("mixed phases:", nummixedphases)
        end = time.time()
        print ("Finished reading tracks file and initializing data structures", end - start)
        print ('excluded_car_tracks', len(excluded_car_tracks), excluded_car_tracks)
        print ('excluded_ped_tracks', len(excluded_ped_tracks), excluded_ped_tracks)
        print ('track_with_two_points', len(track_with_two_points), track_with_two_points)
        print ('blobs', len(blobs), blobs)
        print ('track_with_length_smaller_than_159', len(track_with_length_smaller_than_200), track_with_length_smaller_than_200)
        print ('track_with_length_smaller_than_150_after_zoning', len(track_with_length_smaller_than_150_after_zoning), track_with_length_smaller_than_150_after_zoning)
        print ('tweezer_track', len(tweezer_track), tweezer_track)
        print ('high_speed_ped_track_130', len(high_speed_ped_track_130), high_speed_ped_track_130)
        print ('low_speed_veh_track_100', len(low_speed_veh_track_100), low_speed_veh_track_100)
        print ('signal_violate', len(signal_violate), signal_violate)
        print ('weird_shape', len(weird_shape), weird_shape)
        print ('track_lane_change', len(track_lane_change), track_lane_change)
        print ('track_shape', len(track_shape), track_shape)
        startIndex = 1
        lobj = len(objects)
        #lobj = 1
        odf = pd.DataFrame()
        fid=0
        tid=0
        #change later 1->0 in range
        prefix = objects.copy()
        lph = len(phases)
    
        #fftracks = {}
        finalCentroids = []
        finalCentroidName = []
        tidlist = []
        xlist = []
        ylist = []
        classlist = []
        iidlist = []
        phaselist = []
        namelist = []
        centroid_list = []
        f = open("car_track_clusters.txt", "w")
        for i in range(0,lobj):
            if (objects[i] == 'car'):
                f.write('12\n')
            for j in range(0,lph):
                aHash = HashTable(1000)
                anotracksthreshold = 0
                if (ntrackobj[i][j] > 1000):
                    anotracksthreshold = int(ntrackobj[i][j] * 0.005)
                print ('Begin clustering', objects[i], 'for phase', phases[j])
                ptracks = gtracks[i][j]
                lptracks = len(ptracks)
                if (lptracks == 0):
                    continue
                A = []
                X = []
                F = []
                L = []
                for p in range(0,lptracks):
                    xd = []
                    d = []
                    fd = []
                    ld = []
                    A.append(d)
                    X.append(xd)
                    F.append(fd)
                    L.append(ld)
                    trackp = ptracks[p]
                    for q in range(0,lptracks):
                        trackq = ptracks[q]
                        if (trackp == trackq):
                            #d.append(0)
                            d.append(1)
                            xd.append(0)
                            fd.append(0)
                            ld.append(0)
                        else:
                            if ((trackp, trackq) in distancedict):
                                (tdist, firstdiff, lastdiff) = distancedict[(trackp, trackq)]
                            else:
                                tdist, maxdiff, dist, lenratio, firstdiff, lastdiff = distance(tracksH[trackp], tracksH[trackq], 200)
                                distancedict[(trackp, trackq)] = (tdist, firstdiff, lastdiff)
                                distancedict[(trackq, trackp)] = (tdist, firstdiff, lastdiff)
                            if (tdist > 25 or firstdiff > 15 or lastdiff > 15):
                                d.append(0)
                            else:
                                d.append(1) 
                            xd.append(tdist)
                            fd.append(firstdiff)
                            ld.append(lastdiff)
                for key,val in distancedict.items():
                    print (key, '->', val)
                #print(A)
                #print(X)
                #print(F)
                #print(L)
                D = np.diag(np.sum(A, axis=1))
                # graph laplacian
                L = D-A
                # eigenvalues and eigenvectors
                vals, vecs = np.linalg.eig(L)
                # sort these based on the eigenvalues
                vecs = vecs[:,np.argsort(vals)].real
                vals = vals[np.argsort(vals)].real
                print(vals)

                #limiteigen = 0.1 * vals[vals.size-1]
                #count = len([eigen for eigen in vals if eigen < limiteigen])
                count = 1
                if (vals.size > 1):
                    if (phases[j] % 2 == 0):
                        limiteigen = 0.05 * vals[vals.size-1]
                        if (limiteigen > 0):
                            count = len([eigen for eigen in vals if eigen < limiteigen])
                        else:
                            count = vals.size
                    else:
                        limiteigen = 0.1 * vals[vals.size-1]
                        if (limiteigen > 0):
                            count = len([eigen for eigen in vals if eigen < limiteigen])
                        else:
                            count = vals.size
                        #count = np.argmax(np.diff(vals)) + 1
                if (count > 1):
                    # kmeans on first three vectors with nonzero eigenvalues
                    kmeans = KMeans(n_clusters=count)
                    kmeans.fit(vecs[:,0:count-1])
                    colors = kmeans.labels_

                    findCentroidsOfMultipleClusters(ptracks, colors, phases[j], centroid_list)
                    print("Clusters:", count, colors)
                    for p in range(0,lptracks):
                        nt = tracksH[ptracks[p]]
                        nt.clus = objects[i] + str(phases[j]) + '-' + str(colors[p])
                    print(ptracks)
                    if (objects[i] == 'car'):
                        ptracksarray = np.asarray(ptracks)
                        numtracks=0
                        s=0
                        totaltracks = len(ptracks)
                        while (numtracks<totaltracks):
                            subarray = ptracksarray[np.where(colors==s)[0]]
                            f.write(str(len(subarray))+' ')
                            f.write(objects[i] + str(phases[j])+'-' + str(s) + ' ')
                            f.write(str(subarray).strip('[]')+'\n')
                            s = s+1
                            numtracks = numtracks + len(subarray)
                else:
                    print("There is just one cluster")
                    print(ptracks)
                    if (objects[i] == 'car'):
                        f.write(str(lptracks)+' ')
                        f.write(objects[i] + str(phases[j])+' ')
                        f.write(str(ptracks).strip('[]')+'\n')
                    findCentroid(ptracks, phases[j], centroid_list)
                    for p in range(0,lptracks):
                        nt = tracksH[ptracks[p]]
                        nt.clus = objects[i] + str(phases[j])
        f.close()
        writeToCentroidsTable(centroid_list, engine)
        #wait until ctime reaches endof interval time
        while (ctime < endctime):
            #wait for a minute
            ctime = ctime + datetime.timedelta(minutes=1)

def writeToCentroidsTable(centroid_list, engine):
    tidlist = []
    xlist = []
    ylist = []
    classlist = []
    cameraidlist = []
    iidlist = []
    phaselist = []
    namelist = []
    clen = len(centroid_list)
    for c in centroid_list:
        xlen  = c.x.size
        tidlist.extend([c.tid for k in range(xlen)])
        xlist.extend(c.x.tolist())
        ylist.extend(c.arr.tolist())
        classlist.extend([c.classtype for k in range(xlen)])
        cameraidlist.extend([options['cameraid'] for k in range(xlen)])
        iidlist.extend([options['cityintid'] for k in range(xlen)])
        phaselist.extend([c.phasetype for k in range(xlen)])
        namelist.extend([c.clus for k in range(xlen)])
    centdf = pd.DataFrame()
    centdf['track_id'] = tidlist
    centdf['x'] = xlist
    centdf['y'] = ylist
    centdf['class'] = classlist
    centdf['camera_id'] = cameraidlist
    centdf['intersection_id'] = iidlist
    centdf['phase'] = phaselist
    centdf['name'] = namelist
    centdf.to_sql('CentroidsTable', con=engine, if_exists='append', index = False)
    return

def findCentroidsOfMultipleClusters(ptracks, colors, phase, centroid_list):
    ptracksarray = np.asarray(ptracks)
    numtracks=0
    i=0
    totaltracks = len(ptracks)
    while (numtracks<totaltracks):
        subarray = ptracksarray[np.where(colors==i)[0]]
        arraylen = len(subarray)
        numtracks = numtracks + arraylen
        if (arraylen > 1):
            findCentroid(subarray.tolist(), phase, centroid_list)
        i = i+1

def findCentroid(ptracks, phase, centroid_list):
    mindist = 10000
    mintrack = -1
    maxlen = 0
    slen = len(ptracks)
    if (slen == 0):
        return 

    #Compute centroid among 10% of members chosen randomly
    if (slen == 1 or slen == 2):
        centroid_list.append(tracksH[ptracks[0]])
        return 

    for track in ptracks:
        nt1 = tracksH[track]
        dist = 0
        for otrack in ptracks:
            nt2 = tracksH[otrack]
            if (nt1.tid != nt2.tid):
                (val, firstdiff, lastdiff) = distancedict[(track, otrack)]
                dist = dist + val
        if (dist < mindist):
            mindist = dist
            mintrack = track
    centroid_list.append(tracksH[mintrack])
    return

def setSpeed(x, y, ts, nt):
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
        index = np.argmax(sp2 > 100)
        sx = sp[:, [0]].reshape(x.size-1,)
        sy = sp[:, [1]].reshape(x.size-1,)
        slistx.extend(np.round(sx, 2).tolist())
        slisty.extend(np.round(sy, 2).tolist())
    return slistx, slisty, index, check

def getSPaTIndex(spatdf, ts):
    index = 0
    mask = spatdf['timestamp'] > ts
    pos = np.flatnonzero(mask)
    if (pos.size==0):
        index = len(spatdf.index)-1
    elif (pos[0] > 0):
        index = pos[0]-1
    return index

def getDistanceOfPointFromLine(x1, y1, x2, y2, x0, y0):
    dist = np.absolute((y2-y1) * x0 - (x2 - x1) * y0 + x2*y1 - y2*x1)/(np.sqrt(np.add(np.square(y2-y1), np.square(x2-x1))))
    return dist

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

def distanceNew(nt, nt1, boxDim):
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
            #distance, path = fastdtw (ntpair, nt1pair, dist=euclidean)
    lenratio = 0
    distance = 0
    maxdiff = 100
    angled=180
    firstdiff = 0
    lastdiff = 0
    path = []
    if (distance < 35000):
        minsize = min(nt.x.size, nt1.x.size)
        maxsize = max(nt.x.size, nt1.x.size)
        for i in range(0, minsize-1):
            path.append((i, i))
        for j in range(minsize, maxsize-1):
            if (nt.x.size > nt1.x.size):
                path.append((j, i))
            else:
                path.append((i, j))
        first = 1
        prevpair = (0,0)
        longpath = 0
        area = 0
        z1=0
        z=0
        first = 1
        distlist = []
        alist = []
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
        if (zlist[0] + z1list[0] > 0):
            firstdiff = alist[0]/(zlist[0] + z1list[0])
        if (zlist[len(zlist)-1] + z1list[len(z1list)-1] > 0):
            lastdiff = alist[len(alist)-1]/(zlist[len(zlist)-1] + z1list[len(z1list)-1])
        if ((z + z1) > 0):
            tdist = area/(z + z1)
            lenratio = (z + z1)/(z + z1 + el + el1)
        longlength = z
        if (longpath == 1):
            longlength = z1

    return tdist+maxdiff, maxdiff, distance, lenratio, firstdiff, lastdiff

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
    firstdiff = 0
    lastdiff = 0
    if (distance < 35000):
        first = 1
        prevpair = (0,0)
        longpath = 0
        area = 0
        z1=0
        z=0
        first = 1
        distlist = []
        alist = []
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
        if (zlist[0] + z1list[0] > 0):
            firstdiff = alist[0]/(zlist[0] + z1list[0])
        if (zlist[len(zlist)-1] + z1list[len(z1list)-1] > 0):
            lastdiff = alist[len(alist)-1]/(zlist[len(zlist)-1] + z1list[len(z1list)-1])
        if ((z + z1) > 0):
            tdist = area/(z + z1)
            lenratio = (z + z1)/(z + z1 + el + el1)
        longlength = z
        if (longpath == 1):
            longlength = z1

    return tdist+maxdiff, maxdiff, distance, lenratio, firstdiff, lastdiff

def mergeDistance(nt, nt1):
    tdist = 100
    ntpair = [[x, y] for x, y in zip(nt.x, nt.arr)]
    nt1pair = [[x, y] for x, y in zip(nt1.x, nt1.arr)]
    distance = 100000
    [Ax, Ay] = [nt.x[0], nt.arr[0]]
    [Bx, By] = [nt1.x[0], nt1.arr[0]]
    distance, path = fastdtw (ntpair, nt1pair, dist=euclidean)
    lenratio = 1
    maxdiff = 100
    angled=180
    if (distance < 75000):
        #Compute cosine of the angle between the start points:
        #Compute overlap section
        first = 1
        prevpair = (0,0)
        longpath = 0
        #Compute sum of area of overlapping section
        area = 0
        z1=0
        z=0
        distlist = []
        alist = []
        zlist = []
        z1list = []
        el = 0
        el1 = 0
        #for pair in xindexpairs:
        for pair in path:
            if (first):
                first = 0
            elif ((pair [0] != prevpair[0]) and (pair[1] != prevpair[1])):
                #A: (nt.x[pair[0]], nt.arr[pair[0])
                #B = (nt1.x[pair[1]], nt1.arr[pair[1]])
                #C = (nt1.x[prevpair[1]], nt1.arr[prevpair[1]])
                [Ax, Ay] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Bx, By] = [nt1.x[pair[1]], nt1.arr[pair[1]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                a = 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                z1length = math.sqrt((Bx - Cx)**2 + (By - Cy)**2)
                z1 = z1 + z1length
                #z1 = z1 + math.sqrt((nt1.x[pair[1]] - nt1.x[prevpair[1]])**2 + (nt1.arr[pair[1]]-nt1.arr[prevpair[1]])**2)
                #A: (nt.x[prevpair[0]], nt.arr[prevpair[0])
                #B = ((nt.x[pair[0]], nt.arr[pair[0])
                #C = (nt1.x[prevpair[1]], nt1.arr[prevpair[1]])
                [Ax, Ay] = [nt.x[prevpair[0]], nt.arr[prevpair[0]]]
                [Bx, By] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Cx, Cy] = [nt1.x[prevpair[1]], nt1.arr[prevpair[1]]]
                a = a + 0.5*np.absolute(Ax * (By - Cy) + Bx * (Cy - Ay) + Cx * (Ay - By))
                area = area + a
                zlength = math.sqrt((Ax - Bx)**2 + (Ay - By)**2)
                z = z + zlength
                #z = z + math.sqrt((nt.x[prevpair[0]]-nt.x[pair[0]])**2 + (nt.arr[prevpair[0]]-nt.arr[pair[0]])**2)
                [Ax, Ay] = [nt.x[pair[0]], nt.arr[pair[0]]]
                [Bx, By] = [nt1.x[pair[1]], nt1.arr[pair[1]]]
                distlist.append(math.sqrt((Ax-Bx)**2 + (Ay-By)**2))
                alist.append(a)
                z1list.append(z1length)
                zlist.append(zlength)
            prevpair = pair

        distarray = np.asarray(distlist)
        maxdiff = np.amax(distarray)
        #diff = distarray[distarray>10]
        if ((z + z1) > 0):
            tdist = area/(z + z1)
            lenratio = (z + z1)/(z + z1 + el + el1)
        longlength = z
        if (longpath == 1):
            longlength = z1

    #print ("Distance: ", nt.tid, nt1.tid, tdist, lenratio, distance, maxdiff, angled)
    #print ("Distance: ", nt1.tid, nt.tid, tdist)
    return tdist, maxdiff, distance, lenratio

def computeCentroid(trackseries, centroids, currentIndex, boxDim):
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
                val, maxdiff, dtwdist, ovlp = distance(nt1, nt2, boxDim)
                dist = dist + val
        if (dist < mindist):
            mindist = dist
            mintrack = track
    return tracksH[mintrack]

def KMeansAssignInit(cdf, centroids, ptracks, startIndex, distanceThreshold, overlapThreshold, boxDim):
    s = len(centroids)
    for i in range (0,s+1):
        d = []
        if (i < s):
            nt = centroids[i]
            for j in ptracks:
                ntref = tracksH[j]
                if (nt == ntref):
                    val = 0
                    ovlp = 1
                elif (ntref.isAnomalous == 0 and nt != ntref):
                    val, maxdiff, dtwdist, ovlp = distance(nt, ntref, boxDim)
                else:
                    val = 100
                if (val > distanceThreshold or ovlp < overlapThreshold):
                    val = 100
                d.append(val)
        else:
            for j in ptracks:
                d.append(99)
        cdf['distance_from_{}'.format(i)] = (
            d
            )
    centroid_distance_cols = ['distance_from_{}'.format(i) for i in range(0, s+1)]
    cdf['closest'] = cdf.loc[:, centroid_distance_cols].idxmin(axis=1)
    cdf['closest'] = cdf['closest'].map(lambda x: int(x.lstrip('distance_from_')))
    return cdf

def getCosineOfAngle(Axi, Ayi, Bxi, Byi):
    a = 0
    if ((Axi != 0 or Ayi != 0) and (Bxi != 0 or Byi != 0)):
        [nAxi, nAyi] = 1/np.sqrt(Axi*Axi + Ayi*Ayi) * np.array([Axi, Ayi])
        [nBxi, nByi] = 1/np.sqrt(Bxi*Bxi + Byi*Byi) * np.array([Bxi, Byi])
        a = nAxi*nBxi + nAyi*nByi
    return a

def getSweep(nt):
    ntsize = nt.x.size-1
    a = 0 #orthogonal
    [Axi, Ayi] = [nt.x[0] - nt.x[1], nt.arr[0] - nt.arr[1]]
    [Bxi, Byi] = [nt.x[0] - nt.x[ntsize], nt.arr[0] - nt.arr[ntsize]]
    a = getCosineOfAngle(Axi, Ayi, Bxi, Byi)
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
            #dist = np.absolute((y1[i+1]-y1[i]) * x0 - (x1[i+1] - x1[i]) * y0 + x1[i+1]*y1[i] - y1[i+1]*x1[i])/(np.sqrt(np.add(np.square(y1[i+1]-y1[i]), np.square(x1[i+1]-x1[i]))))
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

def KMeansAssign(cdf, centroids, ptracks, startIndex, distanceThreshold, overlapThreshold, boxDim):
    s = len(centroids)
    throughthresh = float(options['through'])
    for i in range (0,s+1):
        d = []
        if (i < s):
            nt = centroids[i]
            for j in ptracks:
                ntref = tracksH[j]
                if (nt == ntref):
                    val = 0
                    ovlp = 1
                #elif (ntref.isAnomalous == 0 and nt != ntref):
                else:
                    direction = checkDirection(nt, ntref)
                    difflen = np.absolute(nt.len-ntref.len)/np.maximum(nt.len, ntref.len)
                    val = -1
                    if (direction > throughthresh and difflen < 0.5):
                        val, maxdiff, dtwdist, ovlp = distance(nt, ntref, boxDim)
                    if (val == -1 or val >= 10):
                        val = 100
                #if (val == -1 and val > distanceThreshold or ovlp < overlapThreshold):
                d.append(val)
        else:
            for j in ptracks:
                d.append(99)
        cdf['distance_from_{}'.format(i)] = (
            d
            )
    centroid_distance_cols = ['distance_from_{}'.format(i) for i in range(0, s+1)]
    cdf['closest'] = cdf.loc[:, centroid_distance_cols].idxmin(axis=1)
    cdf['closest'] = cdf['closest'].map(lambda x: int(x.lstrip('distance_from_')))
    #cdf['color'] = cdf['closest'].map(lambda x: colmap[x])
    #find the centroid distance cols that are all ones
    #for each such index, create a new centroid
    return cdf

def KMeansIterate(cdf, centroids, ptracks, rrange, startIndex, distanceThreshold, overlapThreshold, boxDim):
    j = 0
    while True:
        #print ("iteration:", j)
        j = j+1
        closest_centroids = cdf['closest'].copy(deep=True)
        fname = 'cdf'+str(j)
        update(centroids, cdf, rrange, startIndex, boxDim)
        cdf = KMeansAssign(cdf, centroids, ptracks, startIndex, distanceThreshold, overlapThreshold, boxDim)
        slen = len(centroids)
        print ("Centroids")
        for i in range(0, slen):
            print (i, centroids[i].tid)
        if closest_centroids.equals(cdf['closest']):
            break
        if (j > 1):
            break

def update(centroids, df, rrange, startIndex, boxDim):
    s = len(centroids)
    for i in range (0, s):
        centroids[i] = computeCentroid((df[df['closest'] == i]['track_id']), centroids, i, boxDim)

def printClusterCenter(centroids, cid, rrange, cid_list, x_list, y_list):
    center = centroids[cid]
    y1 = center.arr
    for i in range (0, rrange):
        if (y1[i] > 0):
            cid_list.append(cid)
            x_list.append(i)
            y_list.append(y1[i])

def printRedJumpToFile(signal_violate):
    sv = pd.DataFrame()
    sv['track_id'] = signal_violate
    sv.to_csv('../DataFiles/redjump.csv', index=False)

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
    df_tracks_with_clusters['cycle'] = cdf['cycle']

    df_tracks_with_clusters.dropna(inplace=True)
    df_tracks_with_clusters.to_csv(filename,  float_format='%.f', mode=mode, header=False, index=False)

def KMeansInit(centroids, ptracks, rrange, isVeh, distanceThreshold, overlapThreshold, boxDim):
    print ("KMeans Init stats:")
    throughthresh = float(options['through'])
    lim = int(1.0 * len(ptracks))
    if (rrange == 1):
        lim = len(ptracks)
    listptracks = ptracks[0:lim]
    newlistptracks = []
    while (listptracks):
        newlistptracks.clear()
        nt1 = tracksH[listptracks[0]]
        close = 0
        for i in range (1, len(listptracks)):
            nt2 = tracksH[listptracks[i]]
            direction = checkDirection(nt1, nt2)
            difflen = np.absolute(lengthH[nt1.tid]-lengthH[nt2.tid])/np.maximum(lengthH[nt1.tid], lengthH[nt2.tid])
            if (direction > throughthresh and difflen < 0.5):
                val, maxdiff, dtwdist, ovlp = distance(nt1, nt2, boxDim)
                if (val < 10):
                    continue
                else:
                    newlistptracks.append(listptracks[i])
        centroids.append(nt1)
        listptracks = newlistptracks.copy()

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
        if (isCurved(biggerTrack) and lb > 250):
            threshold = 0.70
            ls = smallertrack.len
            tdistthresh = 60 * lb/ls
        else:
            threshold = 0.75 #getSweep(biggerTrack)
            tdistthresh = 60
    else:
        biggerTrack, smallertrack, status = isLengthDifferent(nt1, nt2)
        lb = biggerTrack.len
        if (isCurved(biggerTrack) and lb > 250):
            threshold = 0.70
            ls = smallertrack.len
            tdistthresh = 60 * lb/ls
        else:
            threshold = 0.97
            tdistthresh = 50

    return threshold, tdistthresh

def getTrackMatchingThreshold(nt1, nt2):
    #compare lengths of nt1 and nt2
    threshold = 1.0
    tdistthresh = 20
    throughthresh = float(options['through'])
    c1 = isCurved(nt1)
    c2 = isCurved(nt2)
    if (c1 == False and c2 == False):
        threshold = throughthresh
    elif (c1 == True and c2 == True):
        biggerTrack, smallertrack, status = isLengthDifferent(nt1, nt2)
        lb = biggerTrack.len
        if (isCurved(biggerTrack) and lb > 250):
            threshold = 0.70
            ls = smallertrack.len
            tdistthresh = 50 * lb/ls
        else:
            threshold = 0.75 #getSweep(biggerTrack)
            tdistthresh = 50
    else:
        biggerTrack, smallertrack, status = isLengthDifferent(nt1, nt2)
        lb = biggerTrack.len
        if (isCurved(biggerTrack) and lb > 250):
            threshold = 0.70
            ls = smallertrack.len
            tdistthresh = 50 * lb/ls
        else:
            threshold = 0.97
            tdistthresh = 50
    return threshold, tdistthresh

def findLongestPath(index, mclus, centroids):
    group = mclus[index]
    #find the longest track
    lenlist = []
    for member in group:
        lenlist.append(centroids[member].len)
    length = np.percentile(lenlist, 85)
    lenlistarray = np.asarray(lenlist)
    idx = (np.abs(lenlistarray - length)).argmin()
    return group[idx]

def isCloselyRelated(nt1, nt2, isPed):
    status = False
    direction = checkDirection(nt1, nt2)
    scale = max(nt1.scale, nt2.scale)
    tdist = 100
    if (isPed):
        if (direction > 0.82):
            #tdist, maxdiff, dtwdist, ovlp = mergeDistance(nt1, nt2)
            tdist, maxdiff, dtwdist, ovlp = distance(nt1, nt2, 1000)
            tdist = tdist * scale
            if (tdist < 50 and dtwdist < 10000):
                status = True
            print(nt1.tid, nt2.tid, tdist, maxdiff, dtwdist, ovlp, direction, status)
    else:
        threshold, tdistthresh = getTrackMatchingThreshold(nt1, nt2)
        if (direction > threshold):
            #check if the tracks start on the same lane
            isSameLane = isOnSameLane(nt1, nt2)
            if (isSameLane):
                tdist, maxdiff, dtwdist, ovlp = distance(nt1, nt2, 1000)
                #if (maxdiff < 50 and tdist < 15 and dtwdist < 12000):
                if (tdist < tdistthresh):
                    status = True
                elif (tdist < 50 and isCollinear(nt1, nt2)):
                    n1 = nt1.x.size
                    n2 = nt2.x.size
                    if (tdist < 25):
                        status = True
                    elif ((nt1.x[n1-1] < nt2.x[0] or nt2.x[n2-1] < nt1.x[0]) and strictly_increasing(nt1.x)):
                        status = True
                    elif ((nt1.x[n1-1] > nt2.x[0] or nt2.x[n2-1] > nt1.x[0]) and strictly_decreasing(nt1.x)):
                        status = True
                    elif (tdist < 45 and isCurved(nt2)):
                        status = True
                #elif (tdist < 25 and maxdiff < 25): #Assuming lane size is 50 pixels
                    #status = True
                print(nt1.tid, nt2.tid, tdist, maxdiff, dtwdist, ovlp, direction, isSameLane, status)
    return tdist, status

def isComparable(len1, len2):
    maxlen = max(len1, len2)
    minlen = min(len1, len2)
    diff = (maxlen - minlen)/maxlen
    status = True
    if (diff > 0.5):
        status = False
    return status

def closeToMembers(centi, centj, isPed):
    #find the longest track
    status = True
    direction = checkDirection(centi, centj)
    scale = 1

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
                    elif (tdist < 45 and isCurved(centj)):
                        status = True
                    else:
                        status = False
                print(centi.tid, centj.tid, tdist, maxdiff, dtwdist, ovlp, direction, status)
    return status

def printCentroidAndMembers(centroids, cdf):
    s = len(centroids)
    for i in range(0,s):
        print ("Centroid:", centroids[i].tid)
        ser = (cdf[cdf['closest'] == i]['track_id'])
        n = 0
        for j, track in ser.iteritems():
            nt1 = tracksH[track]
            print (j, nt1.tid)
            n = n + 1
        print ("Total tracks:",n)

def postProcessClusters(cdf, centroids):
    for i in range(len(centroids) - 1, -1, -1):
        element = centroids[i]
        s = cdf[cdf['closest'] == i]['track_id']
        if (s.size <= 5):
            centroids[i].isAnomalous = 1
            #print ("Replaced", i, -2)
            cdf['closest'].replace(i, -2, inplace=True)
            #del centroids[i]

def mergeCentroids(centroids, cdf, isPed, boxDim):
    mclus = HashTable(1000)
    hclus = HashTable(1000)
    cdist = HashTable(1000)
    cclus = HashTable(1000)
    s = len(centroids)
    idclus = 1
    for i in range(0,s):
        if (centroids[i].isAnomalous == 1): # or hclus[i] != None):
            continue
        print ("Up i for merging", i, centroids[i].tid, hclus[i])
        for j in range(i+1, s):
            if (centroids[j].isAnomalous == 1):
                continue
            centi = centroids[i]
            centj = centroids[j]
            if (hclus[i] != None):
                centi = centroids[findLongestPath(hclus[i], mclus, centroids)]
            if (hclus[j] != None):
                centj = centroids[findLongestPath(hclus[j], mclus, centroids)]
            print ("Up i, j for merging", i, j, centroids[i].tid, centroids[j].tid)
            if (centi == centj):
                print("Tracks are part of the same cluster")
                continue
            currentDistance, isRelated = isCloselyRelated(centi, centj, isPed)
            if (isRelated):
                print("Potential merge:", i, j)
                clus1 = hclus[i]
                clus2 = hclus[j]
                if (clus1 == None and clus2 == None):
                    hclus[i] = idclus
                    hclus[j] = idclus
                    mclus[idclus] = []
                    mclus[idclus].append(i)
                    mclus[idclus].append(j)
                    print ("Created a new cluster", idclus, "for", i, j)
                    idclus = idclus + 1
                elif (clus1 == None):
                    #if (closeToMembers(i, hclus[j], centroids, mclus, isPed, boxDim)):
                    if (closeToMembers(centi, centj, isPed)):
                        hclus[i] = hclus[j]
                        mclus[hclus[i]].append(i)
                        print ("Appended to existing cluster", hclus[i], i)
                elif (clus2 == None):
                    #if (closeToMembers(j, hclus[i], centroids, mclus, isPed, boxDim)):
                    if (closeToMembers(centi, centj, isPed)):
                        hclus[j] = hclus[i]
                        mclus[hclus[j]].append(j)
                        print ("Appended to existing cluster", hclus[j], j)
                elif (clus1 != clus2):
                    centi = centroids[findLongestPath(hclus[i], mclus, centroids)]
                    centj = centroids[findLongestPath(hclus[j], mclus, centroids)]
                    #if (closeToMembers(i, hclus[j], centroids, mclus, isPed, boxDim) and closeToMembers(j, hclus[i], centroids, mclus, isPed, boxDim)):
                    if (isCloselyRelated(centi, centj, isPed)):
                        for cent in mclus[hclus[i]]:
                            if (cent != i):
                                mclus[hclus[j]].append(cent)
                                hclus[cent] = hclus[j]
                        mclus[hclus[i]].clear()
                        mclus[hclus[i]] = []
                        print ("Merging clusters", hclus[i], hclus[j], i, j)
                        hclus[i] = hclus[j]
                        mclus[hclus[j]].append(i)
        if (hclus[i] == None):
            hclus[i] = idclus
            mclus[idclus] = []
            mclus[idclus].append(i)
            print ("Created a new cluster", idclus, "for", i)
            idclus = idclus + 1
        for k in range(0,s):
            track = centroids[k]
            if (hclus[k] != None):
                centk = centroids[findLongestPath(hclus[k], mclus, centroids)]
                tdist, maxdiff, dtwdist, ovlp = distance(track, centk, 1000)
                cdist[k] = tdist

    #Replace in cdf
    for i in range(1,idclus):
        #print (i, mclus[i])
        for j in range(0, len(mclus[i])):
            #print ("Replaced", mclus[i][j], mclus[i][0])
            cdf['closest'].replace(mclus[i][j], mclus[i][0], inplace=True)

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

#    if (isPed):
#        i = si
#        while (i < se):
#            p = Point(x[i], y[i])
#            if (pedexclusion.contains(p)):
#                excl = 1
#                print ("Found ped track in exclusion zone", tid)
#                break
#            i = i + 1

    return excl, si, se

def isStaticTrack(x, y, tid):
    status = False
    if (x.size > 0):
        a = np.asarray([[x, y] for x, y in zip(x, y)])
        mid = 0 #int(x.size/2)
        b = np.asarray([[x[mid], y[mid]]] * x.size)
        d = np.sqrt(np.einsum('ij,ij->i', (a-b), (a-b))) #https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
        df = np.diff(d)
        decreases = df[df < 0]
        if (decreases.size > 0 and np.mean(d) < 30):
            #print ("Blob:", tid, np.mean(d), decreases.size, x.size)
            status = True
    return status

def getMaxSpeed(tid, sx, sy):
    d = 0
    if (len(sx) > 0):
        d = np.max(np.sqrt(np.add(np.square(sx), np.square(sy))))
    return d

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

def getInitialSpeed(sx, sy):
    speed = 0
    i = 0;
    while (i < sx.size and sx[i] == 0 and sy[i] == 0):
        i = i + 1

    speed = np.sqrt(np.add(np.square(sx[i]), np.square(sy[i])))

    return speed

def getPathLength(tid, x, y):
    length = 0
    if (x.size > 0):
        aa = np.asarray([[x, y] for x, y in zip(x, y)])
        a = aa[:x.size-1]
        b = np.copy(aa[1:])
        d = np.sqrt(np.einsum('ij,ij->i', (a-b), (a-b))) #https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
        length = np.sum(d)
    return length

def trackChangesLane(x, y, start, end, rbit):
    status = False

    #[Ax, Ay] = [x[x.size-1]-x[0], y[y.size-1]-y[0]]
    #if (rbit == 2 or rbit == 6):
        #a = np.absolute(getCosineOfAngle(Ax, Ay, dir_6_x - dir_2_x, dir_6_y - dir_2_y))
        ##if (a > 0.975 and a < 0.99):
        ##if (a > 0.995 and a < 0.999):
        #if (a > start and a < end):
            #status = True
    #elif (rbit == 4):
        #a = np.absolute(getCosineOfAngle(Ax, Ay, dir_8_x_8 - dir_4_x_8, dir_8_y_8 - dir_4_y_8))
        #if (a > start and a < end):
            #status = True
    #elif (rbit == 8):
        #a = np.absolute(getCosineOfAngle(Ax, Ay, dir_4_x_4 - dir_8_x_4, dir_4_y_4 - dir_8_y_4))
        #if (a > start and a < end):
            #status = True

    return status

if __name__=="__main__":
    main()
