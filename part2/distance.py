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
from scipy.cluster.hierarchy import fclusterdata
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
import copy

#@dataclass
class Track:
    def __init__(self, tid, xstart, xend, ystart, yend, rrange, detx, p, c):
        self.tid = tid        #track id
        self.xstart = xstart  #leftmost x coordinate of path
        self.xend = xend      #rightmost x coordinate of path
        self.ystart = ystart  #bottommost y coordinate of path
        self.yend = yend      #topmost y coordinate of path
        self.detx = detx      #vertical trajectory
        self.rrange = rrange
        self.arr = np.empty(shape=rrange)
        self.cid = -1
        self.phasetype = p #major, minor, leftmajor, leftminor
        self.classtype = c #other types are bus, truck, ped
        self.cycles = []      # cycle numbers for which the vehicle is present in the intersection

del_file = '../DataFiles/cluster.csv'
if (os.path.isfile(del_file)):
    os.remove(del_file)
#file_to_rem = pathlib.Path('cluster.csv')
#file_to_rem.unlink()

df = pd.read_csv('../DataFiles/combo.csv')
#df = pd.read_csv('combo.csv',names=['frame_id', 'track_id', 'x', 'y', 'class', 'timestamp','phase'], index_col='timestamp')
#columns = ['frame_id', 'track_id', 'x', 'y', 'class', 'timestamp','phase']

spatdf=pd.read_csv('../DataFiles/SPAT.csv',names=['timestamp', 'hexphase', 'binphase', 'twobitphase', 'cyclenumber'], index_col='timestamp')

np.seterr(all='raise')
fps = 20
xvertBox = 15
maxPedSpeed = 30 #pixels per second
boxOverlapThreshold = 100
tracksH = HashTable(50000)
aHash = HashTable(50000)

#Find the two array indexes between which val lies
def setindexes(val, array):
    r = np.size(array)
    index1 = 0
    index2 = 0
    for i in range(0, r-1):
        if (array[0] < array[r-1] and val <= array[i+1] and val >= array[i]):
            index1 = i
            index2 = i+1
            break
        elif (array[0] > array[r-1] and val >= array[i+1] and val <= array[i]):
            index1 = i
            index2 = i+1
            break
    return index1, index2
        
def linearInterpolation(val, x, y, detx):
    if (detx):
        index1, index2 = setindexes(val, y)
        xval = x[index1]
        if ((y[index2] - y[index1])!=0):
            xval = x[index1] + (x[index2] - x[index1])/(y[index2] - y[index1]) * (val - y[index1])
        return xval
    else:
        index1, index2 = setindexes(val, x)
        yval = y[index1]
        if ((x[index2] - x[index1])!=0):
            yval = y[index1] + (y[index2] - y[index1])/(x[index2] - x[index1]) * (val - x[index1])
        return yval

def getPhase(x, y, t, phases):
    p = t[0] #Check
    if (p not in phases):
        phases.append(p)

    return phases.index(p)

def printTrajectoryFileFor2008Paper(filename, ptracks):
    s = len(ptracks)
    fh = open(filename, 'w') 

    fh.write("2\n")
    fh.write(str(s))
    fh.write("\n")

    for i in range(0,s-1):
        carray = ptracks[i]
        asize = len(carray)/2
        fh.write(str(i))
        fh.write(" ")
        fh.write(str(int(asize)))
        fh.write(" ")
        for item in carray:
            fh.write("%s " % item)
        fh.write("\n");

    fh.close()

nhclus =HashTable(50000)
hclus =HashTable(50000)
dclus = HashTable(50000)
def closeToMembers(tid, refclus, df_, mclus, threshold):
    group = mclus[refclus]
    #print ("closeToMembers", group)
    for member in group:
        #print ("comparing", tid, member, df_.at[tid, member])
        if(df_.at[tid, member] > threshold):
            return 0
    #print ("closeToMembers successful")
    return 1

def distance(nt, nt1):
    maximumdistance = 450000.0
    dist = 1.0
    y1 = nt.arr
    y2 = nt1.arr
    if (np.size(y1) == np.size(y2)):
        process = 1
        #if ((np.absolute(np.subtract(nt1.xstart,nt.xstart)) < boxOverlapThreshold) and
                #(np.absolute(np.subtract(nt1.xend,nt.xend)) < boxOverlapThreshold)):
            #if ((np.absolute(np.subtract(nt1.ystart,nt.ystart)) < boxOverlapThreshold) or
                #(np.absolute(np.subtract(nt1.yend,nt.yend)) < boxOverlapThreshold)):
                #process = 1
        #elif ((np.absolute(np.subtract(nt1.ystart,nt.ystart)) < boxOverlapThreshold) and
                #(np.absolute(np.subtract(nt1.yend,nt.yend)) < boxOverlapThreshold)):
            #if ((np.absolute(np.subtract(nt1.xstart,nt.xstart)) < boxOverlapThreshold) or
                #(np.absolute(np.subtract(nt1.xend,nt.xend)) < boxOverlapThreshold)):
                #process = 1
        y12 = np.absolute(np.subtract(y1, y2))
        dist = np.sum(y12)/maximumdistance

    return dist

def computeCentroid(trackseries, centroids, rrange):
    p = 0
    c = 0
    for i, track in trackseries.iteritems():
        nt1 = tracksH[track]
        p = nt1.phasetype
        c = nt1.classtype
        break
    nt = Track(-1, 0,0,0,0,rrange,0, p, c)
    nt.tid = -1 #reserved for centroid tracks
    numtracks = 0
    centroidVal =  np.zeros(shape=rrange)
    xStart =  0
    xEnd =  0
    yStart =  0
    yEnd =  0
    for i, track in trackseries.iteritems():
        nt1 = tracksH[track]
        centroidVal = np.add(centroidVal,nt1.arr)
        xStart = xStart + nt1.xstart
        xEnd = xEnd + nt1.xend
        yStart = yStart + nt1.ystart
        yEnd = yEnd + nt1.yend
        numtracks = numtracks + 1
    if (numtracks > 0):
        centroidVal = np.divide(centroidVal,numtracks)
        nt.xstart = xStart/numtracks
        nt.xend = xEnd/numtracks
        nt.ystart = yStart/numtracks
        nt.yend = yEnd/numtracks
    nt.arr = centroidVal
    return nt

def KMeansAssign(cdf, centroids, ptracks, startIndex):
    s = len(centroids)
    for i in range (0,s):
        d = []
        nt = centroids[i]
        for j in ptracks:
            ntref = tracksH[j]
            d.append(distance(nt, ntref))
        #cdf['distance_from_{}'.format(i+startIndex)] = (
        cdf['distance_from_{}'.format(i)] = (
            d
            )
    #centroid_distance_cols = ['distance_from_{}'.format(i+startIndex) for i in range(0, s)]
    centroid_distance_cols = ['distance_from_{}'.format(i) for i in range(0, s)]
    cdf['closest'] = cdf.loc[:, centroid_distance_cols].idxmin(axis=1)
    cdf['closest'] = cdf['closest'].map(lambda x: int(x.lstrip('distance_from_')))
    #cdf['color'] = cdf['closest'].map(lambda x: colmap[x])
    return cdf

def KMeansIterate(cdf, centroids, ptracks, rrange, startIndex):
    j = 0
    while True:
        print ("iteration:", j)
        j = j+1
        printClusters(cdf)
        closest_centroids = cdf['closest'].copy(deep=True)
        update(centroids, cdf, rrange, startIndex)
        cdf = KMeansAssign(cdf, centroids, ptracks, startIndex)
        if closest_centroids.equals(cdf['closest']):
            break
        if (j > 30):
            break

def update(centroids, df, rrange, startIndex):
    s = len(centroids)
    for i in range (0, s):
        centroids[i] = computeCentroid((df[df['closest'] == i]['track_id']), centroids, rrange)
        #centroids[i] = computeCentroid((df[df['closest'] == i+startIndex]['track_id']), centroids, rrange)
        #for i, row in df[df['closest'] == i]['x'].iteritems():
            #print (i, row)

def printClusterCenter(centroids, cid, rrange, cid_list, x_list, y_list):
    center = centroids[cid]
    y1 = center.arr
    for i in range (0, rrange):
        if (y1[i] > 0):
            cid_list.append(cid)
            x_list.append(i)
            y_list.append(y1[i])

def createCluster(ptracks, cid, nid):
    #Calculate pairwise distance        
    maximumdistance = 450000.0
    mclus = HashTable(50000)
    start = nid

    s = len(ptracks)
    #df_ stores pairwise similarity
    df_ = pd.DataFrame(index=ptracks, columns=ptracks)
    df_.fillna(1.0, inplace=True)

    for i in range(0,s):
        tid1 = ptracks[i]
        nt1 = tracksH[tid1]
        df_.iat[i,i] = 0.0
        for j in range(i+1, s):
            tid2 = ptracks[j]
            nt = tracksH[tid2]
            dist = distance(nt, nt1)
            df_.iat[i,j] = dist
            df_.iat[j,i] = dist
            print (tid1, tid2, dist)
    for i in range(0,s):
        tid1 = ptracks[i]
        done = 0
        for j in range(i+1, s-1):
            tid2 = ptracks[j]
            dist = df_.iat[i,j]
            if (dist < 0.02):
                done = 1
                clus1 = nhclus[tid1]
                clus2 = nhclus[tid2]
                if (clus1 == None and clus2 == None):
                    nhclus[tid1] = nid
                    nhclus[tid2] = nid
                    mclus[nid] = []
                    mclus[nid].append(tid1)
                    mclus[nid].append(tid2)
                    nid = nid + 1
                elif (clus1 == None):
                    if (closeToMembers(tid1, nhclus[tid2], df_, mclus, 0.09)):
                        nhclus[tid1] = nhclus[tid2]
                        mclus[nhclus[tid1]].append(tid1)
                elif (clus2 == None):
                    if (closeToMembers(tid2, nhclus[tid1], df_, mclus, 0.09)):
                        nhclus[tid2] = nhclus[tid1]
                        mclus[nhclus[tid2]].append(tid2)
        if (done == 0 and nhclus[tid1] == None):
            nhclus[tid1] = nid
            mclus[nid] = []
            mclus[nid].append(tid1)
            nid = nid + 1

    for j in range(start, nid):
        print (j, mclus[j])

    return nid

def postProcessClusters(cdf, ptracks, centroids, filename, mode, rrange):
    agg = cdf.groupby('closest')
    nidclus = 1
    cid_list = []
    for gid, group in agg:
        ltrack = []
        trackseries = group['track_id']
        for i, track in trackseries.iteritems():
            ltrack.append(track)
        nidclus = createCluster(ltrack, gid, nidclus)
    for track in ptracks:
        cid_list.append(nhclus[track])
    
    df_tracks_andh_clusters = pd.DataFrame()
    df_tracks_andh_clusters['track_id'] = ptracks
    df_tracks_andh_clusters['cluster'] = cid_list
    df_tracks_andh_clusters.to_csv(filename, mode=mode)
    

def printClusterCentersToFile(cdf, centroids, filename, mode, rrange):
    df_clus = pd.DataFrame()
    cid_list = []
    x_list = []
    y_list = []
    agg = cdf.groupby('closest')
    for gid, group in agg:
        trackseries = group['track_id']
        printClusterCenter(centroids, gid, rrange, cid_list, x_list, y_list)
    df_clus['cid'] = cid_list
    df_clus['x'] = x_list
    df_clus['y'] = y_list

    df_clus.to_csv(filename, mode=mode, header=False)

def assignClustersToTracks(cdf, centroids):
    j = 0
    clusters = cdf['closest']
    tracks = cdf['track_id']
    for track in tracks:
        nt = tracksH[track]
        nt.cid = clusters[j]
        j = j + 1
    
def printOutliersToFile(cdf, centroids, filename, mode):
    agg = cdf.groupby('closest')
    for gid, group in agg:
        trackseries = group['track_id']
        numtracks = 0
        for track in trackseries:
            numtracks = numtracks + 1
        if (numtracks < 50):
            print ("outlier cluster: ", gid, numtracks)
        else:
            print ("valid cluster: ", gid, numtracks)

def printClustersToFile(cdf, filename, mode):
    df_tracks_with_clusters = pd.DataFrame()
    df_tracks_with_clusters['track_id'] = cdf['track_id'].astype(int)
    df_tracks_with_clusters['closest'] = cdf['closest'].astype(int)
    df_tracks_with_clusters['SPAT'] = cdf['SPAT']
    df_tracks_with_clusters['class'] = cdf['class']
    df_tracks_with_clusters['isAnomalous'] = cdf['isAnomalous']
    df_tracks_with_clusters.dropna(inplace=True)
    df_tracks_with_clusters.to_csv(filename,  float_format='%.f', mode=mode, header=False, index=False)

def printClusters(cdf):
    j = 0
    clusters = cdf['closest']
    tracks = cdf['track_id']
    for track in tracks:
        nt = tracksH[track]
        nt.cid = clusters[j]
        print (nt.tid, nt.cid)
        j = j + 1

def KMeansInit(ptracks, rrange, idclus, isVeh):
    #Calculate pairwise distance        
    maximumdistance = 450000.0
    mclus = HashTable(50000)

    hardS = 8
    if (isVeh == 0):
        hardS = 2
    s = min(len(ptracks), hardS)
    #s = len(ptracks)
    #df_ stores pairwise similarity
    df_ = pd.DataFrame(index=ptracks[0:s], columns=ptracks[0:s])
    df_.fillna(1.0, inplace=True)

    #array storing clusters
    aclus = np.zeros(shape=s)
    startIndex = idclus

    data = np.empty(shape=(s,2))
    directionalarray = np.zeros(shape=s)
    for i in range(0,s):
        tid1 = ptracks[i]
        nt1 = tracksH[tid1]
        df_.iat[i,i] = 0.0
        for j in range(i+1, s):
            tid2 = ptracks[j]
            nt = tracksH[tid2]
            dist = distance(nt, nt1)
            df_.iat[i,j] = dist
            df_.iat[j,i] = dist
    for i in range(0,s):
        tid1 = ptracks[i]
        done = 0
        for j in range(i+1, s-1):
            tid2 = ptracks[j]
            dist = df_.iat[i,j]
            #if (dist < 0.09):
            if (dist < 0.02):
                done = 1
                clus1 = hclus[tid1]
                clus2 = hclus[tid2]
                if (clus1 == None and clus2 == None):
                    hclus[tid1] = idclus
                    hclus[tid2] = idclus
                    aclus[i] = hclus[tid1]
                    aclus[j] = hclus[tid2]
                    mclus[idclus] = []
                    mclus[idclus].append(tid1)
                    mclus[idclus].append(tid2)
                    idclus = idclus + 1
                elif (clus1 == None):
                    #if (closeToMembers(tid1, hclus[tid2], df_, mclus, 0.15)):
                    if (closeToMembers(tid1, hclus[tid2], df_, mclus, 0.09)):
                        hclus[tid1] = hclus[tid2]
                        mclus[hclus[tid1]].append(tid1)
                    #else:
                        #hclus[tid1] = idclus
                        #mclus[idclus] = []
                        #mclus[idclus].append(tid1)
                        #idclus = idclus + 1
                    aclus[i] = hclus[tid1]
                elif (clus2 == None):
                    #if (closeToMembers(tid2, hclus[tid1], df_, mclus, 0.15)):
                    if (closeToMembers(tid2, hclus[tid1], df_, mclus, 0.09)):
                        hclus[tid2] = hclus[tid1]
                        mclus[hclus[tid2]].append(tid2)
                    #else:
                        #hclus[tid2] = idclus
                        #mclus[idclus] = []
                        #mclus[idclus].append(tid2)
                        #idclus = idclus + 1
                    aclus[j] = hclus[tid2]
        if (done == 0 and hclus[tid1] == None):
            hclus[tid1] = idclus
            mclus[idclus] = []
            mclus[idclus].append(tid1)
            idclus = idclus + 1
            aclus[i] = tid1

    #for j in range(startIndex,idclus):
    for j in range(1,idclus):
        print (j, mclus[j])
    if (rrange == 1280):
        df_.to_csv('similarity.csv')
    else:
        df_.to_csv('vsimilarity.csv')

    # centroids[i] = track index in ptracks
    centroids = []
    #for i in range(startIndex, idclus):
    for i in range(1, idclus):
        #if (len(mclus[i]) > 1):
        centroids.append(tracksH[mclus[i][0]])

    return idclus, centroids

def main():
    delta = 1
    xindex = 2
    yindex = 3
    cindex = 4
    tindex = 5
    sindex = 7
    rangex = 1280
    rangey = 960

    #Create x and y axes with ticks at delta
    x = np.linspace(0,rangex,delta,endpoint=False)
    y = np.linspace(0,rangey,delta,endpoint=False)
    
    #Data preprocessing: identify pedestrian tracks and interpolate tracks
    tracks = df['track_id'].unique()

    ntrackobj = np.zeros(shape=2)
    alist=[]
    phases = df['RedStatus'].unique().tolist()
    #objects = df['class'].unique().tolist()
    objects = ['ped', 'non-ped']
    gtracks = [[[] for j in range(len(phases))] for i in range(len(objects))]
    vtracks=[]

    p=0
    for tid in tracks:
        #if (tid > 600):
            #break
        ndf = df.loc[df['track_id'] == tid]
        x = ndf.iloc[:,2].values
        y = ndf.iloc[:,3].values
        c = ndf.iloc[:,4].values
        t = ndf.iloc[:,8].values
        #throw away the tracks that have insufficient information
        r = np.size(x)
        if (r < fps*2):
            continue
        if (np.absolute(x[r-1] - x[0]) < 50 and np.absolute(y[r-1]-y[0]) < 50):
            continue
        step = int(fps*0.5)
        deltax = np.absolute(x[r-1] - x[0])
        detx = 0
        p = int(len(t)-fps/2)
        pindex=phases.index(t[p])
        if (c[0] == 'pedestrian'):
            cindex = 0
        else:
            cindex = 1
        #cindex=objects.index(c[0])
        if (deltax < xvertBox):  #trajectory is vertical for all practical purposes
            rangexy = rangey
            detx = 1
            vtracks.append(tid)
        else:
            rangexy = rangex
            gtracks[cindex][pindex].append(tid)
            ntrackobj[cindex] = ntrackobj[cindex] + 1
        #Preprocessing data for pedestrians
        if (np.absolute(x[fps-1]-x[0]) < maxPedSpeed and np.absolute(y[fps-1]-y[0]) < maxPedSpeed):
            k = 0
            xnsize = int(r/step)
            if (r % step != 0):
                xnsize = xnsize + 1
            xn = np.empty(shape=xnsize)
            yn = np.empty(shape=xnsize)
            for i in range(0, r, step):
                xn[k] = x[int (i/step)*step]
                yn[k] = y[int (i/step)*step]
                k = k+1
            x = xn
            y = yn
            r = np.size(x)
        #Create interpolated versions of each track
        istart = 0
        iend = 0
        d = delta
        rstart = 0
        rend = 0
        rrange = rangex
        if (detx):
            rrange = rangey
        arr = np.empty(shape=rrange)
        if (detx):
            istart = (int) (y[0]/delta)
            iend = (int) (y[r-1]/delta)
            rstart = int(min(y[0], y[r-1]))
            rend = int(max(y[0], y[r-1]))
            for i in range(0, rstart, delta):
                arr[i] = 0
            for i in range(rend, rangexy, delta):
                arr[i] = 0
        else:
            istart = (int) (x[0]/delta)
            iend = (int) (x[r-1]/delta)
            rstart = int(min(x[0], x[r-1]))
            rend = int(max(x[0], x[r-1]))
            for i in range(0, rstart, delta):
                arr[i] = 0
            for i in range(rend, rangexy, delta):
                arr[i] = 0
        if (istart > iend):
            d = -delta
            iend = iend-1
        else:
            iend = iend+1
        for i in range(istart, iend, d):
            arr[i] = linearInterpolation(i, x, y, detx)

        nt = Track(tid, x[0], x[r-1], y[0], y[r-1], rrange, detx, pindex, cindex)
        nt.arr = arr
        tracksH[tid] = nt
 
    #print ("Starting to print trajectory file")
    #printTrajectoryFileFor2008Paper('jan152pm.tra', alist)

    idclus = 1
    startIndex = 1
    #lobj = len(objects)-1
    lobj = 2
    lph = len(phases)
    odf = pd.DataFrame()
    for i in range(0,lobj):
        anotracksthreshold = int(ntrackobj[i]*0.01)
        if (anotracksthreshold < 1):
            anotracksthreshold = 1
        for j in range(0,lph):
            #startIndex = idclus
            idclus = 1
            ptracks = gtracks[i][j]
            if (len(ptracks) == 0):
                continue
            idclus, centroids = KMeansInit(ptracks, 1280, idclus, i)
            cdf = pd.DataFrame()
            cdf['track_id'] = ptracks

            KMeansIterate(KMeansAssign(cdf, centroids, ptracks, startIndex), centroids, ptracks, 1280, startIndex)

    #postProcessClusters(cdf, ptracks, centroids, 'cluster.csv', 'w', 1280)

    #printClusterCentersToFile(cdf, centroids, 'centers.csv', 'w', 1280)
    #printOutliersToFile(cdf, centroids, 'outliers.csv', 'w')
            #cdf = cdf.dropna(inplace=True)
            #cdf = cdf.astype(int)
            cdf['SPAT'] = [phases[j]] * len(ptracks)
            cdf['class'] = [objects[i]] * len(ptracks)
            agg = cdf.groupby('closest')
            for gid, group in agg:
                trackseries = group['track_id']
                numtracks = 0
                for track in trackseries:
                    numtracks = numtracks + 1
                if (numtracks < anotracksthreshold):
                    aHash[track] = 1
            alist = []
            for k in range (0, len(ptracks)):
                anom = 0
                if (aHash[ptracks[k]] != None):
                    anom = 1
                alist.append(anom)
            cdf['isAnomalous'] = alist
            printClustersToFile(cdf, del_file, 'a')
    
    print (ntrackobj[0])
    print (ntrackobj[1])
    #idclus, vcentroids = KMeansInit(vtracks, 960, idclus)
    #cdfv = pd.DataFrame()
    #cdfv['track_id'] = vtracks

    #KMeansIterate(KMeansAssign(cdfv, vcentroids, vtracks), vcentroids, vtracks, 960)

    #printClusterCentersToFile(cdfv, vcentroids, 'centers.csv', 'a', 960)
    #printOutliersToFile(cdfv, vcentroids, 'outliers.csv', 'a')
    #printClustersToFile(cdfv, 'cluster.csv', 'a')
    
if __name__=="__main__":
        main()
