import pandas as pd
import numpy as np
import math
import sys
import datetime
import binascii
from pathlib import Path
import glob

#@dataclass
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


isBaseSet = 0
i = 0
cycle_no = 0
cycle_mark = -1
reportTime = pd.to_datetime('today')
j = -1
timeIndex = 0
eventCIndex = 1
eventPIndex = 2

def processEventList(eventList, currentTime):
    global isBaseSet
    global i
    global cycle_no
    global cycle_mark
    global reportTime
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
                newdf.loc[j]['System_time'] = reportTime
                newdf.loc[j]['Phase'] = '0' * paddedzeroes + hexstr
                newdf.loc[j]['BinPhase'] = sigstatus.decode()
                newdf.loc[j]['RedPhase'] = fhexstr[4:6]
                newdf.loc[j]['Cycle'] = cycle_no
            reportTime = reportTime + step
    for row in eventList:
        if (row[eventCIndex]  == 1):     #Phase begin green 
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

baseTime = pd.to_datetime('today')
step = datetime.timedelta(seconds=1)

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

data_folder = Path('./')
files_to_open = options['atspm']
#files_to_open = '5525_Events_20190710T0600.csv'
#files_to_open = '5525_Events_20190715T0630.csv'
#files_to_open = '5525_Events_20190604T0000.csv'
#files_to_open = 'TRAF_07356_2019_03_22_14*.dat.txt'
all_files  = sorted(data_folder.glob(files_to_open))

#all_files='logs/TRAF_07356_2019_03_22_1400.dat.txt'
newdf = pd.DataFrame(columns=['System_time','Phase','BinPhase','RedPhase','Cycle'], index = range(800000))
sigstatus = bytearray(b'000000000000000011111111')
phaseInfo = []
numPhases = 16
k = 0
while k <= numPhases:      #assign and initialize a phase cycle list for each phase
    cycleList = []
    cycleList.insert(0, PhaseCycle(0, 0, 0, 0, 0, 0, 0, 0, 1, 0)) #Wont the first location have the last phase information? 
    phaseInfo.append(cycleList)
    k = k + 1
  
eventList = []
previousTime = -1
df2 = pd.DataFrame()
timelist = []
eventlist = []
paramlist = []
for file_ in all_files:        #process controller CSVs one by one
    if __name__ == "__main__":
        #file_ = sys.argv[1]
        print ("Processing", file_)
        lineno = 0
        df1 = pd.read_csv(file_, skiprows=[0,1,2,3,4,5], header=None)  # read controller CSV
        #sigstatus = '-'* 24
        #sigstatus = bytearray(b'------------------------')
        for row_id, row in enumerate(df1.values):
            lineno = lineno + 1
            #print (lineno)
            currentTime = pd.to_datetime(row[timeIndex], format = '%m/%d/%Y %H:%M:%S.%f')
            timelist.append(currentTime)
            eventlist.append(row[eventCIndex])
            paramlist.append(row[eventPIndex])
df2['intersection_id'] = [options['cityintid'] for i in range(0,len(timelist))]
df2['time'] = timelist
df2['code'] = eventlist
df2['param'] = paramlist

df2.to_csv(options['reformat'], index=False)
