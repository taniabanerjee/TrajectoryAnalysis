import pandas as pd
import numpy as np
import math
import sys
import datetime
from pathlib import Path
import glob
import os
import csv
import pymysql
import fileinput

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


##Connect to an existing database
mydb = pymysql.connect(host=options['host'], user=options['user'], passwd=options['passwd'], db=options['testdb'], port=int(options['port']))
mycursor = mydb.cursor()
#query="select time,frame_id,unique_ID,center_x,center_y,class,timestamp,intersection_id,w,h from OnlineDualTrackInfo where intersection_id=%(name)s and timestamp between %(name2)s and %(name3)s order by frame_id asc;"
query1="select * from OnlineDisplayInfo where intersection_id=%(name1)s;"
query2="select * from TrackProperties where intersection_id=%(name1)s;"

#df = pd.read_sql(query, mydb, params={'dbname' : options['tracksDB'], 'name' : options['id'], 'name2' : options['start'], 'name3' : options['end']})
df = pd.read_sql(query1, mydb, params={'name1' : options['cityintid']})
df.to_csv(options['display'], index=False)
df = pd.read_sql(query2, mydb, params={'name1' : options['cityintid']})
df.to_csv(options['trackprop'], index=False)
#df.to_csv('../DataFiles/example.csv', index=False)

