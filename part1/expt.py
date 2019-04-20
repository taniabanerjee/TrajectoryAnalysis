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
from frechetdist import frdist

##Connect to an existing database
mydb = pymysql.connect(host='trafficdatarepo.cc8u1njcqcdl.us-east-2.rds.amazonaws.com', user='dhruv', passwd='11235813', db='testdb', port=3306)

mycursor = mydb.cursor()

#query="select track_id, x, y, w, h  from TrackInfo where date='2019-01-15' and time='02-00-PM' and intersection_id=2 order by track_id asc, frame_id asc;"
query="select a.track_id, a.x+0.5*w, a.y+0.5*h, b.BinPhase from TrackInfo a, SPAT b where a.timestamp = b.timestamp order by track_id asc, frame_id asc;"

#query = "SELECT * FROM orders WHERE date_time BETWEEN ? AND ?"
#df = pd.read_sql(query, mydb,  params=(start_date, end_date))

df = pd.read_sql(query, mydb)
df.to_csv('outputxy.csv')


