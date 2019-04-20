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

mydb = pymysql.connect(host='trafficdatarepo.cc8u1njcqcdl.us-east-2.rds.amazonaws.com', user='dhruv', passwd='11235813', db='testdb', port=3306)

mycursor = mydb.cursor()
query="select frame_id, track_id, center_x, center_y, class, timestamp, intersection_id, SPAT, cycle, cluster from DisplayInfo where cycle=371 order by track_id asc, frame_id asc;"
df = pd.read_sql(query, mydb)
df.to_csv('cycle371.csv', index=False)

