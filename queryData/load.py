import sys
import csv
import pymysql

pymysql.install_as_MySQLdb()

#Create a brand new database..
#mydb = pymysql.connect(host='localhost', user='root', passwd='Sbtb11Tb.')
#mycursor = mydb.cursor()
#mycursor.execute("CREATE DATABASE trackDatabase")

#Connect to an existing database
mydb = pymysql.connect(host='localhost', user='root', passwd='Sbtb11Tb.', db='trackDatabase', autocommit='true')
mycursor = mydb.cursor()
mycursor.execute("CREATE TABLE IF NOT EXISTS TrackInfo(time FLOAT, frame_id INT, track_id INT, class varchar(20), x INT, y INT, w INT, h INT, center_x FLOAT, center_y FLOAT, intersection_id INT, timestamp DATETIME)")
mycursor.execute("CREATE TABLE IF NOT EXISTS SPAT(timestamp DATETIME, Phase varchar(6), BinPhase varchar(24), Cycle INT)")

