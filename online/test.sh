#/bin/bash 

rm ../7356/TrackProperties.csv
python3 match.py ../7356/intersection3.txt
python3 display.py ../7356/intersection3.txt
diff ../7356/TrackProperties.csv ../7356/goldenTrackProperties.csv
