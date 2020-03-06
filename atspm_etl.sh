#!/bin/bash

cd /home/ubuntu/data;
curr_date=$(date -d "$date -1 days" +"%Y_%m_%d")
#curr_date="2019_11_27"
#curr_date=$1
mkdir -p $curr_date
cd $curr_date
scp teller@128.227.170.246:ATSPM/*${curr_date}* .
#scp teller@128.227.170.246:ATSPM/TRAF* .
shopt -s nullglob
files=(*)
if [ -z "$files" ]; then
    echo "No files found for date ${curr_date}!"
    return 1
fi
echo "${#files[@]} files downloaded!"
for i in "${files[@]}"; do
  intersection_id=${i:5:5}
  sed -i '1,6d' ${i}
  sed -i "s/\r//g" ${i}
  sed -i "s/$/, $intersection_id, Gainesville, FL/" -i ${i}
  PGPASSWORD=12345678 psql -h atspm.cc8u1njcqcdl.us-east-2.rds.amazonaws.com -p 5432 -U kunwardeep atspm -c "\copy atspm(\"received_date\", \"event_code\", \"event_parameter\", \"intersection_id\", \"city\", \"state\") FROM '${i}' DELIMITER ',' CSV;"
done
