#!/bin/bash

cat << EOS
########################################
I DO buck up fortran source programs.

Process number : $$

########################################
EOS

echo "If you proceed, type 1."
read judge
if [ $judge -eq 1 ]; then
  date=`date "+%Y-%m-%d"`  
  if [ -e "../ver$date" ]; then
    echo "ver$date already exists, so I copied the programs without making a directory."
    cp -r ./ ../ver$date

    echo ""
    echo "Please write the update info."
    read info
    echo "$info" > ../ver$date/Update_info

  else
    echo "ver$date does not exist, so I made a directory to store them, and copied them to the dir."
    cp -r ./ ../ver$date

    echo ""
    echo "Please write the update info."
    read info
    echo "$info" > ../ver$date/Update_info
  fi
else
  echo "Do not BACK UP!"
fi 
