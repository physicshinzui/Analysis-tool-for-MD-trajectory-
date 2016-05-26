#!/bin/bash 
set -e

while [ -e ./output.dat ]; do
  echo -e "Progress (running now): `wc output.dat | awk '{print $1}'` \r\c"
  sleep 1.0
done

