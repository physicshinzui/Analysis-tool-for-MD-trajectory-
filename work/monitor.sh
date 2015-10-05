#!/bin/bash 
set -e

while [ -e ./output.dat ]; do
  wc output.dat | awk '{print $1}'
  sleep 0.5
done
