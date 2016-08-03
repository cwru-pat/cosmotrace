#!/bin/bash

echo "Running some tests."

g++ -std=c++11 -O3 -ffast-math raytrace.cc -lz -o raytrace
if [ $? -ne 0 ]; then
    echo "Error: unable to compile."
    exit 1
fi

for i in `seq 1 6`;
do
  ./raytrace $i
  if [ $? -ne 0 ]; then
      echo "Error running raytrace code."
      exit 1
  fi
done
