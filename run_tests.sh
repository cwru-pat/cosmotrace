#!/bin/bash

echo "Running some tests."

g++ -std=c++11 -O3 -ffast-math raytrace.cc -lz -o raytrace
if [ $? -ne 0 ]; then
    echo "Error: unable to compile."
    exit 1
fi

./raytrace
