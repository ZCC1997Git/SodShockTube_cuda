#!/bin/bash
cmake --build ../build --target SodWave_CPU
cmake --build ../build --target SodWave_GPU
cp ../build/SodWave_CPU ./SodWave_CPU
cp ../build/SodWave_GPU ./SodWave_GPU
./SodWave_CPU
./SodWave_GPU
python3 Plot.py

