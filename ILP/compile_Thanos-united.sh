#!/bin/bash
./dos2unix Thanos-United.cpp
g++ -std=c++11 -O3 -Wall Thanos-United.cpp -o Thanos-United -I/opt/gurobi811/linux64/include -L/opt/gurobi811/linux64/lib -lgurobi_c++ -lgurobi81
