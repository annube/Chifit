#!/bin/zsh
g++ dynload.cpp `R CMD config --cppflags` `R CMD config --ldflags` -lgslcblas -lgsl -lginac -L../../src -lchifit -o dynload
