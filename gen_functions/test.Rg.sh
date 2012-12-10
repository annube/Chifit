#!/bin/zsh
g++ -c test.tmFS.R.g.cpp
g++ -O3 -c tmFS.R.g.cpp
g++ test.tmFS.R.g.o tmFS.R.g.o -o test.tmFS.R.g
./test.tmFS.R.g
