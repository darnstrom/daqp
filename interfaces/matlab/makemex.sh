#!/bin/bash
mex -O ../../build/libdaqpstat.a -I../../include daqpmex.c
mex -O ../../build/libdaqpstat.a -I../../include daqpmexeq.c
mex -O ../../build/libdaqpstat.a -I../../include daqpproxmex.c
