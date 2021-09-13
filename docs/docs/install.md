---
layout: page
title: Installation
permalink: /install/
nav: 2 
---

## Installing from source 
In the cloned directory (`<path were cloned>/daqp`), run the following commands

```
mkdir build
cd build
cmake ..
cmake --build .
```

## Installing MATLAB interface
To build mex-files, run the following commands from the `daqp` directory
```
cd interfaces/matlab
source ./makemex.sh 
```
See `interfaces/test.m` for a MATLAB example 
