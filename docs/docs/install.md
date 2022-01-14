---
layout: page
title: Installation
permalink: /install/
nav: 2 
parent: Getting started 
---
<details open markdown="block">
<summary>
Table of contents
</summary>
{: .text-delta }
1. TOC
{:toc}
</details>


## Installing from source 
The installation requires [CMake](https://cmake.org/) and [GCC](https://gcc.gnu.org/).

In the cloned directory (i.e, `<path were cloned>/daqp`), run the following commands

```shell
mkdir build
cd build
cmake ..
cmake --build .
```

To copy the header files and libraries into `includedir` and `libdir`, respectively, run the following command in the `build` folder 
```shell
cmake --build . --target install
```


## Installing the MATLAB interface
If you also want to generate the mex-file used in the MATLAB interface when building from source, pass the flag `MATLAB` when setting up cmake. That is, replace the above steps with
```shell
mkdir build
cd build
cmake .. -DMATLAB=True
cmake --build .
```

## Installing the Julia interface
In the REPL run the command 
```julia
] add (path_to_DAQP)
```
where `path_to_DAQP` is the location of the `DAQP.jl` subdirectory

## Installing the Python interface
Move to the `daqp-python` subdirectory and call pip from the shell:
```shell
pip install .
```
