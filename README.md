 
# DAQP
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![](https://img.shields.io/badge/docs-online-brightgreen)](https://darnstrom.github.io/daqp)

DAQP is a dual active-set solver that solves convex quadratic programs of the form
```
minimize        0.5 x' H x + f' x

subject to      l  <=  x  <= u
		        bl <=  Ax <= bu.
```
Binary constraints of the form $A x \in \lbrace b_l, b_u \rbrace$ are also supported, allowing for mixed-integer quadratic programs to be solved.

The code is written in C and is *library free*. DAQP can be interfaced to C, Julia, MATLAB, and Python. 

## Quick installation

### Julia  <img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/julia/julia-original.svg" height="20"/>
```bash
] add DAQP
```

### MATLAB  <img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/matlab/matlab-original.svg" height="20"/> 
```bash
websave('install_daqp','https://raw.githubusercontent.com/darnstrom/daqp/master/interfaces/daqp-matlab/install_daqp.m')
install_daqp
```

### Python  <img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/python/python-original.svg" height="20" />
```bash
pip install daqp
```

### Details 
See the [Documentation](https://darnstrom.github.io/daqp/start) for more detailed installation instructions and basic use of the interfaces. 

## Citing DAQP
```
@article{arnstrom2022dual,
  author={Arnström, Daniel and Bemporad, Alberto and Axehill, Daniel},
  journal={IEEE Transactions on Automatic Control},
  title={A Dual Active-Set Solver for Embedded Quadratic Programming Using Recursive {LDL}$^{T}$ Updates},
  year={2022},
  volume={67},
  number={8},
  pages={4362-4369},
  doi={10.1109/TAC.2022.3176430}
}
```
