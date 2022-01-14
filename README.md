## DAQP 
DAQP is a dual active-set solver that solves convex quadratic programs in the form 
```
minimize        0.5 x' H x + f' x

subject to      l  <=  x  <= u
		bl <=  Ax <= bu.
```

The code is written in C and is *library free*. DAQP can be interfaced to C, Julia, MATLAB, and Python. 

See [Documentation](https://darnstrom.github.io/daqp/) for an installation guide and basic use of the interfaces. 

## Citing DAQP 
```
@article{arnstrom2021dual,
  title={A Dual Active-Set Solver for Embedded Quadratic Programming Using Recursive LDL' Updates},
  author={Arnstr{\"o}m, Daniel and Bemporad, Alberto and Axehill, Daniel},
  journal={arXiv preprint arXiv:2103.16236},
  year={2021}
}
```

