## Init
using Pkg
Pkg.activate(joinpath(dirname(@__FILE__), ".."))
import DAQP
## Setup QP data
H = [1.0 0; 0 1];
f = [2.0;2];
A = [1.0 0 ; 0 1];
bupper = [1.0;1];
blower = [-1.0,-1];
sense = [0;0];
## Solve QP
xstar,fval,exitflag,info= DAQP.quadprog(H,f,A,bupper,blower,sense);
@info "Optimal solution:" xstar
@info "Exit flag:" exitflag 
@info "Solver info:" info
