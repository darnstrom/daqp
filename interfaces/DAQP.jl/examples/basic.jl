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
sense = Vector{Cint}([0;0]);
## Solve QP
xstar,fval,exitflag,info= DAQP.quadprog(H,f,A,bupper,blower,sense);
@info "Optimal solution:" xstar
@info "Exit flag:" exitflag 
@info "Solver info:" info
## Solve persistent 
d = DAQP.Model();
DAQP.setup(d,H,f,A,bupper,blower,sense); 
xstar,fval,exitflag,info = DAQP.solve(d);

## Update model
f = -1*f;
DAQP.update(d,H,f,A,bupper,blower,sense);
xstar,fval,exitflag,info = DAQP.solve(d);
