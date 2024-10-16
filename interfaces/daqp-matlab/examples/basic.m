%% Setup QP data
H = [1.0 0; 0 1];
f = [2.0;2];
A = [1.0 0 ; 0 1];
bupper = [1.0;1];
blower = [-1.0,-1];
sense = int32([0;0]);
%% Solve QP (quadprog interface)
[xstar,fval,exitflag,info] = daqp.quadprog(H,f,A,bupper,blower,sense);
%% Solve QP class interface 
addpath utils
n = 5;
ms = 0;
m = 5;
[xref,f,A,bupper,blower,sense]=generate_test_LP(n,m,ms);
d = daqp();
d.setup([],f,A,bupper,blower,sense);
d.settings('eps_prox',1);
[xstar,fval,exitflag,info] = d.solve();
norm(xstar-xref)
norm(f+[eye(ms,n);A]'*info.lambda)
Aext = [eye(ms,n);A]
AS = [4,5,1,3];
AAS = Aext(AS,:);
Af = AAS*f
lam = (AAS*AAS')\(-Af)
p = -f-AAS'*lam
