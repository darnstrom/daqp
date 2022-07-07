%% Setup QP data
n = 200;
m = 500;
ms = 20;
M = randn(n,n);
H = M'*M;
%H(1:ms,1:ms) = 0.01*eye(ms);
%H(ms+1:end,1:ms) = 0; 
%H(1:ms,ms+1:end) = 0;
%H = eye(4);
f = 100*randn(n,1); 
A = randn(m,n);
bupper = 20*rand(m,1);
blower = -20*rand(m,1);
bupper_tot = [ones(ms,1);bupper];
blower_tot = [zeros(ms,1);blower];
sense = int32(zeros(m+ms,1));
sense(1:ms) = sense(1:ms)+16;
%% Solve QP (quadprog interface)
[xstar,fval,exitflag,info] = daqp.quadprog(H,f,A,bupper_tot,blower_tot,sense);
%xstar
%bupper-A*xstar
%-blower+A*xstar
fval = 0.5*xstar'*H*xstar+f'*xstar

%% Gurobi
model.Q = 0.5*sparse(H);
model.A = sparse([A;-A]);
model.rhs = [bupper;-blower];
model.obj = f;
model.sense='<';
model.vtype=repelem('BC',[ms,n-ms]);
model.lb = -inf(n,1);

params.Threads=0;
results = gurobi(model,params)

xgrb = results.x;
fgrb = 0.5*xgrb'*H*xgrb+f'*xgrb

%% Compare 
norm(xstar-xgrb)
fgrb-fval
