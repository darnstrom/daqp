%% Init 
addpath ../interfaces/matlab/
%% Random qp  bounds
%Generate problem 
addpath ../interfaces/matlab/
n=100;
m = 500;
ms = 50;

Mr = randn(n);
H =  Mr'*Mr;
f = randn(n,1);
A = randn(m,n);
bupper = rand(m,1);
blower = -rand(m,1);
ub = rand(ms,1);
lb = -rand(ms,1);

R = chol(H);
M = [eye(ms,n)/R;A/R];
v = R'\f;
dupper = [ub;bupper]+M*v;
dlower= [lb;blower]+M*v;

% Solve and compare with quadprog solution
[xref,fval_ref] = quadprog(H,f,[A;-A],[bupper;-blower],[],[],lb,ub);

daqp_opts = daqp_options();
daqp_opts.progress_tol =1e-6;

sense = zeros(m+ms,1,'int32') ;
[x_daqp,fval_x_daqp, flag_daqp, info_daqp] =  daqpmex_quadprog(H',f,A',[ub;bupper],[lb;blower],[],[],sense,daqp_opts);
err=norm(x_daqp-xref)
info_daqp
fval_x_daqp-fval_ref
min([ub;-lb;bupper;-blower]+[-eye(ms,n);eye(ms,n);-A;A]*x_daqp)

%% Temp dissmiss
daqp_opts.eps_prox = 1e-6;
%[x_prox,fval_x_prox, flag_prox, info_prox] =  daqpmex_quadprog(H',f,A',[ub;bupper],[lb;blower],[],[],sense,daqp_opts);
%
%
%fval_daqp = 0.5*x_daqp'*H*x_daqp+f'*x_daqp;
%fval_prox = 0.5*x_prox'*H*x_prox+f'*x_prox;
%
%fval_daqp-fval_prox
err
%% Random qp 
%Generate problem 
n=10;
m = 100;

Mr = randn(n);
H =  Mr'*Mr;
f = randn(n,1);
A = randn(m,n);
b = rand(m,1);

R = chol(H);
M = A/R;
v = R'\f;
d = b+M*v;

% Solve and compare with quadprog solution

daqp_opts = daqp_options();

sense = zeros(m,1,'int32') ;
[x_daqp,fval_x_daqp, flag, info_daqp] =  daqpmex_quadprog(H',f,A',b,sense,daqp_opts);
[xref,fval_ref] = quadprog(H,f,A,b);
err=norm(x_daqp-xref)

%daqp_opts.eps_prox = 1e-6;
%[x_prox,fval_x_prox, flag, info_prox] =  daqpmex_quadprog(H',f,A',b,sense,daqp_opts);


%fval_daqp = 0.5*x_daqp'*H*x_daqp+f'*x_daqp;
%fval_prox = 0.5*x_prox'*H*x_prox+f'*x_prox;
%
%fval_daqp-fval_prox
info_daqp
%% Random lp 
%Generate problem 
n=10;
m = 100;

f = randn(n,1);
A = randn(m,n);
b = rand(m,1);


daqp_opts = daqp_options();
daqp_opts.eps_prox=1;
sense = zeros(m,1,'int32') ;
[x_daqp,fval_x_daqp, flag, daqp_info] =  daqpmex_quadprog([],f,A',b,sense,daqp_opts);
[xref,fval_ref] = linprog(f,A,b);
err=norm(x_daqp-xref)
daqp_info
