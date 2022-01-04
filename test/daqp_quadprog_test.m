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
[x_daqp,fval_daqp, flag_daqp, info_daqp] =  daqpmex_quadprog(H',f,A',[ub;bupper],[lb;blower],[],[],sense,daqp_opts);
err=norm(x_daqp-xref)
info_daqp
fval_daqp-fval_ref
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
%% Random qp  bounds (soft)
%Generate problem 
addpath ../interfaces/matlab/
n=50;
m = 100;
ms = 25;

Mr = randn(n);
H =  Mr'*Mr;
f = randn(n,1);
A = randn(m,n);
bupper = randn(m,1);
blower = bupper-rand(m,1);
ub = rand(ms,1);
lb = -rand(ms,1);
S = ones(m,1);

R = chol(H);
M = [eye(ms,n)/R;A/R];
v = R'\f;
dupper = [ub;bupper]+M*v;
dlower= [lb;blower]+M*v;

% Solve and compare with quadprog solution

daqp_opts = daqp_options();
daqp_opts.rho_soft = 1e-3;
daqp_opts.progress_tol =1e-6;
daqp_opts.iter_limit=2500;

Asoft_upper = [A,-S];
Asoft_lower= [A,S];
Hsoft = blkdiag(H,(1/daqp_opts.rho_soft)^2);
fsoft = [f;0];
[xref,fval_ref,flag_ref] = cplexqp(Hsoft,fsoft,[Asoft_upper;-Asoft_lower],[bupper;-blower],[],[],[lb;-inf(n-ms+1,1)],[ub;inf(n-ms+1,1)]);

%Ref
sense = zeros(2*m+ms,1,'int32') ;
sense2 = zeros(m+ms,1,'int32') ;
sense2(ms+1:end)=8*S;
[x_daqp_plain,fval_daqp, flag_daqp_plain, info_daqp_plain] =  daqpmex_quadprog(Hsoft',fsoft,[Asoft_upper;Asoft_lower]',[ub;bupper;1e9*ones(m,1)],[lb;-1e9*ones(m,1);blower],[],[],sense,daqp_opts);

[x_daqp,fval_daqp, flag_daqp, info_daqp] =  daqpmex_quadprog(H',f,A',[ub;bupper],[lb;blower],[],[],sense2,daqp_opts);
daqp_opts.eps_prox=1e4;
[x_daqp_prox,fval_daqp_prox, flag_daqp_prox, info_daqp_prox] =  daqpmex_quadprog(H',f,A',[ub;bupper],[lb;blower],[],[],sense2,daqp_opts);
min_slack = min([ub;-lb;bupper;-blower]+[-eye(ms,n);eye(ms,n);-A;A]*x_daqp);
x_daqp_ext = [x_daqp;info_daqp.soft_slack*daqp_opts.rho_soft];
x_daqp_prox_ext = [x_daqp_prox;info_daqp_prox.soft_slack*daqp_opts.rho_soft];
fval_daqp_ext = 0.5*x_daqp_ext'*Hsoft*x_daqp_ext+fsoft'*x_daqp_ext;

fval_daqp_prox_ext = 0.5*x_daqp_prox_ext'*Hsoft*x_daqp_prox_ext+fsoft'*x_daqp_prox_ext;

err=norm(x_daqp_ext-xref);
fprintf('\n\n============= Result: ============ \n flag: %d|fval diff: %f|norm_diff : %f |min slack: %f\n',flag_daqp,fval_ref-fval_daqp_ext,err,min_slack);
info_daqp
info_daqp_prox
norm(x_daqp-x_daqp_prox)
fval_daqp_ext-fval_daqp_prox_ext


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

