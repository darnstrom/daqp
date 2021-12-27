%% Init 
addpath ../interfaces/matlab/
%% Random qp 
%Generate problem 
n=10;
m = 100;

eps_prox = 1e-3;
Mr = randn(n);
Hnom= Mr'*Mr;
H = Hnom+eps_prox*eye(n);
f = randn(n,1);
A = randn(m,n);
b = rand(m,1);

R = chol(H);
M = A/R;
v = R'\f;
d = b+M*v;
Rp = R+diag(-diag(R)+1./diag(R)); % Invert R diags..
Rp= Rp'; 
Rp=Rp(tril(true(size(Rp)))); % matlab colmajor, daqp rowmajor...

% Solve and compare with quadprog solution
sense = zeros(m,1,'int32') ;
[x_daqpprox,~,exitflag_daqp,cpuTime,iters] = daqpproxmex(M',b,Rp,f,eps_prox,sense);
fval = 0.5*x_daqpprox'*Hnom*x_daqpprox+f'*x_daqpprox;
[xref,fval_ref] = quadprog(Hnom,f,A,b);
err=norm(x_daqpprox-xref)
fval-fval_ref
