%% Init 
addpath ../interfaces/matlab/
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

% Solve and compare with quadprog solition
sense = zeros(m,1,'int32') ;
[u_daqp,fval_u_daqp, flag, time] =  daqpmex(M',d,sense);
x_daqp = -(R\(u_daqp+v));
[xref,fval_ref] = quadprog(H,f,A,b);
err=norm(x_daqp-xref)
