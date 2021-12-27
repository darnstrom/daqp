%% Init 
addpath ../interfaces/matlab/
%% Random lp 
%Generate problem 
n=10;
m = 100;

f = randn(n,1);
A = randn(m,n);
b = rand(m,1);

% Solve and compare with linprog solution
sense = zeros(length(b),1,'int32');
[x_dalp,~,flag,CPUtime,iter] = daqpproxmex(A',b,zeros(0,0),f,1,sense);
[xref,fval_ref] = linprog(f,A,b);
err=norm(x_dalp-xref)
