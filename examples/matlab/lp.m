%% Setup LP data
n = 100;
m = 500;
ms = 50;
f = randn(n,1);
A = randn(m,n);
bupper = 10*rand(ms+m,1);
blower = -10*rand(ms+m,1);
sense = zeros(ms+m,1,'int32'); 
cplex_opts = cplexoptimset('cplex');
%% Solve LP (linprog interface)
cplex_opts.lpmethod=1;
[x,fval,exitflag,info] = daqp.linprog(f,A,bupper,blower,sense);
[xref,~,~,output_cplex] = cplexlp(f,[A;-A],[bupper(ms+1:end);-blower(ms+1:end)],[],[],[blower(1:ms);-inf(n-ms,1)],[bupper(1:ms);inf(n-ms,1)],[],cplex_opts);
norm(x-xref)
info
output_cplex
