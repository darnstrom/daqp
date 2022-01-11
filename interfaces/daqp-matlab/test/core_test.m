classdef core_test < matlab.unittest.TestCase
% Sovle 1000 randomly generated QPs
  methods (Test)
	function random_feasible_QPs(testCase)
	  % Test on randomly generated feasible QPs
	  rng('default');
	  nQPs = 100;
	  n = 100; m = 500; ms = 50;
	  nAct = 80;
	  kappa = 1e2;
	  tol = 1e-4;
	  solve_times = zeros(nQPs,1);
	  for i = 1:nQPs
		[xref,H,f,A,bupper,blower,sense]=generate_test_QP(n,m,ms,nAct,kappa);
		d = daqp();
		d.setup(H,f,A,bupper,blower,sense);
		[x,fval,exitflag, info] = d.solve(); 

		testCase.verifyEqual(exitflag,int32(1));
		testCase.verifyLessThan(norm(x-xref),tol);
		solve_times(i) = info.solve_time;
	  end
	  fprintf('========================== DAQP =============================\n')
	  fprintf('Solve times [s]: |avg: %2.6f| max: %2.6f| min %2.6f|\n',mean(solve_times),max(solve_times),min(solve_times))
	  fprintf('=============================================================\n')
	end
	function prox_random_feasible_QPs(testCase)
	  % Test on randomly generated feasible QPs
	  rng('default');
	  nQPs = 100;
	  n = 100; m = 500; ms = 50;
	  nAct = 80;
	  kappa = 1e2;
	  tol = 1e-4;
	  solve_times = zeros(nQPs,1);
	  for i = 1:nQPs
		[xref,H,f,A,bupper,blower,sense]=generate_test_QP(n,m,ms,nAct,kappa);
		d = daqp();
		d.settings('eps_prox',1e-2);
		d.setup(H,f,A,bupper,blower,sense);
		[x,fval,exitflag, info] = d.solve(); 

		testCase.verifyEqual(exitflag,int32(1));
		testCase.verifyLessThan(norm(x-xref),tol);
		solve_times(i) = info.solve_time;
	  end
	  fprintf('\n======================== DAQP PROX ==========================\n')
	  fprintf('Solve times [s]: |avg: %2.6f| max: %2.6f| min %2.6f|\n',mean(solve_times),max(solve_times),min(solve_times))
	  fprintf('=============================================================\n')
	end
	function infeasible_QP(testCase) 
	  H = eye(2);
	  f = zeros(2,1);
	  A = [1 1];
	  bupper = [1;1;20];
	  blower = [1;1;19];
	  sense = zeros(3,1,'int32');
	  d = daqp();
	  d.setup(H,f,A,bupper,blower,sense);
	  [~,~,exitflag,~] = d.solve();
	  testCase.verifyEqual(exitflag,int32(-1));
	  % Retry but soften the constraint
	  sense(3)=8; % soften
	  d = daqp();
	  d.setup(H,f,A,bupper,blower,sense);
	  [~,~,exitflag,info] = d.solve();
	  testCase.verifyEqual(exitflag,int32(2));
	end
	
	function equality_QP(testCase) 
	  H = eye(2); 
	  f = 10*ones(1,1);
	  A = [4 1];
	  bupper = [1;1;0];
	  blower = -[1;1;0];
	  sense = zeros(3,1,'int32');
	  sense(3)=5; % Set third constraint to equality 
	  d = daqp();
	  d.setup(H,f,A,bupper,blower,sense);
	  [x,~,exitflag,~] = d.solve();
	  testCase.verifyEqual(exitflag,int32(1));
	  testCase.verifyLessThan(norm(x-[-0.25;1]),1e-8);
	  % Distort the Hessian and ensure that the same result holds 
	  R=rand(2);
	  d = daqp();
	  d.setup(R'*R,R'*f,[eye(2);A]*R,bupper,blower,sense);
	  [x,~,exitflag,~] = d.solve();
	  testCase.verifyEqual(exitflag,int32(1));
	  testCase.verifyLessThan(norm(R*x-[-0.25;1]),1e-8);

	  % Make sure the solver detects an overdetermined set of constraints 
	  sense(1:2) = 5;
	  [x,~,exitflag,~] = d.quadprog(R'*R,R'*f,[eye(2);A]*R,bupper,blower,sense);
	  testCase.verifyEqual(exitflag,-6);

	end
	
	function random_feasible_LPs(testCase)
	  % Test on randomly generated feasible LPs
	  rng('default');
	  nQPs = 100;
	  n = 100; m = 500; ms = 50;
	  tol = 1e-4;
	  solve_times = zeros(nQPs,1);
	  for i = 1:nQPs
		[xref,f,A,bupper,blower,sense]=generate_test_LP(n,m,ms);
		d = daqp();
		d.setup([],f,A,bupper,blower,sense);
		d.settings('eps_prox',1);
		[x,fval,exitflag, info] = d.solve(); 

		testCase.verifyEqual(exitflag,int32(1));
		testCase.verifyLessThan(norm(x-xref),tol);
		if(norm(x-xref)>tol)
		  x_linprog = linprog(f,[A;-A],[bupper(ms+1:end);-blower(ms+1:end)],[],[],[blower(1:ms);-inf(n-ms,1)],[bupper(1:ms);inf(n-ms,1)]);
		  fprintf('Linprog error: %f\n',norm(xref-x_linprog));
		  disp(info)
		end
		solve_times(i) = info.solve_time;
	  end
	  fprintf('\n========================== DALP =============================\n')
	  fprintf('Solve times [s]: |avg: %2.6f| max: %2.6f| min %2.6f|\n',mean(solve_times),max(solve_times),min(solve_times))
	  fprintf('=============================================================\n')
	end
  end
end
