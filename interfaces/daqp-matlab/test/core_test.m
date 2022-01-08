classdef core_test < matlab.unittest.TestCase
% Sovle 1000 randomly generated QPs
  methods (Test)
	function random_feasible_QPs(testCase)
	  % Test on randomly generated feasible QPs
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
	function infeasible_QPs(testCase) 
	  H = eye(2);
	  f = zeros(2,1);
	  A = [1 1];
	  bupper = [1;1;20];
	  blower = [1;1;19];
	  sense = zeros(3,1);
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
  end
end

