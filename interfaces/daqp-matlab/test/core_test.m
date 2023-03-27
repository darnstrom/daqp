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
            tol = 1e-5;
            solve_times = zeros(nQPs,1);
            primal_errors = zeros(nQPs,1);
            dual_errors = zeros(nQPs,1);
            for i = 1:nQPs
                [xref,H,f,A,bupper,blower,sense]=generate_test_QP(n,m,ms,nAct,kappa);
                d = daqp();
                d.setup(H,f,A,bupper,blower,sense);
                [x,fval,exitflag, info] = d.solve();

                testCase.verifyEqual(exitflag,int32(1));
                primal_errors(i) = norm(x-xref);
                dual_errors(i) = norm(H*x+f+[eye(ms,n);A]'*info.lambda);
                testCase.verifyLessThan(primal_errors(i),tol);
                testCase.verifyLessThan(dual_errors(i),tol);
                testCase.verifyLessThan(0.5*x'*H*x+f'*x-fval,tol);
                solve_times(i) = info.solve_time;
                solve_error(i) = norm(x-xref);
            end
            fprintf('========================== DAQP =============================\n')
            fprintf('Solve times [s]: |avg: %2.6f| max: %2.6f| min %2.6f|\n',mean(solve_times),max(solve_times),min(solve_times))
            fprintf('Primal Errors  : |avg: %2.2e| max: %2.2e| min %2.2e|\n',mean(primal_errors),max(primal_errors),min(primal_errors))
            fprintf('Dual Errors    : |avg: %2.2e| max: %2.2e| min %2.2e|\n',mean(dual_errors),max(dual_errors),min(dual_errors))
            fprintf('=============================================================\n')
        end
        function prox_random_feasible_QPs(testCase)
            % Test on randomly generated feasible QPs
            rng('default');
            nQPs = 100;
            n = 100; m = 500; ms = 50;
            nAct = 80;
            kappa = 1e2;
            tol = 1e-5;
            solve_times = zeros(nQPs,1);
            solve_errors = zeros(nQPs,1);
            for i = 1:nQPs
                [xref,H,f,A,bupper,blower,sense]=generate_test_QP(n,m,ms,nAct,kappa);
                d = daqp();
                d.settings('eps_prox',1e-2);
                d.setup(H,f,A,bupper,blower,sense);
                [x,fval,exitflag, info] = d.solve();

                testCase.verifyEqual(exitflag,int32(1));
                testCase.verifyLessThan(norm(x-xref),tol);
                testCase.verifyLessThan(norm(H*x+f+[eye(ms,n);A]'*info.lambda),tol);
                testCase.verifyLessThan(0.5*x'*H*x+f'*x-fval,tol);
                solve_times(i) = info.solve_time;
                solve_errors(i) = norm(x-xref);
            end
            fprintf('\n======================== DAQP PROX ==========================\n')
            fprintf('Solve times [s]: |avg: %2.6f| max: %2.6f| min %2.6f|\n',mean(solve_times),max(solve_times),min(solve_times))
            fprintf('Solution Errors: |avg: %2.2e| max: %2.2e| min %2.2e|\n',mean(solve_errors),max(solve_errors),min(solve_errors))
            fprintf('=============================================================\n')
        end

        function random_feasible_LPs(testCase)
            % Test on randomly generated feasible LPs
            rng('default');
            nQPs = 100;
            n = 100; m = 500; ms = 50;
            tol = 1e-5;
            solve_times = zeros(nQPs,1);
            solve_errors = zeros(nQPs,1);
            for i = 1:nQPs
                [xref,f,A,bupper,blower,sense]=generate_test_LP(n,m,ms);
                d = daqp();
                d.setup([],f,A,bupper,blower,sense);
                d.settings('eps_prox',1);
                [x,fval,exitflag, info] = d.solve(); 

                testCase.verifyEqual(exitflag,int32(1));
                testCase.verifyLessThan(norm(x-xref),tol);
                testCase.verifyLessThan(norm(f'*xref-fval),tol);
                testCase.verifyLessThan(norm(f+[eye(ms,n);A]'*info.lambda),tol);
                if(norm(x-xref)>tol)
                    x_linprog = linprog(f,[A;-A],[bupper(ms+1:end);-blower(ms+1:end)],[],[],[blower(1:ms);-inf(n-ms,1)],[bupper(1:ms);inf(n-ms,1)]);
                    fprintf('Linprog error: %f\n',norm(xref-x_linprog));
                    disp(info)
                end
                solve_times(i) = info.solve_time;
                solve_errors(i) = norm(x-xref);
            end
            fprintf('\n========================== DALP =============================\n')
            fprintf('Solve times [s]: |avg: %2.6f| max: %2.6f| min %2.6f|\n',mean(solve_times),max(solve_times),min(solve_times))
            fprintf('Solution Errors: |avg: %2.2e| max: %2.2e| min %2.2e|\n',mean(solve_errors),max(solve_errors),min(solve_errors))
            fprintf('=============================================================\n')
        end

        function feasibility(testCase)
            n = 25; m = 500;
            nQPs = 50;
            r = 1e-7;
            %r = -0.5e-5;
            for i = 1:nQPs
                [A,bupper,blower,sense] =generate_test_feasibility(n,m,r);
                [x,fval,exitflag,daqp_feas_info] = daqp.quadprog([],[],A,bupper,blower,sense);
                testCase.verifyEqual(exitflag,int32(1));
            end

            r =  -1e-6-1e-7
            for i = 1:nQPs
                [A,bupper,blower,sense] =generate_test_feasibility(n,m,r);
                [x,fval,exitflag,daqp_infeas_info] = daqp.quadprog([],[],A,bupper,blower,sense);
                testCase.verifyEqual(exitflag,int32(-1));
            end
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
            [~,~,exitflag,infeasible_info] = d.solve();
            testCase.verifyEqual(exitflag,int32(-1));
            infeasible_info
            % Retry but soften the constraint
            sense(3)=8; % soften
            d = daqp();
            d.setup(H,f,A,bupper,blower,sense);
            [~,~,exitflag,soft_info] = d.solve();
            testCase.verifyEqual(exitflag,int32(2));
            soft_info
        end

        function trivial_infeasible_QP(testCase)
            H = eye(2);
            f = zeros(2,1);
            A = [1 1];
            bupper = [1;1;1];
            blower = [2;-1;-1];
            sense = zeros(3,1,'int32');
            d = daqp();
            setup_flag = d.setup(H,f,A,bupper,blower,sense);
            testCase.verifyEqual(setup_flag,-1);
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
            [x,~,exitflag,equality_info] = d.solve();
            testCase.verifyEqual(exitflag,int32(1));
            testCase.verifyLessThan(norm(x-[-0.25;1]),1e-8);
            equality_info
            % Distort the Hessian and ensure that the same result holds
            R=rand(2);
            d = daqp();
            d.setup(R'*R,R'*f,[eye(2);A]*R,bupper,blower,sense);
            [x,~,exitflag,equality_info2] = d.solve();
            testCase.verifyEqual(exitflag,int32(1));
            testCase.verifyLessThan(norm(R*x-[-0.25;1]),1e-8);
            equality_info2

            % Make sure the solver detects an overdetermined set of constraints
            sense(1:2) = 5;
            [x,~,exitflag,info_overdetermined] = d.quadprog(R'*R,R'*f,[eye(2);A]*R,bupper,blower,sense);
            testCase.verifyEqual(exitflag,-6);

        end


        function unbounded_LP(testCase)
            f = [1; 1];
            A = [1  0];
            bupper = [1];
            blower = [-1];
            sense= int32([0]);
            d = daqp();
            d.setup([],f,A,bupper,blower,sense);
            d.settings('eps_prox',1);
            [x,fval,exitflag, unb_lp_info] = d.solve();
            testCase.verifyEqual(exitflag,int32(-3));
            unb_lp_info
        end

        function random_bnb(testCase)
            % generate and solve with daqp
            rng('default');
            n = 150; m = 300; ms = 20; me = 0;
            tol = 1e-5;
            M = randn(n,n);
            H = M'*M;
            f = 100*randn(n,1);
            f(1:ms) = -sign(f(1:ms)).*f(1:ms);
            A = randn(m,n);
            bupper = 20*rand(m,1); blower = -20*rand(m,1);
            bupper(ms+1:ms+me)=0; blower(ms+1:ms+me)=0;
            bupper_tot = [ones(ms,1);bupper];
            blower_tot = [zeros(ms,1);blower];
            sense = int32(zeros(m+ms,1));
            sense(1:ms) = sense(1:ms)+16;
            sense(ms+1:ms+me) = sense(ms+1:ms+me)+5;
            [x,fval,exitflag,daqp_bnb_info] = daqp.quadprog(H,f,A,bupper_tot,blower_tot,sense);
            testCase.verifyEqual(exitflag,int32(1));
            display(daqp_bnb_info)

            % compare with Gurobi (comparison skipped if Gurobi is not available)
            model.Q = 0.5*sparse(H);
            model.A = sparse([A;-A(me+1:end,:)]);
            model.rhs = [bupper;-blower(me+1:end)];
            model.obj = f;
            model.sense=repelem('=<',[me,2*(m-me)]);
            model.vtype=repelem('BC',[ms,n-ms]);
            model.lb = -inf(n,1);
            params.OutputFlag=0;
            if(exist('gurobi'))
                gurobi_result = gurobi(model,params)
                xref = gurobi_result.x;
                testCase.verifyLessThan(norm(x-xref),tol);
            end
        end

        function random_bnb_subopt(testCase)
            % generate and solve with daqp
            rng('default');
            n = 150; m = 300; ms = 20; me = 0;
            tol = 1e-5;
            M = randn(n,n);
            H = M'*M;
            f = 100*randn(n,1);
            f(1:ms) = -sign(f(1:ms)).*f(1:ms);
            A = randn(m,n);
            bupper = 20*rand(m,1); blower = -20*rand(m,1);
            bupper(ms+1:ms+me)=0; blower(ms+1:ms+me)=0;
            bupper_tot = [ones(ms,1);bupper];
            blower_tot = [zeros(ms,1);blower];
            sense = int32(zeros(m+ms,1));
            sense(1:ms) = sense(1:ms)+16;
            sense(ms+1:ms+me) = sense(ms+1:ms+me)+5;
            d = daqp();
            d.setup(H,f,A,bupper_tot,blower_tot,sense);
            d.settings('rel_subopt',0.1)
            [x,fval,exitflag,daqp_bnb_info_subopt] = d.solve();
            testCase.verifyEqual(exitflag,int32(1));
            display(daqp_bnb_info_subopt)

            % compare with Gurobi (comparison skipped if Gurobi is not available)
            model.Q = 0.5*sparse(H);
            model.A = sparse([A;-A(me+1:end,:)]);
            model.rhs = [bupper;-blower(me+1:end)];
            model.obj = f;
            model.sense=repelem('=<',[me,2*(m-me)]);
            model.vtype=repelem('BC',[ms,n-ms]);
            model.lb = -inf(n,1);
            params.OutputFlag=0;

            if(exist('gurobi'))
                gurobi_result = gurobi(model,params)
                rel_error = (fval-gurobi_result.objval)/fval;
                fprintf('fval_daqp:%f fval_grb:%f, rel_error:%f',fval,gurobi_result.objval,rel_error);
                testCase.verifyLessThan((fval-gurobi_result.objval)/fval,0.1);
            end

        end
    end
end
