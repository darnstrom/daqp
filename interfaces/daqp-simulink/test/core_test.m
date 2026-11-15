classdef core_test < matlab.unittest.TestCase
    % Tests of the DAQP S-function in Simulink

    methods (TestClassSetup)
        function addUtils(testCase)
            here = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, '..')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, '..', 'utils')));
        end
    end

    methods (Test)
        function random_feasible_QPs(testCase)
            % A different random QP at every time step
            rng('default');
            T = 20; n = 50; m = 40; ms = 10; nAct = 10; kappa = 1e2;
            tol = 1e-5;
            for t = 1:T
                [xref{t},data(t).H,data(t).f,data(t).A,data(t).ub,data(t).lb,sense] = ...
                    generate_test_QP(n,m,ms,nAct,kappa); %#ok<AGROW>
                data(t).sense = double(sense);
            end
            for warm_start = 0:1
                out = simulate_daqp(data, sprintf('1000,%d', warm_start));
                testCase.verifyEqual(out.exitflag, ones(1,T));
                for t = 1:T
                    x = out.x(:,t); lam = out.lambda(:,t);
                    d = data(t);
                    testCase.verifyLessThan(norm(x-xref{t}), tol);
                    testCase.verifyLessThan(norm(d.H*x+d.f+[eye(ms,n);d.A]'*lam), tol);
                    testCase.verifyLessThan(abs(0.5*x'*d.H*x+d.f'*x-out.fval(t)), tol);
                end
            end
        end

        function warm_start(testCase)
            % Slowly varying linear term => warm starting should reduce
            % the number of iterations, without changing the solution
            rng(1);
            T = 30; n = 30; m = 50; ms = 10;
            [~,H,f,A,ub,lb] = generate_test_QP(n,m,ms,15,1e2);
            for t = 1:T
                data(t) = struct('H',H,'f',f*cos(0.05*(t-1)),'A',A, ...
                    'lb',lb,'ub',ub,'sense',zeros(m,1)); %#ok<AGROW>
            end
            cold = simulate_daqp(data, '1000,0');
            warm = simulate_daqp(data, '1000,1');
            fixed = simulate_daqp(data, '1000,2');
            testCase.verifyEqual(cold.exitflag, ones(1,T));
            testCase.verifyEqual(warm.exitflag, ones(1,T));
            testCase.verifyEqual(fixed.exitflag, ones(1,T));
            testCase.verifyEqual(warm.x, cold.x, 'AbsTol', 1e-6);
            testCase.verifyEqual(fixed.x, cold.x, 'AbsTol', 1e-6);
            testCase.verifyEqual(warm.lambda, cold.lambda, 'AbsTol', 1e-6);
            testCase.verifyLessThan(sum(warm.iter), sum(cold.iter));
            testCase.verifyLessThan(sum(fixed.iter), sum(cold.iter));
        end

        function unconnected_inputs(testCase)
            H = 2*eye(3); f = -4*ones(3,1);
            % Only A and ub: x1+x2+x3 <= 1
            d = struct('H',H,'f',f,'A',ones(1,3),'lb',[],'ub',1,'sense',[]);
            out = simulate_daqp(d, '1000', [true true true false true false]);
            testCase.verifyEqual(out.exitflag, 1);
            testCase.verifyEqual(out.x, ones(3,1)/3, 'AbsTol', 1e-9);
            testCase.verifyEqual(out.lambda, 10/3, 'AbsTol', 1e-9);
            % Only simple bounds: -1 <= x <= 1
            d = struct('H',H,'f',[-4;0;4],'A',[],'lb',-ones(3,1),'ub',ones(3,1),'sense',[]);
            out = simulate_daqp(d, '1000', [true true false true true false]);
            testCase.verifyEqual(out.exitflag, 1);
            testCase.verifyEqual(out.x, [1;0;-1], 'AbsTol', 1e-9);
            testCase.verifyEqual(out.lambda, [2;0;-2], 'AbsTol', 1e-9);
            % No constraints
            d = struct('H',H,'f',[-4;0;4],'A',[],'lb',[],'ub',[],'sense',[]);
            out = simulate_daqp(d, '1000', [true true false false false false]);
            testCase.verifyEqual(out.exitflag, 1);
            testCase.verifyEqual(out.x, [2;0;-2], 'AbsTol', 1e-9);
        end

        function equality_and_infeasible(testCase)
            % x1+x2 = 1 (sense = 5), then x1+x2 >= 3 with |x| <= 1
            % (infeasible), and then back to the equality constraint
            H = eye(2); f = zeros(2,1); A = [1 1];
            feas = struct('H',H,'f',f,'A',A,'lb',[-1;-1;1],'ub',[1;1;1],'sense',[0;0;5]);
            infeas = struct('H',H,'f',f,'A',A,'lb',[-1;-1;3],'ub',[1;1;4],'sense',[0;0;0]);
            data = [feas, infeas, feas];
            for warm_start = 0:2
                out = simulate_daqp(data, sprintf('1000,%d', warm_start));
                testCase.verifyEqual(out.exitflag, [1 -1 1]);
                testCase.verifyEqual(out.x(:,[1 3]), 0.5*ones(2,2), 'AbsTol', 1e-9);
            end
        end

        function hierarchical(testCase)
            % From https://github.com/darnstrom/daqp/issues/161
            % Simple bounds (hard) and three equality constraints with
            % decreasing priority
            A = [1 -1 0 -5 -5 0 0 0 0; -1 1 0 0 0 0 0 0 0; -0.5 -0.5 0 0 0 0 0 0 0];
            lb = [0 0 0 -1.8*ones(1,6) 0 0 0]';
            ub = [4 4 4 1.8*ones(1,6) 0 0 0]';
            rhs = [0 10 20 26 30 26 10 0 -15 26];
            for t = 1:numel(rhs)
                lb(10) = rhs(t); ub(10) = rhs(t);
                data(t) = struct('H',eye(9),'f',zeros(9,1),'A',A,'lb',lb, ...
                    'ub',ub,'sense',[zeros(9,1); 5*ones(3,1)]); %#ok<AGROW>
            end
            cold = simulate_daqp(data, '1000,0,[9 10 11 12]');
            warm = simulate_daqp(data, '1000,1,[9 10 11 12]');
            for out = [cold, warm]
                testCase.verifyGreaterThan(out.exitflag, 0);
                Ax = A*out.x;
                % The first level can reach [-22, 22]
                testCase.verifyEqual(Ax(1,:), max(-22, min(22, rhs)), 'AbsTol', 1e-6);
                % ... and when it cannot be fulfilled, the lower levels are violated
                testCase.verifyEqual(Ax(2:3, rhs > 22), repmat([-4; -2], 1, sum(rhs > 22)), 'AbsTol', 1e-6);
            end
            testCase.verifyEqual(warm.x, cold.x, 'AbsTol', 1e-6);
        end

        function rtwmakecfg_sources(testCase)
            % The build of generated code should include all DAQP sources
            here = fileparts(mfilename('fullpath'));
            info = rtwmakecfg();
            srcs = dir(fullfile(here, '..', '..', '..', 'src', '*.c'));
            testCase.verifyEqual(sort(info.sources), ...
                sort([{srcs.name}, {'daqp_sfunc_core.c'}]));
            for p = [info.sourcePath, info.includePath]
                testCase.verifyTrue(isfolder(p{1}), p{1});
            end
            for f = info.sources
                found = cellfun(@(d) isfile(fullfile(d, f{1})), info.sourcePath);
                testCase.verifyTrue(any(found), f{1});
            end
        end
    end
end
