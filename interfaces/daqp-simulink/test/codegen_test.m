classdef codegen_test < matlab.unittest.TestCase
    % Tests of code generation (Simulink Coder, generic real-time target)
    % for the DAQP S-function: the generated executable should give the
    % same results as the simulation

    methods (TestClassSetup)
        function setup(testCase)
            here = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, '..')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, '..', 'utils')));
            % Build in a temporary folder
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
        end
    end

    methods (Test)
        function warm_start(testCase)
            rng(1);
            T = 10; n = 10; m = 20; ms = 5;
            [~,H,f,A,ub,lb] = generate_test_QP(n,m,ms,5,1e2);
            for t = 1:T
                data(t) = struct('H',H,'f',f*cos(0.1*(t-1)),'A',A, ...
                    'lb',lb,'ub',ub,'sense',zeros(m,1)); %#ok<AGROW>
            end
            for warm_start = 0:2
                params = sprintf('1000,%d', warm_start);
                testCase.verifySameAsSim(data, params, true(1,6));
            end
        end

        function unconnected_inputs(testCase)
            % Only H, f, A, and ub
            d = struct('H',2*eye(3),'f',-4*ones(3,1),'A',ones(1,3),'lb',[],'ub',1,'sense',[]);
            out = testCase.verifySameAsSim(d, '1000', [true true true false true false]);
            testCase.verifyEqual(out.x, ones(3,1)/3, 'AbsTol', 1e-9);
        end

        function hierarchical(testCase)
            % From https://github.com/darnstrom/daqp/issues/161
            A = [1 -1 0 -5 -5 0 0 0 0; -1 1 0 0 0 0 0 0 0; -0.5 -0.5 0 0 0 0 0 0 0];
            lb = [0 0 0 -1.8*ones(1,6) 0 0 0]';
            ub = [4 4 4 1.8*ones(1,6) 0 0 0]';
            rhs = [0 20 30 10];
            for t = 1:numel(rhs)
                lb(10) = rhs(t); ub(10) = rhs(t);
                data(t) = struct('H',eye(9),'f',zeros(9,1),'A',A,'lb',lb, ...
                    'ub',ub,'sense',[zeros(9,1); 5*ones(3,1)]); %#ok<AGROW>
            end
            out = testCase.verifySameAsSim(data, '1000,1,[9 10 11 12]', true(1,6));
            testCase.verifyEqual(A(1,:)*out.x, max(-22, min(22, rhs)), 'AbsTol', 1e-6);
        end
    end

    methods
        function out = verifySameAsSim(testCase, data, params, connected)
            ref = simulate_daqp(data, params, connected, 'sim');
            out = simulate_daqp(data, params, connected, 'codegen');
            testCase.verifyEqual(out.exitflag, ref.exitflag);
            testCase.verifyGreaterThan(out.exitflag, 0);
            testCase.verifyEqual(out.iter, ref.iter);
            testCase.verifyEqual(out.x, ref.x, 'AbsTol', 1e-9);
            testCase.verifyEqual(out.lambda, ref.lambda, 'AbsTol', 1e-9);
            testCase.verifyEqual(out.fval, ref.fval, 'AbsTol', 1e-9);
        end
    end
end
