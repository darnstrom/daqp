classdef core_test < matlab.unittest.TestCase
    % Sovle 1000 randomly generated QPs
    methods (Test)
        function random_feasible_QPs(testCase)
            % Test on randomly generated feasible QPs
            rng('default');
            nQPs = 20;
            n = 5; m = 10; ms = 5;
            nAct = 3;
            kappa = 1e2;
            tol = 1e-5;
            
            % init Vectors
            primal_errors = zeros(nQPs,1);
            dual_errors = zeros(nQPs,1);
            f_val_errors = zeros(nQPs,1);
            
            % Init QP
            xref = cell(nQPs,1);
            H = cell(nQPs,1);
            f = cell(nQPs,1);
            A = cell(nQPs,1);
            bupper = cell(nQPs,1);
            blower = cell(nQPs,1);
            sense = cell(nQPs,1);
            for i = 1:nQPs
                [xref{i},H{i},f{i},A{i},bupper{i},blower{i},sense{i}] = generate_test_QP(n,m,ms,nAct,kappa);
            end
            
            % Call DAQP in Simulink with each input combination 
            inputs = struct('H',H,'f',f,'A',A,'bupper',bupper,'blower',blower,'sense',sense);
            data = simulate_model('DAQP_simulink_test',inputs);
            
            % Calculate the Model Errors
            for i = 1:nQPs
                primal_errors(i) = norm(data(i).x-xref{i});
                dual_errors(i) = norm(H{i}*data(i).x+f{i}+[eye(ms,n);A{i}]'*data(i).lambda);
                f_val_errors(i) = norm(0.5*data(i).x'*H{i}*data(i).x+f{i}'*data(i).x-data(i).fval);
            end
            
            % Test            
            testCase.verifyEqual([data.exitflag],ones(1,nQPs));
            testCase.verifyLessThan(primal_errors,tol);
            testCase.verifyLessThan(dual_errors,tol);
            testCase.verifyLessThan(f_val_errors,tol);
            
            % Print
            fprintf('========================== DAQP =============================\n')
            fprintf('Primal Errors  : |avg: %2.2e| max: %2.2e| min %2.2e|\n',mean(primal_errors),max(primal_errors),min(primal_errors))
            fprintf('Dual Errors    : |avg: %2.2e| max: %2.2e| min %2.2e|\n',mean(dual_errors),max(dual_errors),min(dual_errors))
            fprintf('Fval Errors    : |avg: %2.2e| max: %2.2e| min %2.2e|\n',mean(f_val_errors),max(f_val_errors),min(f_val_errors))
            fprintf('=============================================================\n')
        end
    end
end
