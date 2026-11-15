function simOut = daqp_example(save_model)
% DAQP_EXAMPLE Example of using the DAQP S-function in Simulink.
%
%   simOut = daqp_example() builds the model daqp_example, opens it (if
%   MATLAB runs with a desktop), and simulates it.
%   daqp_example(true) also saves the model as daqp_example.slx in the
%   current folder.
%
%   At each time step, the QP
%       min  ||x - r(t)||^2  =  x'*x - 2 r(t)'*x  (+ constant)
%       s.t. -1 <= x1 <= 1, -1 <= x2 <= 1   (simple bounds)
%            x1 + x2 <= 1.2                  (general constraint)
%   is solved, i.e., x is the point in the feasible set that is closest to
%   the reference r(t), which moves along a circle with radius 1.5.
%   In terms of the inputs of the DAQP block:
%       H = 2*eye(2), f = -2*r(t), A = [1 1],
%       lb = [-1; -1; -inf], ub = [1; 1; 1.2]
%   The sense input is not connected (=> only inequality constraints).
%   The block parameters are "100, 1": at most 100 iterations, and warm
%   start (the working set of the previous time step is reused).
%
%   The S-function must be compiled (make_sfunc) and on the MATLAB path.

if nargin < 1
    save_model = false;
end
if exist('DAQP_sfunc', 'file') ~= 3
    error('daqp_example:sfunc', ['The S-function DAQP_sfunc was not found. ' ...
        'Compile it with make_sfunc and add interfaces/daqp-simulink to the path.']);
end

model = 'daqp_example';
if bdIsLoaded(model)
    close_system(model, 0);
end
new_system(model);
set_param(model, 'Solver', 'FixedStepDiscrete', 'FixedStep', '0.05', ...
    'StopTime', '10', 'UnconnectedInputMsg', 'none');

% The solver
add_block('simulink/User-Defined Functions/S-Function', [model '/DAQP'], ...
    'FunctionName', 'DAQP_sfunc', 'Parameters', '100, 1', ...
    'Position', [400 100 520 300]);

% Inputs: H, f = -2*r(t), A, lb, ub (sense is not connected)
add_block('simulink/Sources/Constant', [model '/H'], 'Value', '2*eye(2)');
add_block('simulink/Sources/Sine Wave', [model '/r(t)'], ...
    'Amplitude', '1.5', 'Frequency', '1', 'Phase', '[pi/2; 0]', ...
    'SampleTime', '0.05');
add_block('simulink/Math Operations/Gain', [model '/Gain'], 'Gain', '-2');
add_block('simulink/Sources/Constant', [model '/A'], 'Value', '[1 1]', ...
    'VectorParams1D', 'off');
add_block('simulink/Sources/Constant', [model '/lb'], 'Value', '[-1; -1; -inf]');
add_block('simulink/Sources/Constant', [model '/ub'], 'Value', '[1; 1; 1.2]');
add_line(model, 'H/1', 'DAQP/1');
add_line(model, 'r(t)/1', 'Gain/1');
add_line(model, 'Gain/1', 'DAQP/2');
add_line(model, 'A/1', 'DAQP/3');
add_line(model, 'lb/1', 'DAQP/4');
add_line(model, 'ub/1', 'DAQP/5');

% Outputs: x and r(t) in a scope, exitflag and iterations in displays,
% and x, exitflag, and iterations logged to the workspace
add_block('simulink/Sinks/Scope', [model '/x and r(t)'], 'NumInputPorts', '2');
add_line(model, 'DAQP/1', 'x and r(t)/1');
add_line(model, 'r(t)/1', 'x and r(t)/2');
add_block('simulink/Sinks/Terminator', [model '/lambda']);
add_line(model, 'DAQP/2', 'lambda/1');
add_block('simulink/Sinks/Terminator', [model '/fval']);
add_line(model, 'DAQP/3', 'fval/1');
add_block('simulink/Sinks/Display', [model '/exitflag']);
add_line(model, 'DAQP/4', 'exitflag/1');
add_block('simulink/Sinks/Display', [model '/iterations']);
add_line(model, 'DAQP/5', 'iterations/1');
logged = {'x', 1; 'exitflag', 4; 'iter', 5};
for i = 1:size(logged, 1)
    block = [model '/log_' logged{i,1}];
    add_block('simulink/Sinks/To Workspace', block, ...
        'VariableName', logged{i,1}, 'SaveFormat', 'Array');
    add_line(model, sprintf('DAQP/%d', logged{i,2}), ['log_' logged{i,1} '/1']);
end
Simulink.BlockDiagram.arrangeSystem(model);

if usejava('desktop')
    open_system(model);
end
if save_model
    save_system(model, fullfile(pwd, [model '.slx']));
end

% Simulate and check the result
simOut = sim(model, 'ReturnWorkspaceOutputs', 'on');
x = simOut.get('x');
exitflag = simOut.get('exitflag');
iter = simOut.get('iter');
if any(exitflag <= 0) || any(abs(x(:)) > 1+1e-9) || any(sum(x, 2) > 1.2+1e-9)
    error('daqp_example:result', 'Unexpected result of the example.');
end
fprintf('DAQP example: %d QPs solved, %.2f iterations per QP on average.\n', ...
    numel(exitflag), mean(iter));
end
