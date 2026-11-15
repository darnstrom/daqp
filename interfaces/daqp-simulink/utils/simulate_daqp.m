function out = simulate_daqp(data, params, connected, mode)
% SIMULATE_DAQP Solve a sequence of QPs with the DAQP S-function in Simulink.
%
%   out = simulate_daqp(data, params, connected, mode)
%
%   data:      struct array (one element per time step) with the fields
%              H, f, A, lb, ub, sense (only the connected ones are used)
%   params:    parameters of the S-function, e.g., '1000,1'
%   connected: logical [1 x 6], which of the inputs H, f, A, lb, ub, sense
%              to connect (default: all)
%   mode:      'sim' (default) to simulate the model, or 'codegen' to
%              generate code for the model (with the generic real-time
%              target), build it, and run the resulting executable
%              (in the current folder)
%   out:       struct with the fields x [n x T], lambda [m x T], fval,
%              exitflag, and iter [1 x T]
%
%   A model with From Workspace blocks for the inputs and To Workspace
%   blocks for the outputs is created (and closed) for each call. The data
%   is stored in the model workspace, so that it is part of generated code.

inputs = {'H', 'f', 'A', 'lb', 'ub', 'sense'};
outputs = {'x', 'lambda', 'fval', 'exitflag', 'iter'};
if nargin < 3 || isempty(connected)
    connected = true(1, 6);
end
if nargin < 4
    mode = 'sim';
end
T = numel(data);

model = 'daqp_sfunc_test_model';
if bdIsLoaded(model)
    close_system(model, 0);
end
new_system(model);
cleanup = onCleanup(@() close_system(model, 0));
set_param(model, 'Solver', 'FixedStepDiscrete', 'FixedStep', '1', ...
    'StopTime', num2str(T-1), 'UnconnectedInputMsg', 'none');

add_block('simulink/User-Defined Functions/S-Function', [model '/DAQP'], ...
    'FunctionName', 'DAQP_sfunc', 'Parameters', params);

mws = get_param(model, 'ModelWorkspace');
for i = find(connected)
    name = inputs{i};
    add_block('simulink/Sources/From Workspace', [model '/' name], ...
        'VariableName', [name '_in'], 'SampleTime', '1', ...
        'Interpolate', 'off', 'OutputAfterFinalValue', 'Holding final value');
    add_line(model, [name '/1'], sprintf('DAQP/%d', i));
    assignin(mws, [name '_in'], to_signal({data.(name)}));
end
for i = 1:numel(outputs)
    add_block('simulink/Sinks/To Workspace', [model '/' outputs{i} '_out'], ...
        'VariableName', outputs{i}, 'SaveFormat', 'Array');
    add_line(model, sprintf('DAQP/%d', i), [outputs{i} '_out/1']);
end

switch mode
    case 'sim'
        simOut = sim(model, 'ReturnWorkspaceOutputs', 'on');
        for i = 1:numel(outputs)
            out.(outputs{i}) = simOut.get(outputs{i}).'; % [width x T]
        end
    case 'codegen'
        % The To Workspace blocks are logged (with the prefix rt_) to
        % <model>.mat by the executable
        set_param(model, 'SystemTargetFile', 'grt.tlc', ...
            'GenerateMakefile', 'on', 'MatFileLogging', 'on', ...
            'GenerateReport', 'off');
        slbuild(model);
        exe = fullfile(pwd, model);
        if ispc
            exe = [exe '.exe'];
        end
        [status, msg] = system(['"' exe '"']);
        if status ~= 0
            error('simulate_daqp:codegen', 'Running %s failed:\n%s', exe, msg);
        end
        logged = load(fullfile(pwd, [model '.mat']));
        for i = 1:numel(outputs)
            out.(outputs{i}) = logged.(['rt_' outputs{i}]).'; % [width x T]
        end
    otherwise
        error('simulate_daqp:mode', 'Unknown mode "%s"', mode);
end
end

function s = to_signal(values)
% From Workspace data (structure with time) for a sequence of values
T = numel(values);
values{end+1} = values{end}; % (avoid the ambiguity of a single sample)
s.time = (0:T)';
v = values{1};
if iscolumn(v) % Vector => 1-D signal
    s.signals.values = cell2mat(cellfun(@(x) x(:)', values(:), 'UniformOutput', false));
    s.signals.dimensions = numel(v);
else
    s.signals.values = cat(3, values{:});
    s.signals.dimensions = size(v);
end
end
