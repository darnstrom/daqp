function results = run_all_tests(varargin)
% Compile the S-function and run the Simulink tests
%   run_all_tests            - simulation tests
%   run_all_tests('codegen') - also test code generation (requires Simulink Coder)
%   run_all_tests('example') - also run the example (examples/daqp_example.m)
this_dir = fileparts(mfilename('fullpath'));
make_sfunc();
tests = {fullfile(this_dir, 'test', 'core_test.m')};
if any(strcmp(varargin, 'codegen'))
    tests{end+1} = fullfile(this_dir, 'test', 'codegen_test.m');
end
results = runtests(tests);
disp(table(results));
assertSuccess(results);
if any(strcmp(varargin, 'example'))
    addpath(fullfile(this_dir, 'examples'));
    daqp_example();
end
end
