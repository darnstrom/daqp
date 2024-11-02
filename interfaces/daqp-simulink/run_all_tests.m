
% Change the path to the folder where this file is located
cd(fileparts(mfilename('fullpath')))
addpath('test')
addpath('utils')

% Compile the s-function
doDebug = false;
make_sfunc(doDebug);

runtests('core_test')