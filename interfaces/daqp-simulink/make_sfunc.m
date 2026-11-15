function make_sfunc(do_debug)
%-------------------------------------------------------------------------%
%
%   Description : MEX-compilation of the daqp solver for simulink with code
%                 generation.
%   Author      : Christopher Schulte, RWTH Aachen University, IRT
%                 christopher.schulte@rwth-aachen.de
%
%-------------------------------------------------------------------------%

if nargin == 0
    do_debug = false;
end

% Get path of this file and the daqp sources
this_dir  = fileparts(mfilename('fullpath'));
daqp_root = fullfile(this_dir, '..', '..');
include_dir = fullfile(daqp_root, 'include');
source_files = dir(fullfile(daqp_root, 'src', '*.c'));
source_files = fullfile({source_files.folder}, {source_files.name});

args = {'-outdir', this_dir, ['-I' include_dir], ['-I' this_dir]};
if do_debug
    args = [args, {'-g'}];
end

mex(args{:}, fullfile(this_dir, 'DAQP_sfunc.c'), ...
    fullfile(this_dir, 'daqp_sfunc_core.c'), source_files{:});

end
