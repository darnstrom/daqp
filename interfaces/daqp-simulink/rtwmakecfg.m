function makeInfo = rtwmakecfg()
% RTWMAKECFG Add include and source directories to RTW make files.
% makeInfo = RTWMAKECFG returns a structured array containing following fields:


fprintf('DAQP: Adding source and include directories to make process.\n')


% get the daqp root and src/include directories
daqpRootDir = fileparts(fileparts(fileparts(mfilename('fullpath')))); % get the "../../" path
daqpSourcesDir = fullfile(daqpRootDir, 'src');
daqpIncludeDir = fullfile(daqpRootDir, 'include');

% Setting up the return structure with
makeInfo.sourcePath = {daqpSourcesDir};
makeInfo.includePath = {daqpSourcesDir, daqpIncludeDir};
makeInfo.sources = {'daqp.c', 'factorization.c', 'utils.c' 'api.c', 'auxiliary.c', 'daqp_prox.c', 'bnb.c'};

% Display contents.
fprintf(' - additional source directories:\n');
fprintf('   - "%s"\n', makeInfo.sourcePath{:});
fprintf(' - additional include directories:\n');
fprintf('   - "%s"\n', makeInfo.includePath{:});
fprintf(' - additional source files:\n');
fprintf('   - "%s"\n', makeInfo.sources{:});

end