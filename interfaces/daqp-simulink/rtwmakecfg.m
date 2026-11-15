function makeInfo = rtwmakecfg()
% RTWMAKECFG Add include and source directories to RTW make files.
% makeInfo = RTWMAKECFG returns a structured array containing following fields:

fprintf('DAQP: Adding source and include directories to make process.\n')

% get the daqp root and src/include directories
sfuncDir = fileparts(mfilename('fullpath'));
daqpRootDir = fileparts(fileparts(sfuncDir)); % get the "../../" path
daqpSourcesDir = fullfile(daqpRootDir, 'src');
daqpIncludeDir = fullfile(daqpRootDir, 'include');
daqpSources = dir(fullfile(daqpSourcesDir, '*.c'));

% Setting up the return structure with
makeInfo.sourcePath = {daqpSourcesDir, sfuncDir};
makeInfo.includePath = {daqpIncludeDir, sfuncDir};
makeInfo.sources = [{daqpSources.name}, {'daqp_sfunc_core.c'}];

% Display contents.
fprintf(' - additional source directories:\n');
fprintf('   - "%s"\n', makeInfo.sourcePath{:});
fprintf(' - additional include directories:\n');
fprintf('   - "%s"\n', makeInfo.includePath{:});
fprintf(' - additional source files:\n');
fprintf('   - "%s"\n', makeInfo.sources{:});

end
