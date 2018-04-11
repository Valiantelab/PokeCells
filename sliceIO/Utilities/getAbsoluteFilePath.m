function [aPath] = getAbsoluteFilePath(fRoot, rpath)

% Get Absolute file path depending on inputs

rPath = java.io.File(rpath);

if isempty(fRoot)
    aPath = rpath;
elseif rPath.isAbsolute
        aPath = rpath;
else
    aPath = fullfile(fRoot, rpath);
end
        
    