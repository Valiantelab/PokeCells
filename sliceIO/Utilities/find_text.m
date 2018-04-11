function [cindex] = find_text(text_cell, str, findSubStr)

% Finds the index of a specific string in a cell array of strings

if nargin < 3; findSubStr = 0; end

cindex = 0;

if ~findSubStr
    for i=1:numel(text_cell)
        if strcmpi(text_cell{i},str)
            cindex = i;
            break;
        end
    end
else
    for i=1:numel(text_cell)
        if strfind(text_cell{i},str)
            cindex = i;
            break;
        end
    end
end