function [ss, analysisParams] = excel_read(Dir, fName)

% Function to read in a excel spreadsheet and parse the data out of it.
% The filename field must be specified to be a text column or this will not
% work.  The comment field is optional.  The excel file must follow a
% specific convention:

%  'Dir'          'Filename'    'Cond1'    'Cond1'    'Cond2'    'Cond2' ... CH to ANALYZE  'SFs'(optional #)...
%                                start        end       start      end
%         

%ss: spreadsheet
%sf: special field

sf_names = {'Comment', 'ExportDir', 'Tags', 'GroupInc'};

full_path = fullfile(Dir,[fName '.xlsx']);

if ~exist(full_path, 'file')
    display(full_path);
    error('Unable to open specified file');
end

%% Load the files and their condtions etc etc

% Read in SHEET 1
display(sprintf('EXCEL_READ: Opening %s',fullfile(Dir, fName)));
[ss.num, ss.txt, ss.raw] = xlsread(full_path, 1);

% Read in SHEET 2 - channel assignments
[ch.num, ch.txt, ch.raw] = xlsread(full_path, 2);

if numel(ss.txt) < 3
    % There must be more than just a Dir and fname field
    error('Not enough fields.');
end

% Get number of conditions and files
nfiles = size(ss.txt,1)-1;
ncond = (size(ss.txt,2)-6)/3; % I guess 6 is hard coded

% Get the condition names
for i=1:ncond
    cond_names{i} = ss.txt{1,3*i}; %start/end/anchan for each
end

% Find any special fields
nsf = numel(sf_names);
disp(nsf)
for i=1:nsf
    sf_index(i) = find_text( {ss.txt{1,:}}, sf_names{i});
end   

% Populate a cell of file parameters for each of the files

if max(cellfun(@isempty,  {ss.txt{2:end,2}})) == 1
    % Not all filenames filled, there are blanks
    error('Error in reading excel spreadsheet, filename fields are not all completely text');
end

for i=1:nfiles
    analysisParams(i) = setAnalysisParams();
    analysisParams(i).Dir = ss.txt{i+1,1};
    analysisParams(i).fname = ss.txt{i+1,2};
    analysisParams(i).cond.names = {};
    
    for j=1:ncond
        t_limits = ss.num(i,(3*j-2):(3*j-1));
        analysisParams(i).cond.times(:,j) = t_limits;
        analysisParams(i).cond.names{j} = cond_names{j};
        analysisParams(i).cond.fname{j} = ss.txt{i+1,3*j-1}(1:end-1);
    end
    
    % Get the special fields
    c = 0;
    analysisParams(i).sf.names = [];  % Clear the names
    for j=1:nsf
        if sf_index(j)
            c = c + 1;
            analysisParams(i).sf.names{c} = sf_names{j};
            analysisParams(i).sf.vals{c} = ss.txt{i+1,sf_index(j)};
        end
    end
    % Set the channel to analyze which must be just after the condition
    % timings
    analysisParams(i).ch = ss.num(i,3*ncond);
    
    % Get the channel labels from the table.  If there are multiple files
    % for the same slice it is assumed that the channel assigments did not
    % change during the experiment
    
    for l=1:(size(ch.txt,1)-1)
        chfnames{l} = ch.txt{l+1,1}(1:end-1);
    end
    
    rindex = find_text(chfnames, analysisParams(i).cond.fname{1});
    if ~rindex
        display(sprintf('%s has no channel assigment.', analysisParams(i).cond.fname{1}));
    else
        analysisParams(i).chlabels = {ch.txt{rindex+1,2:end}};
    end
        
end