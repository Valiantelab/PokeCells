% Part 1: AHP for layer2/3 comparison
% Part 2: Active & Passive properties DC injection
% Part 3: Ramp

%% Part 1: Load in all of the data
clc;
clear all; 
close all;

global DATA_DIR;
DATA_DIR = 'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\AHPMay7\';
CODE_DIR =  'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\';
addpath(genpath(CODE_DIR));

excelName = 'AHP05102018';

[~, ap] = excel_read(DATA_DIR, excelName);

numCells = length(ap);
for i = 1:numCells
    if strcmp(ap(i).Dir, 'Total_5_Homeira')
        %Homeira's protocol
        %disp('---- Homeira Protocol ----')
        ap(i).io.pampstartIdeal = -400;
        ap(i).io.pampstepIdeal = 50;
        
        ap(i).io.pulsestartIdeal = 160; 
        ap(i).io.pulsedurIdeal = 600;

    %elseif strcmp(ap(i).Dir, 'Total2Homeira')
     %   continue
        %too many different protocols...
    elseif strcmp(ap(i).Dir, 'TotaL_5_Lihua') || strcmp(ap(i).Dir, 'Total2lihua')
        %Lihua
        %disp('---- Lihua Protocol ----')
        ap(i).io.pampstartIdeal = -400;
        ap(i).io.pampstepIdeal = 25;
        
        ap(i).io.pulsestartIdeal = 65; 
        ap(i).io.pulsedurIdeal = 600; 
    else
        disp(ap(i).Dir)
        error('Silly goose, you forgot this directory');
    end
    
    ap(i).io.pulsestart = ap(i).io.pulsestartIdeal; 
    ap(i).io.pulsedur = ap(i).io.pulsedurIdeal;
    ap(i).io.pampstart = ap(i).io.pampstartIdeal;
    ap(i).io.pampstep = ap(i).io.pampstepIdeal;
end

% Run the above, then save workspace if your little heart so desires
%%
%data is acquired as t x 2chan x episodes

close all;
badCell = zeros(1, numCells);
startCell = 154;
numPruned = 18;

for i = startCell:numCells
    disp('Now analyzing slice..')
    disp(i)

    io(i-numPruned) = sliceIO(ap(i));
%     hold = input('Is everything good? Blank for YES, anything else for no');
%     
%     if ~strcmp(hold, '')
%         badCell(i) = 1;
%     end
%     
    close all;
end

%% Print to excel file
% AHP only
xlrow = 1; %track which row we are on

outputxlname = 'AHPvalues05102018_R4.xlsx';

xltgt =  strcat('A', num2str(xlrow));
    
headings = {'Data Dir', 'File name', 'Layer', 'Current injection (pA)', 'AHP AUC (ms * mV)', 'AHP min (mV)', 'Cell#', 'Episode#'};
xlswrite(outputxlname, headings, 1, xltgt);

numCells = length(io);
for i = 1:numCells
    currStart = ap(i).io.pampstartIdeal;
    currStep = ap(i).io.pampstepIdeal;
    currEnd = currStart + currStep * (io(i).numEpisodes{1} - 1);
    
    currLevels = currStart:currStep:currEnd;
    
    if length(currLevels) ~= io(i).numEpisodes{1}
        error('You messed up!!!1!!1! >:) ')
    end
    if strcmp(ap(i).Dir, 'Total_5_Homeira') || strcmp(ap(i).Dir, 'TotaL_5_Lihua')
        layer = 5;
    elseif  strcmp(ap(i).Dir, 'Total2Homeira') || strcmp(ap(i).Dir, 'Total2lihua')
        layer = 2;
    end
    
    for j = 1:io(i).numEpisodes{1}
        if ~isempty(io(i).ahp{j}.auc)
            xlrow = xlrow + 1;
            dataArr = { ap(i).Dir,  ap(i).fname, layer, currLevels(j), io(i).ahp{j}.auc{1}, io(i).ahp{j}.min{1}, i, j};
        %else
        %    dataArr = { ap(i).Dir,  ap(i).fname, layer, currLevels(j), 'None', 'None', i, j};
        
            xltgt =  strcat('A', num2str(xlrow));
            xlswrite(outputxlname, dataArr, 1, xltgt);
        end
    end
    disp(i)
end
%% Print to excel file - active and passive props
xlrow = 1; %track which row we are on

outputxlname = 'ActivePassiveProps.xlsx';

xltgt =  strcat('A', num2str(xlrow));
    
headings = {'Data Dir', 'File name', 'Layer', 'Current injection (pA)', 'AHP AUC (ms * mV)', 'AHP min (mV)', 'Cell#', 'Episode#'};
xlswrite(outputxlname, headings, 1, xltgt);

for i = 1:numCells
    currStart = ap(i).io.pampstartIdeal;
    currStep = ap(i).io.pampstepIdeal;
    currEnd = currStart + currStep * (io(i).numEpisodes{1} - 1);
    
    currLevels = currStart:currStep:currEnd;
    
    if length(currLevels) ~= io(i).numEpisodes{1}
        error('You messed up!!!1!!1! >:) ')
    end
    if strcmp(ap(i).Dir, 'Total_5_Homeira') || strcmp(ap(i).Dir, 'TotaL_5_Lihua')
        layer = 5;
    elseif  strcmp(ap(i).Dir, 'Total2Homeira') || strcmp(ap(i).Dir, 'Total2lihua')
        layer = 2;
    end
    
    for j = 1:io(i).numEpisodes{1}
        if ~isempty(io(i).ahp{j}.auc)
            xlrow = xlrow + 1;
            dataArr = { ap(i).Dir,  ap(i).fname, layer, currLevels(j), io(i).ahp{j}.auc{1}, io(i).ahp{j}.min{1}, i, j};
            xltgt =  strcat('A', num2str(xlrow));
            xlswrite(outputxlname, dataArr, 1, xltgt);
        end
    end
end