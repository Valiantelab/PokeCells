% todo: 1) Check for onechannel setting DONE
%       2) Check current levels (different conditions) DONE
%       3) Print to Excel

clc;
clear all; 
close all;

global DATA_DIR;
DATA_DIR = 'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\AHP\';
CODE_DIR =  'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\';
addpath(genpath(CODE_DIR));

excelName = 'AHP';

[~, ap] = excel_read(DATA_DIR, excelName);

numCells = length(ap);
for i = 1:numCells
    if strcmp(ap(i).Dir, 'Total_5_Homeira') || strcmp(ap(i).Dir, 'Total2Homeira')
        %Homeira's protocol
        disp('---- Homeira Protocol ----')
        ap(i).io.pampstartIdeal = -400;
        ap(i).io.pampstepIdeal = 50;
        
        ap(i).io.pulsestartIdeal = 160; 
        ap(i).io.pulsedurIdeal = 600;
        %13o21027 - 11th cell, different start and duration
    elseif strcmp(ap(i).Dir, 'TotaL_5_Lihua') || strcmp(ap(i).Dir, 'Total2lihua')
        %Lihua
        disp('---- Lihua Protocol ----')
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
%skip = [];
for i = 11:numCells
 %   if isempty(find(skip == i, 1))
        disp('Now analyzing slice..')
        disp(i)

        io(i) = sliceIO(ap(i));
        hold = input('Take a look at stuff, then press enter'); 
        close all;
  %  end
end

%% Print to excel file
xlrow = 1; %track which row we are on

outputxlname = 'AHP_values.xlsx';

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
