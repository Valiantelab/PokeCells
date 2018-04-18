clc;
clear all; 
close all;

global DATA_DIR;
DATA_DIR = 'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\';

addpath(genpath(DATA_DIR));

excelName = 'sampleIO';

[~, ap] = excel_read(strcat(DATA_DIR, '\Total5\'), excelName);

numCells = length(ap);
for i = 1:numCells
    ap(i).io.pulsestartIdeal = 160; %you overwrite this anyways!!
    ap(i).io.pulsestart = ap(i).io.pulsestartIdeal; %you overwrite this anyways!!
    ap(i).io.pulsedurIdeal = 600;
    ap(i).io.pulsedur = ap(i).io.pulsedurIdeal;
    ap(i).io.pampstartIdeal = -400;
    ap(i).io.pampstepIdeal = 50;
end

%data is acquired as t x 2chan x episodes

sliceIO(ap(2)); % There will be three for three different slices

