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
    %this is a bandaid solution that I've applied
    %need a cure
    ap(i).io.pulsedurIdeal = 600;
    ap(i).io.pampstartIdeal = -400;
    ap(i).io.pampstepIdeal = 50;
    ap(i).io.pulsestartIdeal = 160; 
    
    ap(i).io.pulsestart = ap(i).io.pulsestartIdeal; 
    ap(i).io.pulsedur = ap(i).io.pulsedurIdeal;
    ap(i).io.pampstart = ap(i).io.pampstartIdeal;
    ap(i).io.pampstep = ap(i).io.pampstepIdeal;
end

%data is acquired as t x 2chan x episodes
for i = 7:7
    sliceIO(ap(i)); % There will be three for three different slices
end
