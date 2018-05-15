%% Part 1: Load in all of the data
clc;
clear all; 
close all;

global DATA_DIR;
DATA_DIR = 'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\C3202018\';
CODE_DIR =  'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\';
addpath(genpath(CODE_DIR));

excelName = 'RAMP05152018';

[~, ap] = excel_read(DATA_DIR, excelName);
numCells = length(ap);

for i = 1:1
    io(i) = rampAnalysis(ap(i));
end

%% for printing
for i = 1:numCells
    if strcmp(ap(i).Dir, 'Case1\Ramp')
    elseif strcmp(ap(i).Dir, 'Case2\Ramp')
    else
        error('You''re soo dumb, wrong directory string!');
    end
end