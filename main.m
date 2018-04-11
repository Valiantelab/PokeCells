clc;
clear all; 
close all;

global DATA_DIR;
DATA_DIR = 'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\';

addpath(genpath(DATA_DIR));

excelName = 'sampleIO';
[~, ap] = excel_read(strcat(DATA_DIR, '\Total5\'), excelName);

for i = 1:length(ap)
    ap(i).io.pulsestart = 160;
    ap(i).io.pulsedur = 600;
    ap(i).io.pampstart = -400;
    ap(i).io.pampstep = 50;
end

sliceIO(ap(1)); % There will be three for three different slices

