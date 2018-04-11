clc;
clear all; 
close all;

global DATA_DIR;
DATA_DIR = 'C:\Users\Saima\Documents\REL Projects\Homeira Heterogeneity\';

addpath(genpath(DATA_DIR));

[~, ap] = excel_read(strcat(DATA_DIR, '\Total5\'), 'sampleIO');
sliceIO(ap(3)); % There will be three for three different slices

