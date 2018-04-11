function [fcount, figs] = figure_set(fcount, figs, fname)
% Keeps track on new figures and their handles
fcount = fcount + 1;
figs(fcount).h = figure(fcount);
figs(fcount).name = fname;
set(figs(fcount).h, 'Name', figs(fcount).name)
clf;