function [] = axes_text_style(curr_ax, fsize, fname)

if nargin < 1; curr_ax = gca; end
if nargin < 2; fsize = 8; end
if nargin < 3; fname = 'Calibri Light'; end

if isempty(curr_ax)
    allAxes = findall(0,'type','axes');
    set(allAxes, 'FontSize', fsize, 'FontName', fname);
    set(allAxes, 'TickDir', 'out');
else
    set(curr_ax, 'FontSize', fsize, 'FontName', fname);
    set(curr_ax, 'TickDir', 'out');
end
