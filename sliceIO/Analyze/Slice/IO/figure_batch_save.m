function [] = figure_batch_save(figs, eDir, dosave)
% Save for making figures
if dosave
    display('Saving all figures...')
    for i=1:numel(figs)
        save_figure(figs(i).h, eDir, figs(i).name);
    end
end