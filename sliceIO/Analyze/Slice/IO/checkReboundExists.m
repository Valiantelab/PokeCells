function [ReboundSpikes] = getReboundSpikes(dataPostInjection, analysisParams)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

ReboundSpikes = {};

[~, nepochs] = size(dataPostInjection);

warning('off');

for j=1:nepochs
    ts = squeeze(dataPostInjection(:,j));
    [ReboundSpikes.pks{j},ReboundSpikes.locs{j}] = findpeaks(ts, 'MINPEAKHEIGHT', analysisParams.io.minpeakheight);
    
    if ~isempty(ind)
        ReboundSpikes.pks{j} = ReboundSpikes.pks{j}(ind);
        ReboundSpikes.locs{j} = ReboundSpikes.locs{j}(ind);
        ReboundSpikes.num{j} = [length(ind)];
    else
        ReboundSpikes.pks{j} = [];
        ReboundSpikes.locs{j} = [];
        ReboundSpikes.num{j} = [0];
    end
    
    warning('on');
end

