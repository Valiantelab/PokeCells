function [ReboundSpikes] = getReboundSpikes(dataPostInjection, analysisParams)
%getReboundSpikes:
%   dataPostInjection is the data after the current injection period
%   returns ReboundSpikes .pks, .locs, and .num (number of spikes)

ReboundSpikes = {};
[~, nepochs] = size(dataPostInjection);

warning('off');

for j=1:nepochs
    ts = squeeze(dataPostInjection(:,j));
    [ReboundSpikes.pks{j},ReboundSpikes.locs{j}] = findpeaks(ts, 'MINPEAKHEIGHT', analysisParams.io.minpeakheight);
    
    if ~isempty(ReboundSpikes.locs{j})
        ReboundSpikes.num{j} = [ReboundSpikes.locs{j}];
    else
        ReboundSpikes.num{j} = [0];
    end
    
    warning('on');
end

