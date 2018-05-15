function [Spikes] = findSpikes(data, tvec, analysisParams)

[~, nepochs] = size(data);
Spikes = [];

warning('off');

for j=1:nepochs
    ts = squeeze(data(:,j));
    [Spikes.pks{j},Spikes.locs{j}] = findpeaks(ts, 'MINPEAKHEIGHT', analysisParams.io.minpeakheight);
    
    %check if spikes are within current injection period
    ind = find(tvec(Spikes.locs{j}) > analysisParams.io.pulsestart  & tvec(Spikes.locs{j}) < (analysisParams.io.pulsestart+analysisParams.io.pulsedur));
    
    if ~isempty(ind)
        Spikes.pks{j} = Spikes.pks{j}(ind);
        Spikes.locs{j} = Spikes.locs{j}(ind);
        Spikes.num{j} = [length(ind)];
    else
        Spikes.pks{j} = [];
        Spikes.locs{j} = [];
        Spikes.num{j} = [0];
    end
    
end

warning('on');
