function [F, S, SR] = collect_features(R, ap)
% Collects the various features of spikes in epochs and organizes them into
%
% USAGE: [F, S, SR] = collect_features(R)
%
% Input R (structure with field obtained from slice_io.m
%   R.spikes = raw traces of the spikes
%   R.S = cell array of spike times etc for spikes in each epoch
%   R..mp = various membrane parameters
% Output:
%   F = cell array of features
%   S = cell array of spikes
%   SR = cell array of sampling rates - (needed to reconstruct the time
%       series

if nargin < 2; ap = setAnalysisParams(); end;
%you're calling this thing too much, it makes the code rlly slow...

%somehow MATLAB was not assigning things...
F = {};
S = {};
SR = {};

% Interate over each epoch
for e=1:numel(R.spikes)        
    % Interate over each spike
    count = 0;
    if min(size(R.spikes{e}))
        for s = 1:min(size(R.spikes{e}))
            count = count + 1;
            F{e}{count} = spike_features(R.spikes{e}(s,:),  R(1).S.sampRate, ap);
            S{e}{count} = R.spikes{e}(s,:);
            SR{e}{count} = R.S.sampRate; % Keep the sampling rate so can plot each spike properly
        end       
    end
end
    
