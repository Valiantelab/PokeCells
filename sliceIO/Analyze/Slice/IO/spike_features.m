function [F] = spike_features(V, sr, ap)
%Function to collect compute various spike features from a raw data trace
% Input:
%   V = voltage trace
%   sr = sampling rate
% Output:
% F - feature structure:
%   maxRise -   maximal rate of rise
%   maxFall -   maximal repolarization rate
%   ahp -       afterhyper-polarization
%   w -         spike witdth at half maximal amplitude
%   amp -       spike amplitude (baseline to peak)
%   peak -      maximal voltage attained

t = (0:(length(V)-1))/sr*1e3; % convert to ms
t = t-t(end)/2;
p_sample = (length(t)-1)/2;

dV = diff(V)*sr; % not all data sampled at the same rate
loc = (length(V)-1)/2;
F.maxRise = max(dV(1:loc));
F.maxFall = max(abs(dV((loc+1):end)));

loc_maxRise = find(dV == F.maxRise);
baseline = mean(V(find(t<-2)));
F.amp = max(V)-baseline;
F.sp = V;
F.peak = max(V);

halfV = find(abs((V(find(t>-1 & t< 1))-baseline)) > F.amp/2);
if isempty(halfV)
    F.w = NaN;
else
    F.w = (halfV(end) - halfV(1))/sr*1e3;
end


tahp = find(t<0 & t> -1);
% [~, locs] = findpeaks(diff(dV(tahp)),'MINPEAKHEIGHT', ap.io.spike_threshold);
% if isempty(locs)
%     locs(1) = 1;
% end

%sp_start = tahp(1)+locs(1)-2;
sp_start =  p_sample-numel(halfV)*2;
if sp_start < 0
    sp_start =1;
end
F.ahp = V(sp_start)-min(V(find(t>0)));
if F.ahp < 0
    F.ahp = 0;
end

% figure(11);clf;plot(t, V);
% hold on;
% plot(t(sp_start), V(sp_start), '.r', 'MarkerSize', 15);

