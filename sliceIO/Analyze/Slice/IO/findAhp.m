function [ahp] = findAhp(V, tvec, analysisParams, wfend, sampRate, episode)
%   findAHP 
%   mode - 'AUC' will return AHP as area under the curve in units msec * mV
%          'MIN' will return AHP

%two methods for finding AHP
ahp.auc = {};
ahp.min = {};
numPoints = size(V, 1);

baselineRange = 1:analysisParams.io.pulsestart*1e-3*sampRate; %differed
%from post current injection bline

baselineVec = ones(1, numPoints);

% using starting for resting potential
restingPotential = mean(V(baselineRange));
stdRestingPotential = std(V(baselineRange));

% resting potential from end
%restingPotential = mean(V(end-2500:end));
%stdRestingPotential = std(V(end-2500:end));

threshBaseline = [restingPotential - 2*stdRestingPotential, restingPotential + 2*stdRestingPotential];
%lower, upper
%use a small sliding window to account for noise
windowSize = 21; %odd value to center at this point


ahpRange = wfend:length(V); %approximate range for AHP
ahpStart = find(V(ahpRange) <= restingPotential); % check when no longer depolarized, works
ahpStart = ahpStart(1) + wfend;

if ahpStart > wfend + length(baselineRange)*1.5;
    disp('Bad start to AHP')
    return
end
ahpRange = ahpStart:length(V);

meanVec = movmean(V(ahpRange), windowSize);

%need to forward seek until the lower bound is passed
seek = find(meanVec <= threshBaseline(1));
if ~isempty(seek)
    seek = seek(1) + ahpStart;
else
    disp('cannot find the End of AHP');
    return
end
%now search for end of AHP

ahpEnd = find(meanVec((seek-ahpStart):end) >= threshBaseline(1));   % check when no longer hyperpolarized, returned to baseline
if ~isempty(ahpEnd)
    ahpEnd = ahpEnd(1) + seek;
    if ahpEnd <= (seek + windowSize*2) %index is too close together
        ahpEnd = find(meanVec((seek - ahpStart + windowSize*2):end) >= threshBaseline(1));   % check when no longer hyperpolarized, returned to baseline
        ahpEnd = ahpEnd(1) + seek + windowSize*2;
    end
else
    disp('cannot find the End of AHP');
    return
end
%found the range of ahp, ahpStart:ahpEnd
%integrate using trapz

centeredV = V - restingPotential; %want the area wrt baseline
movingAvV = movmean(V(ahpStart:ahpEnd), windowSize); %compute a moving average to find minimum
%can get standard deviation for this too

ahp.auc = {trapz(tvec(ahpStart:ahpEnd), centeredV(ahpStart:ahpEnd))};
ahp.min = {min(movingAvV)};
    
% plotting
currFig = figure();
hold on
plot(tvec, V);
plot(tvec, baselineVec * restingPotential);
plot(tvec(ahpStart), V(ahpStart), 'b*', 'MarkerSize', 10)
%plot(tvec, baselineVec * restingPotential + 2 * stdRestingPotential);
plot(tvec, baselineVec * threshBaseline(1));
plot(tvec(ahpEnd), V(ahpEnd), 'r*', 'MarkerSize', 10)
hold off
title(strcat('Episode ', num2str(episode)))

%comment

end

