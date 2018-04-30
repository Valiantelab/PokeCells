function [R] = sliceIO(dataAnalysisParams, doSave)

% Function to analyze some rudimentary input output relationship for
% neurons

% USAGE [S, spikes, mp] = sliceIO(dap, ap, dosave)
%
% INPUT:
%   dataAnalysisParams - data analysis params, obtained from the excel_read function
%   doSave - save stuff

% OUTPUT:
% R -   Results structure (see below)

%  Taufik A Valiante (2018)

global DATA_DIR;


if nargin < 3; doSave = false; end

dataChannel = dataAnalysisParams.io.ch_to_analyze;
fn = dataAnalysisParams.fname(1:end-1); % FileName


%if nargin < 2; 
%analysisParams = setAnalysisParams; end

% Get the data
% -------------------------------------------------------------
dpath = getAbsoluteFilePath(DATA_DIR, dataAnalysisParams.Dir);
fpath = fullfile(dpath, [fn '.abf']);

if ~exist(fpath, 'file')
    disp(fpath);
    error('File does not exist!');
end

[data,sampInterval,hdr]=abfload(fpath); % Load the data
numEpisodes = hdr.lActualEpisodes;
if numEpisodes >= 40
    disp('Some other protocol')
    return
end
if hdr.nOperationMode ~= 5
    % Means that data were acquired in some other mode than that required
    % for I/O
    error('Data not acquired in: waveform fixed-length mode');
end

%------------Compute a few constants -----------------------%

[nPoints, nChannels, nEpochs] = size(data);
hdr.lNumSamplesPerEpisode = nPoints;
hdr.nADCNumChannels = nChannels;
hdr.lActualEpisodes = nEpochs; %episodes
dataAnalysisParams.io.maxT = dataAnalysisParams.io.pulsestart + dataAnalysisParams.io.pulsedur;

sampRate = 1/(sampInterval*1e-6); % convert sample interval to frequency
tvec = (sampInterval*1e-3)*(0:(hdr.lNumSamplesPerEpisode-1)); % Time vector

% Check for a voltage offset that at times can be quite significant
if ~isempty(find_text(dataAnalysisParams.cond.names, 'offset'))
    vOffset = dataAnalysisParams.cond.times(1); % Set the offset to the first value the offset field
else
    vOffset = 0;
end

% Keep track of the figures for printing
fcount = 0;
figs = [];

% Plot the raw data -----------------------------------------------------
[fcount, figs] = figure_set(fcount, figs, [upper(fn) '-Raw Data']);
for i=1:nChannels
    subplot(nChannels, 1, i)
    D = squeeze(data(:,i,:));
    plot(tvec,D(:,[1 end]));
    xlabel('Time(ms)');
    ylabel(hdr.recChUnits{i});
    axis([0 tvec(end) [min(min(D))-10 max(max(D))+10]]);
    axes_text_style();
end

% Create the waveform if none exists ---------------------------------

if nChannels == 1
    %%%------------ Let's ignore it for now ------------%%%
    %disp('Only one channel!');
    %disp(fpath);
    %return
    %number of channels
    I = waveform_create(dataAnalysisParams, tvec, size(data,3));
    V = squeeze(data(:,1,:));
    vchan = 1;
    wfend = (dataAnalysisParams.io.pulsestart + 1.01*dataAnalysisParams.io.pulsedur)*sampRate/1000; % pulse timings in millis
else
    % Get voltage channel and current channel in spreadsheet otherwise just
    % default
    disp('Two channels loaded')
    disp(fpath);
    
    %%--Remove this?--%%
    vchan = find_text(dataAnalysisParams.chlabels, 'V');
    ichan = find_text(dataAnalysisParams.chlabels, 'I');
    if isempty(vchan) || isempty(ichan)
        vchan = 1;
        ichan = 2;
    end
    
    V = squeeze(data(:,vchan,:));
    I = squeeze(data(:,ichan,:));
        
    % Timing of current waveform 
    [~, wfstart] = min(diff(I(:,1))); % derivative marks onset & offset
    [~, wfend] = max(diff(I(:,1)));
    dataAnalysisParams.io.pulsestart = wfstart/sampRate*1000; % when the pulse starts
    dataAnalysisParams.io.pulsedur = (wfend - wfstart)/sampRate*1000; % when the pulse ends
    dataAnalysisParams.io.maxT = dataAnalysisParams.io.pulsestart + dataAnalysisParams.io.pulsedur; % end of the paradigm
    
    % Amplitude of the waveform
    steps = mean(I(wfstart:wfend,:));
    dataAnalysisParams.io.pampstep = mean(diff(steps));
    dataAnalysisParams.io.pampstart = min(steps);
end

V = V + vOffset; % V contains all the voltage traces from all episodes..

% Plot the waveform -----------------------------------------------------
[fcount, figs] = figure_set(fcount, figs, [upper(fn) ' Waveform']);
plot(tvec,I(:,[1 end]));
ylim([min(ylim)-100 max(ylim) + 100]);
axes_text_style();

%---------------------- Membrane properties and the like -----------------%
% Get the memobrane properties
[fcount, figs] = figure_set(fcount, figs, [upper(fn) '-Ih Calculation']);
[mp] = membrane_properties(V, hdr, dataAnalysisParams, sampRate);
        
[fcount, figs] = figure_set(fcount, figs, [upper(fn) '-Input resistance']);
plot(mp.pulses(1:length(mp.sag_peak)),mp.sag_peak, '.', 'MarkerSize', 20);
hold on;
plot(mp.xfit, mp.yfit, 'r');
hold off;
title(sprintf('Resistance = %4.0f MOhm Rmp = %4.0fmV', mp.b(2)*1e3, mp.b(1)));
axes_text_style();
xlabel('pA');
ylabel('mV');

%----------------- Spike analysis ---------------------------------------%

% Find all the spikes
[FoundSpikes] = findSpikes(V, tvec, dataAnalysisParams);

% Store some parameters
FoundSpikes.sampRate = sampRate;
FoundSpikes.tvec = tvec;
FoundSpikes.hdr = hdr;

% Get plot layout
[r, c] = rc_plot(numEpisodes+1);

% Maybe at some point 2 cells (or more) will be recorded from??
cmap = colormap(lines);

ftext = sprintf('%s-SPIKES Channel # %d',upper(fn), dataChannel);
[fcount, figs] = figure_set(fcount, figs, ftext);

R.numEpisodes = {numEpisodes};
for episode = 1:numEpisodes
    ax(episode) = subplot(r,c,episode);
    plot(tvec,V(:,episode));
    xlabel('Time(ms)', 'FontSize', 6);
    ylabel(hdr.recChUnits{vchan}, 'FontSize', 6);
    axis([0 tvec(end) dataAnalysisParams.io.yaxis]);
    
    % Plot the peaks
    if ~isempty(FoundSpikes.pks{episode})   
        hold on;
        plot(tvec(FoundSpikes.locs{episode}), FoundSpikes.pks{episode}, '.r', 'LineStyle', 'none');
        hold off;
    end
    lv(episode) = locVar(tvec, FoundSpikes.locs{episode}, dataAnalysisParams);
    
    if lv(episode) >= 0
        title(sprintf('Lv = %6.4f', lv(episode)));
    else
        title(sprintf('Episode #%d', episode));
    end
    set(gca, 'FontSize', 6);
end
linkaxes(ax,'xy');

% Store the Lvs of Shinomoto
FoundSpikes.lv = lv;

ccount = 0;
l_text = {};
% Plot the ISI as a function of spike number
subplot(r,c,[numEpisodes+1 r*c]);
for episode = 1:numEpisodes
    hold on;
    if ~(FoundSpikes.num{episode}< dataAnalysisParams.io.minspikes) && (min(tvec(FoundSpikes.locs{episode}))< dataAnalysisParams.io.maxT)
    %if ~(length(FoundSpikes.locs{episode}) < dataAnalysisParams.io.minspikes) && (min(tvec(FoundSpikes.locs{episode}))< dataAnalysisParams.io.maxT)
        ccount = ccount + 1;
        ISI = diff(tvec(FoundSpikes.locs{episode}));
        plot(ISI, 'Color',cmap(ccount,:));
        l_text{ccount} = sprintf('%d',episode);
        FoundSpikes.isi{episode} = ISI;
    else
        FoundSpikes.isi{episode} = {};
    end 
end
hold off;
xlabel('Spike number');
ylabel('ISI (ms)');
legend(l_text, 'Location', 'EastOutside');
axis([0 max(xlim) 0 max(ylim)]);

set(gca, 'FontSize', 6);

% Keep some paramaters that might have changed
FoundSpikes.ap = dataAnalysisParams;

%----------------- Individual Spike display --------------------------%

% Display all the spikes
spikeWindow = dataAnalysisParams.io.spikewindowtime*1e-3*sampRate;
spiketVec = (-spikeWindow:spikeWindow)/sampRate*1e3;

ftext = sprintf('%s - ISI Channel # %d',upper(fn), vchan);
[~, figs] = figure_set(fcount, figs, ftext);

SpikeWf = {};
for episode = 1:numEpisodes
    ax(episode) = subplot(r,c,episode);
    set(gca, 'FontSize', 6);
    %if ~isempty(length(FoundSpikes.locs{j}))
    if (FoundSpikes.num{episode}) %skip if 0 spikes found
        hold on;
        sp_count = 0;
        for spikeIndex = 1:FoundSpikes.num{episode}
            if tvec(FoundSpikes.locs{episode}(spikeIndex)) < dataAnalysisParams.io.maxT && FoundSpikes.locs{episode}(spikeIndex) > spikeWindow
                sp_count = sp_count +1;
                startIndex = FoundSpikes.locs{episode}(spikeIndex)-spikeWindow;
                endIndex = FoundSpikes.locs{episode}(spikeIndex)+spikeWindow;
                SpikeWf{episode} = data(startIndex:endIndex, dataChannel, episode)';
                
                % Display the first x spikes 
                if (sp_count <= dataAnalysisParams.io.firstspikestodisp)
                    plot(spiketVec, SpikeWf{episode}, 'Color', cmap(sp_count,:));
                else
                    plot(spiketVec,SpikeWf{episode}, 'Color', [.7 .7 .7]);
                end
            end
        end

        hold off;
        %SpikeWf{episode} = sp;

        title(sprintf('Episode #%d', episode));
        axis([spiketVec(1) spiketVec(end) dataAnalysisParams.io.yaxis]);
    end
end
linkaxes(ax,'xy');

ccount = 0;
l_text = {};
% Plot the spike amplitudes
subplot(r,c,[numEpisodes+1 r*c]);
max_peaks = -1;
for episode = 1:numEpisodes
    hold on;
    if ~isempty(FoundSpikes.pks{episode})
        ccount = ccount + 1;
        plot(FoundSpikes.pks{episode}, 'Color',cmap(ccount,:));
        l_text{ccount} = sprintf('%d',episode);
        if length(FoundSpikes.pks{episode}) > max_peaks
            max_peaks = length(FoundSpikes.pks{episode});
        end
    end 
end
hold off;
xlabel('Spike number');
ylabel('Peak (mV)');
legend(l_text, 'Location', 'EastOutside');
if max_peaks > 0
    axis([0 max_peaks ylim]);
end
set(gca, 'FontSize', 6);

%-------------------------- SUMMARY NUMBERS -----------------------------%
if  prod(double(cellfun('isempty', FoundSpikes.pks))) == 1
    disp('NO spikes detected for this cell...no summary stats')
    return;
end

R.spikes = SpikeWf;
R.S = FoundSpikes;
R.mp = mp;

[F, ~, ~] = collect_features(R);

feat{1} = F;
[all_features] = collapse_feature(feat, dataAnalysisParams.io.features, dataAnalysisParams.io.normalize);

%----------------------- Export the plots -------------------------------%
global FIGURE_DIR;

if ~isempty(FIGURE_DIR)
    figure_batch_save(figs, FIGUR_DIR, doSave);
else
    disp('No figure directory specified to export to.');
end



%Rebound Spiking --------------------------
close all;

ReboundSpikes = getReboundSpikes(data(wfend:end, dataChannel, :), dataAnalysisParams); %on all episodes
R.ReboundSpikes = ReboundSpikes;

%AHP calculation --------------------------
%computing which currents are hyperpol and which are depol
pampend = dataAnalysisParams.io.pampstartIdeal + dataAnalysisParams.io.pampstepIdeal * (numEpisodes - 1); 
dataAnalysisParams.io.idealI = dataAnalysisParams.io.pampstartIdeal:dataAnalysisParams.io.pampstepIdeal:pampend; %store the set I values

idealI = dataAnalysisParams.io.idealI;
posI = idealI(idealI>0);
negI = idealI(idealI<0);

%restingPotential = mp.resting; %this is average over all trials... let's
%compute it individually instead
plottingWindow = 500;

ahpFigHandle = figure;
title('AHP');

reboundFigHandle = figure;
title('Rebound');

adpFigHandle = figure;
title('ADP');

R.ahp = {};
for episode = 1:numEpisodes
    % use wf end to find AHP
    %plot(plottingWindow, mp.resting, 'kx', 'markerSize', 12);
    %disp('Episode');
    %disp(episode)
    
    if (ReboundSpikes.num{episode})
        %disp('------------ Rebound spikes, no AHP ------------');
        figure(reboundFigHandle);
        hold on;
        %plot(data(wfend-plottingWindow:end, dataChannel, episode));
        plot(data(:, dataChannel, episode));
        hold off;
        R.ahp{episode}.auc = {};
        R.ahp{episode}.min = {};
    elseif (idealI(episode) > 0) % no rebound spiking, and depolarizing current - do AHP stuff
        %disp('------------ No rebound spikes, computing AHP ------------');
        figure(ahpFigHandle);
        hold on;
        %plot(data(wfend-plottingWindow:end, dataChannel, episode));
        plot(data(:, dataChannel, episode))
        hold off;
        R.ahp{episode} = findAhp(V(:, episode), tvec, dataAnalysisParams, wfend, sampRate, episode);
    elseif (idealI(episode) < 0)
        %disp('------------Hyperpolarizing Current------------');
        %----Can compute ADP----%
        figure(adpFigHandle);
        hold on;
        %plot(data(wfend-plottingWindow:end, dataChannel, episode));
        plot(data(:, dataChannel, episode))
        hold off;
        R.ahp{episode}.auc = {};
        R.ahp{episode}.min = {};
    else
        R.ahp{episode}.auc = {};
        R.ahp{episode}.min = {};
        % 0 pA
    end
end



%---------------------   OTHER FUNCTIONS --------------------------------%

      
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

