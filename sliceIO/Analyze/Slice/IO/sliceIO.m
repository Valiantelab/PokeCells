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
%if nargin < 2; 
%analysisParams = setAnalysisParams; end
fn = dataAnalysisParams.fname(1:end-1); % FileName

% Get the data
% -------------------------------------------------------------
dpath = getAbsoluteFilePath(DATA_DIR, dataAnalysisParams.Dir);
fpath = fullfile(dpath, [fn '.abf']);

if ~exist(fpath, 'file')
    disp(fpath);
    error('File does not exist!');
end

[data,sampInterval,hdr]=abfload(fpath); % Load the data

if hdr.nOperationMode ~= 5
    % Means that data were acquired in some other mode than that required
    % for I/O
    error('Data not acquired in: waveform fixed-length mode');
end

%------------Compute a few constants -----------------------%

[npoints, nchan, nepochs] = size(data);
hdr.lNumSamplesPerEpisode = npoints;
hdr.nADCNumChannels = nchan;
hdr.lActualEpisodes = nepochs;
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
for i=1:hdr.nADCNumChannels
    subplot(hdr.nADCNumChannels, 1, i)
    D = squeeze(data(:,i,:));
    plot(tvec,D(:,[1 end]));
    xlabel('Time(ms)');
    ylabel(hdr.recChUnits{i});
    axis([0 tvec(end) [min(min(D))-10 max(max(D))+10]]);
    axes_text_style();
end

% Create the waveform if none exists ---------------------------------

if hdr.nADCNumChannels == 1
    %number of channels
    I = waveform_create(dataAnalysisParams, tvec, size(data,3));
    V = squeeze(data(:,1,:));
else
    % Get voltage channel and current channel in spreadsheet otherwise just
    % default
    
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

V = V + vOffset;

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
[Spikes] = findSpikes(V, tvec, dataAnalysisParams);

% Store some parameters
Spikes.sampRate = sampRate;
Spikes.tvec = tvec;
Spikes.hdr = hdr;

% Get plot layout
[r, c] = rc_plot(hdr.lActualEpisodes+1);

% Maybe at some point 2 cells (or more) will be recorded from??
cmap = colormap(lines);

ftext = sprintf('%s-SPIKES Channel # %d',upper(fn), dataAnalysisParams.io.ch_to_analyze);
[fcount, figs] = figure_set(fcount, figs, ftext);

for j = 1:hdr.lActualEpisodes
    ax(j) = subplot(r,c,j);
    plot(tvec,V(:,j));
    xlabel('Time(ms)', 'FontSize', 6);
    ylabel(hdr.recChUnits{vchan}, 'FontSize', 6);
    axis([0 tvec(end) dataAnalysisParams.io.yaxis]);
    
    % Plot the peaks
    if ~isempty(Spikes.pks{j})   
        hold on;
        plot(tvec(Spikes.locs{j}), Spikes.pks{j}, '.r', 'LineStyle', 'none');
        hold off;
    end
    lv(j) = locVar(tvec, Spikes.locs{j}, dataAnalysisParams);
    
    if lv(j) >= 0
        title(sprintf('Lv = %6.4f', lv(j)));
    else
        title(sprintf('Episode #%d', j));
    end
    set(gca, 'FontSize', 6);
end
linkaxes(ax,'xy');

% Store the Lvs of Shinomoto
Spikes.lv = lv;

ccount = 0;
l_text = {};
% Plot the ISI as a function of spike number
subplot(r,c,[hdr.lActualEpisodes+1 r*c]);
for j = 1:hdr.lActualEpisodes
    hold on;
    if ~(length(Spikes.locs{j}) < dataAnalysisParams.io.minspikes) && (min(tvec(Spikes.locs{j}))< dataAnalysisParams.io.maxT)
        ccount = ccount + 1;
        ISI = diff(tvec(Spikes.locs{j}));
        plot(ISI, 'Color',cmap(ccount,:));
        l_text{ccount} = sprintf('%d',j);
        Spikes.isi{j} = ISI;
    else
        Spikes.isi{j} = {};
    end 
end
hold off;
xlabel('Spike number');
ylabel('ISI (ms)');
legend(l_text, 'Location', 'EastOutside');
axis([0 max(xlim) 0 max(ylim)]);

set(gca, 'FontSize', 6);

% Keep some paramaters that might have changed
Spikes.ap = dataAnalysisParams;

%----------------- Individual Spike display --------------------------%

% Display all the spikes
w = dataAnalysisParams.io.spikewindowtime*1e-3*sampRate;
t = (-w:w)/sampRate*1e3;

ftext = sprintf('%s - ISI Channel # %d',upper(fn), vchan);
[~, figs] = figure_set(fcount, figs, ftext);

for j = 1:hdr.lActualEpisodes
    ax(j) = subplot(r,c,j);
    set(gca, 'FontSize', 6);
    if ~isempty(length(Spikes.locs{j}))
        hold on;
        sp = [];
        sp_count = 0;
        for s = 1:length(Spikes.locs{j})
            if tvec(Spikes.locs{j}(s)) < dataAnalysisParams.io.maxT && Spikes.locs{j}(s) > w
                sp_count = sp_count +1;
                sp(sp_count,:) = data((Spikes.locs{j}(s)-w):(Spikes.locs{j}(s)+w),dataAnalysisParams.io.ch_to_analyze,j);
                % Display the first x spikes 
                if (sp_count <= dataAnalysisParams.io.firstspikestodisp)
                    plot(t,sp(sp_count,:), 'Color', cmap(sp_count,:));
                else
                    plot(t,sp(sp_count,:), 'Color', [.7 .7 .7]);
                end
            end
        end

        hold off;
        spikes{j} = sp;

        title(sprintf('Episode #%d', j));
        axis([t(1) t(end) dataAnalysisParams.io.yaxis]);
    end
end
linkaxes(ax,'xy');

ccount = 0;
l_text = {};
% Plot the spike amplitudes
subplot(r,c,[hdr.lActualEpisodes+1 r*c]);
max_peaks = -1;
for j = 1:hdr.lActualEpisodes
    hold on;
    if ~isempty(Spikes.pks{j})
        ccount = ccount + 1;
        plot(Spikes.pks{j}, 'Color',cmap(ccount,:));
        l_text{ccount} = sprintf('%d',j);
        if length(Spikes.pks{j}) > max_peaks
            max_peaks = length(Spikes.pks{j});
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
if  prod(double(cellfun('isempty', Spikes.pks))) == 1
    display('NO spikes detected for this cell...no summary stats')
    return;
end

R.spikes = spikes;
R.S = Spikes;
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
%---------------------   OTHER FUNCTIONS --------------------------------%

      
function [Spikes] = findSpikes(data, tvec, analysisParams)

[~, nepochs] = size(data);
Spikes = [];

warning('off');

for j=1:nepochs
    ts = squeeze(data(:,j));
    [Spikes.pks{j},Spikes.locs{j}] = findpeaks(ts, 'MINPEAKHEIGHT', analysisParams.io.minpeakheight);
    ind = find(tvec(Spikes.locs{j}) > analysisParams.io.pulsestart  & tvec(Spikes.locs{j}) < (analysisParams.io.pulsestart+analysisParams.io.pulsedur));
    if ~isempty(ind)
        Spikes.pks{j} = Spikes.pks{j}(ind);
        Spikes.locs{j} = Spikes.locs{j}(ind);
    else
        Spikes.pks{j} = [];
        Spikes.locs{j} = [];
    end
end

warning('on');

