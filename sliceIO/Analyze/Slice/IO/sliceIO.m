function [R] = sliceIO(dap, ap, dosave)

% Function to analyze some rudimentary input output relationship for
% neurons

% USAGE [S, spikes, mp] = sliceIO(dap, ap, dosave)
%
% INPUT:
%   dap - data analysis params (dap), obtained from the excel_read function
%   ap  - analysis params (ap), can be specified for each dap (if
%           required), can be helpful when batch processing stuff
%   dosave - save stuff

% OUTPUT:
% R -   Results structure (see below)

%  Taufik A Valiante (2018)

global DATA_DIR;

if nargin < 3; dosave = false; end
if nargin < 2; ap = sl_sync_params; end
fn = dap.fname(1:end-1); % FileName


% Get the data
% -------------------------------------------------------------
dpath = getAbsoluteFilePath(DATA_DIR, dap.Dir);
fpath = fullfile(dpath, [fn '.abf']);

if ~exist(fpath, 'file')
    disp(fpath);
    error('File does not exits');
end

[d,si,hdr]=abfload(fpath); % Load the data

if hdr.nOperationMode ~= 5
    % Means that data were acquired in some other mode than that required
    % for I/O
    error('Data not acquired in: waveform fixed-length mode');
end

%------------Compute a few constants -----------------------%

[npoints, nchan, nepochs] = size(d);
hdr.lNumSamplesPerEpisode = npoints;
hdr.nADCNumChannels = nchan;
hdr.lActualEpisodes = nepochs;
ap.io.maxT = ap.io.pulsestart + ap.io.pulsedur;

sr = 1/(si*1e-6); % convert sample interval to frequency
T = (si*1e-3)*(0:(hdr.lNumSamplesPerEpisode-1)); % Time vector

% Check for a voltage offset that at times can be quite significant
if ~isempty(find_text(dap.cond.names, 'offset'))
    vOffset = dap.cond.times(1); % Set the offset to the first value the offset field
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
    D = squeeze(d(:,i,:));
    plot(T,D(:,[1 end]));
    xlabel('Time(ms)');
    ylabel(hdr.recChUnits{i});
    axis([0 T(end) [min(min(D))-10 max(max(D))+10]]);
    axes_text_style();
end



% Create the waveform if none exists ---------------------------------

if hdr.nADCNumChannels == 1
    wf = waveform_create(ap, T, size(d,3));
    V = squeeze(d(:,1,:));
else
    % Get voltage channel and current channel in spreadsheet otherwise just
    % default
    
    vchan = find_text(dap.chlabels, 'V');
    ichan = find_text(dap.chlabels, 'I');
    if isempty(vchan) || isempty(ichan)
        vchan = 1;
        ichan = 2;
    end
    
    V = squeeze(d(:,vchan,:));
    I = squeeze(d(:,ichan,:));
    wf = I; % Current Waveform
        
    % Timing of waveform 
    [~, wfstart] = min(diff(wf(:,1))); % derivative marks onset & offset
    [~, wfend] = max(diff(wf(:,1)));
    ap.io.pulsestart = wfstart/sr*1000; % when the pulse starts
    ap.io.pulsedur = (wfend - wfstart)/sr*1000; % when the pulse ends
    ap.io.maxT = ap.io.pulsestart + ap.io.pulsedur; % end of the paradigm
    
    % Amplitude of the waveform
    steps = mean(wf(wfstart:wfend,:));
    ap.io.pampstep = mean(diff(steps));
    ap.io.pampstart = min(steps);
end

V = V + vOffset;

% Plot the waveform -----------------------------------------------------
[fcount, figs] = figure_set(fcount, figs, [upper(fn) ' Waveform']);
plot(T,wf(:,[1 end]));
ylim([min(ylim)-100 max(ylim) + 100]);
axes_text_style();

%---------------------- Membrane properties and the like -----------------%
% Get the memobrane properties
[fcount, figs] = figure_set(fcount, figs, [upper(fn) '-Ih Calculation']);
[mp] = membrane_properties(V, hdr, ap, sr);
        
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
[S] = find_spikes(V, T, ap);

% Store some parameters
S.sr = sr;
S.T = T;
S.hdr = hdr;

% Get plot layout
[r, c] = rc_plot(hdr.lActualEpisodes+1);


% Maybe at some point 2 cells (or more) will be recorded from??
cmap = colormap(lines);

ftext = sprintf('%s-SPIKES Channel # %d',upper(fn), ap.io.ch_to_analyze);
[fcount, figs] = figure_set(fcount, figs, ftext);

for j = 1:hdr.lActualEpisodes
    ax(j) = subplot(r,c,j);
    plot(T,V(:,j));
    xlabel('Time(ms)', 'FontSize', 6);
    ylabel(hdr.recChUnits{vchan}, 'FontSize', 6);
    axis([0 T(end) ap.io.yaxis]);
    
    % Plot the peaks
    if ~isempty(S.pks{j})   
        hold on;
        plot(T(S.locs{j}), S.pks{j}, '.r', 'LineStyle', 'none');
        hold off;
    end
    lv(j) = locVar(T, S.locs{j}, ap);
    
    if lv(j) >= 0
        title(sprintf('Lv = %6.4f', lv(j)));
    else
        title(sprintf('Episode #%d', j));
    end
    set(gca, 'FontSize', 6);
end
linkaxes(ax,'xy');

% Store the Lvs of Shinomoto
S.lv = lv;

ccount = 0;
l_text = {};
% Plot the ISI as a function of spike number
subplot(r,c,[hdr.lActualEpisodes+1 r*c]);
for j = 1:hdr.lActualEpisodes
    hold on;

%     % Subtract off the resting membrane potential
%     if ~isempty(S.pks{j})
%         S.pks{j} = S.pks{j}- mp.resting;
%     end

    if ~(length(S.locs{j}) < ap.io.minspikes) && (min(T(S.locs{j}))< ap.io.maxT)
        ccount = ccount + 1;
        ISI = diff(T(S.locs{j}));
        plot(ISI, 'Color',cmap(ccount,:));
        l_text{ccount} = sprintf('%d',j);
        S.isi{j} = ISI;
    else
        S.isi{j} = {};
    end 
end
hold off;
xlabel('Spike number');
ylabel('ISI (ms)');
legend(l_text, 'Location', 'EastOutside');
axis([0 max(xlim) 0 max(ylim)]);

set(gca, 'FontSize', 6);

% Keep some paramaters that might have changed
S.ap = ap;

%----------------- Individual Spike display --------------------------%

% Display all the spikes
w = ap.io.spikewindowtime*1e-3*sr;
t = (-w:w)/sr*1e3;

ftext = sprintf('%s - ISI Channel # %d',upper(fn), vchan);
[~, figs] = figure_set(fcount, figs, ftext);

for j = 1:hdr.lActualEpisodes
    ax(j) = subplot(r,c,j);
    set(gca, 'FontSize', 6);
    if ~isempty(length(S.locs{j}))
        hold on;
        sp = [];
        sp_count = 0;
        for s = 1:length(S.locs{j})
            if T(S.locs{j}(s)) < ap.io.maxT && S.locs{j}(s) > w
                sp_count = sp_count +1;
                sp(sp_count,:) = d((S.locs{j}(s)-w):(S.locs{j}(s)+w),ap.io.ch_to_analyze,j);
                % Display the first x spikes 
                if (sp_count <= ap.io.firstspikestodisp)
                    plot(t,sp(sp_count,:), 'Color', cmap(sp_count,:));
                else
                    plot(t,sp(sp_count,:), 'Color', [.7 .7 .7]);
                end
            end
        end
%         if sp_count>1
%             plot(t,mean(sp),'Color', [0 0 0]);
%         end
        hold off;
        spikes{j} = sp;
%         ylabel(hdr.recChUnits{ap.io.ch_to_analyze});
%         xlabel('Time (ms)');
        title(sprintf('Episode #%d', j));
        axis([t(1) t(end) ap.io.yaxis]);
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
    if ~isempty(S.pks{j})
        ccount = ccount + 1;
        plot(S.pks{j}, 'Color',cmap(ccount,:));
        l_text{ccount} = sprintf('%d',j);
        if length(S.pks{j}) > max_peaks
            max_peaks = length(S.pks{j});
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
if  prod(double(cellfun('isempty', S.pks))) == 1
    display('NO spikes detected for this cell...no summary stats')
    return;
end

R.spikes = spikes;
R.S = S;
R.mp = mp;

[F, ~, ~] = collect_features(R);

feat{1} = F;
[all_features] = collapse_feature(feat, ap.io.features, ap.io.normalize);

% % Make sure there is something in here to plot lest and error should occur
% if min(size(all_features))
%     ftext = sprintf('%s - FEATURES - Channel # %d',upper(fn), ap.io.ch_to_analyze);
%     [fcount, figs] = figure_set(fcount, figs, ftext);
%     plot_features(all_features,ap.io.firstspikestodisp, ap.io.features, ap.io.normalize);
% end

% Do some stats
%features_stats(all_features,ap);

%----------------------- Export the plots -------------------------------%
global FIGURE_DIR;

if ~isempty(FIGURE_DIR)
    figure_batch_save(figs, FIGUR_DIR, dosave);
else
    disp('No figure directory specified to export to.');
end
%---------------------   OTHER FUNCTIONS --------------------------------%

      
function [S] = find_spikes(d, T, ap)

[~, nepochs] = size(d);
S = [];

warning('off');

for j=1:nepochs
    ts = squeeze(d(:,j));
    [S.pks{j},S.locs{j}] = findpeaks(ts, 'MINPEAKHEIGHT', ap.io.minpeakheight);
    ind = find(T(S.locs{j}) > ap.io.pulsestart  & T(S.locs{j}) < (ap.io.pulsestart+ap.io.pulsedur));
    if ~isempty(ind)
        S.pks{j} = S.pks{j}(ind);
        S.locs{j} = S.locs{j}(ind);
    else
        S.pks{j} = [];
        S.locs{j} = [];
    end
end

warning('on');

