function [analysisParams] = setAnalysisParams()
% ---- Directory information initialized
analysisParams.Dir = '';
analysisParams.ExportDir = [];
analysisParams.fname = '';
analysisParams.comment = {};

% ---- Condition infomation
% For fname it is explicitly assumed that the files for the different
% conditons exist in the same directory

analysisParams.cond.fname = {};     % Filenames for the various conditons
analysisParams.cond.names = {};     % Names of the various conditions
analysisParams.cond.times = [];     % Times of the various conditons
analysisParams.chlabels = {};


% ---- Channel information
analysisParams.ch = 1;              % Channel to analyze

% Important that this is an integer divider of the original sampling rate
analysisParams.resample = false;
analysisParams.srate = 2500;        % Sampling rate to decimate data to - (why?)

% Factor by which to oversample the data - thus is ap.srate = 5000, the low
% pass filter will be set to 1000;
analysisParams.over_sample = 5;     


% Data source
analysisParams.load_mat = 1;

% Preprocessing
analysisParams.prep_all = 1;                % Preprocess all the files in the spreadsheet, otherwise just those with GroupInc = 'yes';

% ---- Filtering
analysisParams.notches = [60];        % Principal frequncies for notch filtering
analysisParams.notch = 0;                   % Do notch (1) or not (0)
analysisParams.nharm = 5;                   % Number of harmonics to use for notch filtering
analysisParams.bstop = 2;                   % Notch width in Hz

analysisParams.filter.ftype = 'firnotch';
analysisParams.iirnotch.order = 10;
analysisParams.iirnotch.ftype = 'butter';
analysisParams.firnotch.order = 1000;

analysisParams.notch_ch_list = [1 0];

% ---- Special fields
analysisParams.sf.names = {'Comment', 'ExportDir', 'Tags', 'GroupInc'};

%ap.sf.names = {};
analysisParams.sf.vals = {};

%% PLOTTING
analysisParams.pl.show_axes = 1;
analysisParams.pl.colorbar = 1;
analysisParams.pl.ranges = 1;
analysisParams.pl.zeroline = 1;
analysisParams.yaxis = [0 0.6];
analysisParams.plotsync = 0;
analysisParams.pl.axis = 'linear';

analysisParams.pl.textprop = {'FontName', 'FontSize', 'TickDir'};
analysisParams.pl.textpropval = {'Times', 7, 'out'};

analysisParams.pl.axprop = {'TickDir'};
analysisParams.pl.axpropval = {'out'};

% High pass filtering for display
analysisParams.disp_filter.on = 1;  % Turn on the filer
analysisParams.disp_filter.Fc = 4;  % Hz
analysisParams.disp_filter.order = 1000;  % Order of the filter

% for displaying raw data traces
analysisParams.raw.yaxis_range = [];
analysisParams.raw.xaxis_link = 0;

%-----------I/O--------------------%
analysisParams.io.minpeakheight = 0;
analysisParams.io.yaxis = [-210 120];
analysisParams.io.ch_to_analyze = 1;
analysisParams.io.minspikes = 3;
analysisParams.io.maxT = 765; % in ms
analysisParams.io.spikewindowtime = 5;
analysisParams.io.firstspikestodisp = 3;
analysisParams.io.spike_threshold = 5000; % factor increase in rate of rise that signals start of spike

analysisParams.io.sstate_dur = 200;
analysisParams.io.features = {'maxRise', 'maxFall', 'ahp', 'w', 'amp', 'peak'};
analysisParams.io.alpha = 0.05;
analysisParams.io.normalize = 0; % Normalization of the features: 1- within epoch, 2- within celll
analysisParams.io.min_isi = 3;

analysisParams.io.nudge = 10; % Nudge the start of the pulse due to measurement jitter
analysisParams.io.fit_dur = 50;

% LV
analysisParams.io.lv_minisi = 2;
analysisParams.io.lv_min_n_toplot = 5;

analysisParams.io.Ih_delta = 20;
analysisParams.io.Ih_threshold_peak = 1;
analysisParams.io_Ih_threshold_n = 3;
analysisParams.io_Ih_minpointstofit = 3;
analysisParams.io.Ih_pulsestop = 1;

%Homeira
% ap.io.pulsestart = 160;
% ap.io.pulsedur = 600;
% ap.io.pampstart = -400;
% ap.io.pampstep = 50;

% Plotting
analysisParams.io.isi_axis = [0 15 0 100];
analysisParams.io.lv_axis = [0 11 0 0.6];
analysisParams.io.Ih_yaxis = [-130 -50 0 10];

