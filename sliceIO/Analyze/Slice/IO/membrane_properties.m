function [mp] = membrane_properties(data, hdr, analysisParams, sampRate)

% compute the resting membrane potential
mp.resting = mean(mean(mean(data(1:analysisParams.io.pulsestart*1e-3*sampRate,:))));
pulses = analysisParams.io.pampstart + analysisParams.io.pampstep*(0:(hdr.lActualEpisodes-1));
ind = find(pulses == max(pulses(find(pulses <= 0))));

% Points for state value to use for IV
sstime = [analysisParams.io.pulsestart+0.8*analysisParams.io.pulsedur analysisParams.io.pulsestart+analysisParams.io.pulsedur];
srange = fix(sstime*1e-3*sampRate);

% Total time of the trace
T = (1/sampRate*1e3)*(0:(hdr.lNumSamplesPerEpisode-1));

% Steady State values
sstate = squeeze(mean(data(srange(1):srange(2),1:ind)));

% Analyze the sag of Ih
btime = find(T> sstime(1) & T < sstime(2));
subplot(1,2,1);

mp.Ih_peak(1) = 0;

tsag = find(T> (analysisParams.io.pulsestart+analysisParams.io.nudge) & T < (analysisParams.io.pulsestart+analysisParams.io.pulsedur));
mp.Ih_tsag = T(tsag);
for i=1:ind-analysisParams.io.Ih_pulsestop
    
    plot(T,data(:,i), 'k');
    hold on;
    plot(T(tsag),data(tsag,i), '-g');
    
    % Try to refine the fit so it does not explode - this whole aspect of
    % initialization\fitting could be much better done
    
    % beta(1) = A1
    % beta(2) = membrane time constant
    % beta(3) = steady state potential
    % beta(4)= Ih amplitude
    % beta(5) = Ih time constant
    
    
    if (i ==1)
        beta0 = [abs(min(data(tsag,i)))-abs(sstate(i)) 15 sstate(i) 20 50];
    else
        beta0(1) = beta0(1)*abs(pulses(i)/pulses(1));
        beta0(3) = sstate(i);
        beta0(5) = beta0(5)*1.5;
    end
    
    [beta, yfit] = fit_tau(T(tsag),data(tsag,i), beta0);
    
    mp.Ih_beta(i,:) = beta;
    mp.Ih_yfit(i,:) = yfit;    
     
    if any(isnan(beta))
        % the fitting went wrong
        sag_peak(i) = sstate(i);
        mp.Ih_peak(i) = 0;
        
        % Try again with standard starting values
        beta0 = [abs(min(data(tsag,i)))-abs(sstate(i)) 5 sstate(i) 5 75];
    else 
        
        % If a single exp fit make sure beta0 is updated properly
        beta0 = beta;
    
        sag_peak(i) = min(yfit);
        p = yfit(end)-min(yfit);

        if p < 0
            mp.Ih_peak(i) = 0;
            sag_peak(i) = sstate(i);
        else
            mp.Ih_peak(i) = p;
        end
    end
end

axis([analysisParams.io.pulsestart-50 analysisParams.io.pulsestart+analysisParams.io.pulsedur+50 min(min(data))-1 -40]);
axes_my_defaults();
if ~isempty(analysisParams.io.Ih_yaxis)
    ylim(analysisParams.io.Ih_yaxis(1:2));
end
ylabel('mV');
xlabel('Time (ms)');
hold off;

mp.sag_peak = sag_peak;
% Do the IO calculation here
np = size(sag_peak,2);
X = [ones(1,np)' pulses(1:np)'];
b = regress(sag_peak',X);
xfit = pulses(1):0.1:pulses(np);
yfit = b(1) + b(2)*xfit;

mp.resistance = b(2)*1e3;
mp.pulses = pulses;
mp.yfit = yfit;
mp.xfit = xfit;
mp.ind = ind;
mp.sstate = sstate;
mp.b = b;

% Now do the IO since Ih and the sag have been computed

if numel(mp.Ih_peak) >= analysisParams.io_Ih_minpointstofit
    subplot(1,2,2)
    plot(sstate(1:numel(mp.Ih_peak)), mp.Ih_peak, '.-b');
    xlabel('Steady state voltage (mV)');
    ylabel('Peak');

    ssvoltage = sstate(1:numel(mp.Ih_peak));
    X = [ones(size(ssvoltage')) ssvoltage'];
    [b,~,~,~,stats] = regress(mp.Ih_peak',X);
    yfit = b(1) + b(2)*(ssvoltage(1):0.1:ssvoltage(end));
    hold on;
    plot((ssvoltage(1):0.1:ssvoltage(end)), yfit, 'r');
    mp.Ih_b = b;
    title(sprintf('Slope = %6.2f, R^2=%6.2f, F=%6.2f, p=%6.2f', b(2), stats(1), stats(2), stats(3)));
    axes_my_defaults();
    mp.Ih_stats = stats;
    if ~isempty(analysisParams.io.Ih_yaxis)
        axis(analysisParams.io.Ih_yaxis)
    end
else
    mp.Ih_stats = 0;
    display('Ih analysis aborted...no Ih present.')
end

function [beta yfit] = fit_tau(t,y, beta0)
ts = t - t(1);

warning off;
try
    beta = nlinfit(ts,y',@mem_tau_ih, beta0);
catch
    try
        % Try with a single exponential function
        beta = nlinfit(ts,y',@mem_tau_exp1, beta0(1:3));
        yfit = mem_tau(beta, ts);
        plot(ts+t(1),yfit, 'm');
        beta = [beta NaN NaN];
    catch
        beta = [NaN NaN NaN NaN NaN];
        yfit = ones(1,length(ts))*beta0(3);
        plot(ts+t(1),yfit, 'm');
        return;
    end

end
warning off;

yfit = mem_tau_ih(beta, ts);
plot(ts+t(1),yfit, 'm');


