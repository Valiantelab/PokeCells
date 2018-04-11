function [wf] = waveform_create(ap, T, npulses)

pulse_samples = find(T >= ap.io.pulsestart & T < (ap.io.pulsestart + ap.io.pulsedur));

pulse = zeros(1,length(T));
pulse(pulse_samples) = 1;

ps = (0:(npulses-1))*ap.io.pampstep+ap.io.pampstart;
wf = (repmat(pulse, npulses, 1).*repmat(ps,length(T),1)')';




