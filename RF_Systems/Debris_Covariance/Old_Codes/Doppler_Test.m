% =================================================
% Pulse-Doppler Radar: Single Target Example
% =================================================
clear; clc; close all;

% -----------------------------
% 1) Radar & Waveform Parameters
% -----------------------------
c          = 3e8;         % Speed of light (m/s)
fc         = 10e9;        % Carrier frequency (Hz)
lambda     = c / fc;      % Wavelength (m)

PRF        = 30e3;         % Pulse repetition frequency (Hz)
Tp         = 1/PRF;       % Pulse repetition interval (s)
N_pulses   = 32;          % Number of pulses (slow-time samples)

pulseWidth = 2e-6;        % Pulse width (s)
fs         = 10e6;        % Fast-time sampling rate (Hz)
N_fast     = round(pulseWidth * fs);  % # of samples in one pulse

% -----------------------------
% 2) Single Target Definition
% -----------------------------
% Let's specify the target by range, velocity, and amplitude.
r_tgt    = 300;    % target range (m)
v_tgt    = 50;     % target radial velocity (m/s) 
                   %  (positive => receding, negative => approaching)
amp_tgt  = 1.0;    % amplitude (relative scale)

% Round-trip time delay for initial range
%  tau = 2*r / c  (out and back)
tau0     = 2 * r_tgt / c;

% Monostatic Doppler shift for the target 
fD_tgt   = 2 * v_tgt / lambda;  % (Hz)

% We'll store the echoes from each pulse in a 2D array: [fast-time, pulses]
rxData   = zeros(N_fast, N_pulses);

% Define a simple rectangular transmit pulse in fast-time
txPulse  = ones(1, N_fast);   % amplitude 1, width = pulseWidth

% -----------------------------
% 3) Simulate Received Pulses
% -----------------------------
for n = 1:N_pulses
    % Time of the start of the n-th pulse (slow-time)
    t_pulseStart = (n-1)*Tp;
    
    % Range at the n-th pulse (assuming uniform motion)
    % r(n) = r_tgt + v_tgt * (n-1)*Tp
    r_n   = r_tgt + v_tgt*(n-1)*Tp;
    tau_n = 2*r_n / c;  % updated round-trip delay at pulse n
    
    % Time vector for the fast-time samples within this pulse
    t_fast = (0 : N_fast-1)/fs;
    
    % Delayed time relative to echo arrival
    t_delayed = t_fast - tau_n;
    
    % Construct the target echo for this pulse
    echo = zeros(size(txPulse));
    
    % Indices where the echo actually appears within the pulse window
    idxInPulse = (t_delayed >= 0) & (t_delayed < pulseWidth);
    if any(idxInPulse)
        % Incorporate Doppler (baseband approach)
        echo(idxInPulse) = amp_tgt * cos(2*pi * fD_tgt * t_delayed(idxInPulse));
    end
    
    % Store in our data matrix
    rxData(:, n) = echo.';
end

% -----------------------------
% Optional: Add Noise
% -----------------------------
SNR_dB = 30;  % desired SNR in dB
signalPower = mean(abs(rxData(:)).^2);
noisePower  = signalPower / (10^(SNR_dB/10));
rxData = rxData + sqrt(noisePower/2)*(randn(size(rxData)) + 1j*randn(size(rxData)));

% -----------------------------
% 4) Range-Doppler Processing
% -----------------------------
% (a) FFT across fast-time for each pulse -> Range 
% (b) FFT across pulses (slow-time) -> Doppler 
% We'll produce a 2D range-Doppler map.

% (a) Range FFT (zero‐pad optional)
rangeFFTsize = 2*N_fast;  % zero‐pad to improve freq resolution
rangeFFT = fft(rxData, rangeFFTsize, 1);
rangeFFT = rangeFFT(1:N_fast, :);  % keep only positive freq half

% (b) Doppler FFT
rdMap = fftshift( fft(rangeFFT, N_pulses, 2), 2 ); 
% rdMap is now [N_fast x N_pulses], with frequency shift around center

% -----------------------------
% 5) Axes for Range & Doppler
% -----------------------------
% Range resolution ~ c / (2 * fs) for simple rectangular pulses
rangeAxis_m = (0 : N_fast-1) * (c / (2 * fs));

% Doppler bin spacing = PRF / N_pulses
% Convert Doppler frequency to velocity: v = (lambda/2)*fD
dopplerFreqAxis = (-N_pulses/2 : N_pulses/2-1) * (PRF / N_pulses);
velAxis_mps     = (lambda/2) * dopplerFreqAxis;

% -----------------------------
% 6) Plot the Range-Doppler Map
% -----------------------------
figure;
imagesc(velAxis_mps, rangeAxis_m, 20*log10(abs(rdMap)));
axis xy; colorbar;
xlabel('Radial Velocity (m/s)');
ylabel('Range (m)');
title('Range-Doppler Map (dB)');
