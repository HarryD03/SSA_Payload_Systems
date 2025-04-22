%% Set the SNR threshold by the Neymann-Pearson Criterion 
% Ref: Richards, M. A. Fundamentals of Radar Signal Processing. New York: McGraw-Hill, 2005, pp 298–336.
% Used to Define SNR Threshold Requirement
clc; clear; close all
% Requires Phased Array System Toolbox

% Define parameters
N = 3;                                % Number of pulses for coherent integration
pfa = [1e-3, 1e-4, 1e-5, 1e-6, 1e-8];       % Desired Probability of False Alarm

% Compute Probability of Detection using rocpfa
[Pd, snr_dB] = rocpfa(pfa, 'NumPulses', N);

% Plot the result
figure;
plot(snr_dB, Pd, 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Probability of Detection');
title(sprintf('Coherent Integration with 1 Pulses', N));

% Add a legend to show the different Pfa values
legendStr = arrayfun(@(x) sprintf('Pfa = %.1e', x), pfa, 'UniformOutput', false);
legend(legendStr, 'Location', 'SouthEast');
