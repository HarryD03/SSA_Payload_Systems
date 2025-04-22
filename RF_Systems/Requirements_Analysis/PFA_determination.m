clear; clc; close all

% Define parameters
pfa_min = 1e-10;       % minimum probability of false alarm
pfa_max = 1e-1;        % maximum probability of false alarm
num_points = 1000;     % number of points in the linspace vector

% Create a linspace vector for the probability of false alarm
pfa_vec = logspace(-10,-1,100);

% Calculate the number of detection events required for a false alarm.
% For a given pfa, the expected number of detections until a false alarm is 1/pfa.
detection_events = 1 ./ pfa_vec;

% Divide the number of detection events by 125 (the minimum number of detections)
% to get the number of observed debris before a false alarm is raised.
observed_debris_before_FA = detection_events / 125;

% Define the debris population
debris_population = 1100000;

% Calculate the number of failures (i.e., false alarms) till the catalogue is built.
% This is the debris population divided by the number of observed debris per false alarm.
failures_till_catalogue = debris_population ./ observed_debris_before_FA;
% Note: This simplifies to: failures_till_catalogue = debris_population * 125 * pfa_vec;

second_visit_failures= failures_till_catalogue./observed_debris_before_FA;

additional_time_hrs = (failures_till_catalogue - second_visit_failures) * 5400/3600;

base_time_hrs = 5208;

% Plot the results using a logarithmic scale for the x-axis
figure;
semilogx(pfa_vec, failures_till_catalogue, 'LineWidth', 2);
xlabel('Probability of False Alarm (pfa)');
ylabel('Number of Failures Till Catalogue Completion');
title('Catalogue Failures vs. Probability of False Alarm');
grid on;
ylim([0 debris_population])
