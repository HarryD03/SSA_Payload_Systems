% Define radar and platform parameters
v = 100;           % Platform speed in m/s
lambda = 0.03;     % Radar wavelength in meters

% Define range of incidence angles (in degrees) and candidate PRF values (in Hz)
theta_deg = linspace(20, 60, 100);      % Incidence angles from 20° to 60°
prf_values = linspace(1000, 15000, 200);  % PRF values from 1 kHz to 15 kHz

% Create a grid of PRF and incidence angles
[PRF, THETA] = meshgrid(prf_values, theta_deg);

% Compute the expected Doppler frequency at each incidence angle
fd = 2 * v * sind(THETA) / lambda;

% Compute the aliasing order: 
% For each (theta, PRF), aliasing order is the number of times the Doppler frequency 
% exceeds half the PRF (i.e. the number of ambiguities).
alias_order = floor(abs(fd) ./ (PRF/2));

% Plot the aliasing order as a function of PRF and incidence angle
figure;
imagesc(prf_values, theta_deg, alias_order);
set(gca, 'YDir', 'normal');  % Corrects the vertical axis so lower angles appear at the bottom
xlabel('PRF (Hz)');
ylabel('Incidence Angle (°)');
title('Zebra Plot: Aliasing Order vs. Incidence Angle and PRF');
colorbar;
colormap(gray);

% Optional: Overlay a contour for the zero-aliasing region
hold on;
contour(prf_values, theta_deg, alias_order, [0 0], 'LineColor', 'r', 'LineWidth', 2);
hold off;
