%% Minimum Range Calculation for 125 Samples Acquisition
% The sensor takes 125 samples with a time between samples of 1e-3 s.
% Therefore, the total dwell time required is:
dt = 1e-3;                   % Sampling interval [s]
num_samples = 125;           % Number of samples
T_dwell = dt * num_samples;  % Total dwell time [s] = 0.125 s

% Sensor field-of-view (FOV) parameters:
half_cone_angle_deg = 11;              % Half cone angle [deg]
half_cone_angle_rad = deg2rad(half_cone_angle_deg);  % Convert to radians

% Define debris speed range (in km/s)
v_min = 0.1;        % Minimum debris speed [km/s]
v_max = 14;         % Maximum debris speed [km/s]
num_v = 100;        % Number of speed points to evaluate
v_vec = linspace(v_min, v_max, num_v);  % Debris speed vector [km/s]

% Compute minimum range required for the debris to be in view for T_dwell seconds.
% The debris traverses the full FOV diameter (2*R*tan(half_angle)) in time T_dwell.
% Therefore, R_min is given by:
%   T_dwell = (2*R_min*tan(half_cone_angle)) / v
%   R_min = (T_dwell * v) / (2*tan(half_cone_angle))
R_min = (T_dwell .* v_vec) ./ (2 * tan(half_cone_angle_rad));


% Plot the results
figure;
plot(v_vec, R_min, 'b-', 'LineWidth', 2);
xlabel('Debris Speed (km/s)');
ylabel('Minimum Range (km)');
title('Minimum Range for 125 Samples vs. Debris Speed');
legend('Calculated R_{min}', 'Clipped to R_{max} (50 km)', 'Location', 'Best');
grid on;
