%% Main IOD Script: Take Measurements and Propagate Forward for a Revisit Time
clear; close all
% Purpose and Methodology: 
% 1) Input SAR sizing/uncertainty parameters
% 2) Conduct coordinate transformations from spherical LVLH to Cartesian ECI
% 3b) If not Herricks Gibbs; Assume the angle rate can be extracted by finite differences (verified)
% 4) Pass Measurement Uncertainty through EKF to get Final measurement uncertainty/covariance (not verified)
% 5) Estimate the revisit time required based on the Projected Area of the Optical Instrument and covariance propagation (verified)

%% Constants
mu = 3.986e5;         % Earth's gravitational parameter (km^3/s^2)
Re = 6371;            % Earth radius [km]
% SC parameters TBC 
h_SC = 590;           % [km]
e = 0;
i = 90;
Om = 90;
om = 0;
theta = 0;
c = 3e8;
f = 5e9;
lambda = c/f;
SC_kep = [Re + h_SC, e, deg2rad(i), deg2rad(Om), deg2rad(om), deg2rad(theta)]; % [a e i Om om theta]

%% SAR Sizing/Uncertainty Parameters (Inputs)
Nsub = 8;
v_rel = 7.5;                  % Maximum relative velocity of the debris [km/s]
FOV = 10.1103;                % ScanSAR half cone FOV [deg]
FOV_strip = 10.1103 / Nsub;     % FoV per strip [deg]
Range_max = 50.10;            % [km]
Optical_SW = 100;             % [km] Optical swath width (TBD)
prf = 937500;
% Angle Errors at design point
erange = 46.1e-3;             % Range error [km]
eangle = deg2rad(0.1);        % Angle error [rad]
erange_rate = 0.005e-3;
ediff = 0.01;
eaz_rate = sqrt(2*(eangle^2)); % Finite difference error (Taylor expansion)
eel_rate = eaz_rate;

%% Precompute Range, Az, El for Verification
SW_max = 2 * Range_max * tand(FOV_strip);   % Strip swath width [km]
t_total = SW_max / (8 * v_rel);               % Timespan in view [s]
N = floor(prf * t_total);                     % Number of data points 
PRI = 1 / prf;
dt = PRI;                                   % Time between points ~ integration time of a 'scene'
time = linspace(0, t_total, N);

%% Generate Orbit Data That Pass Through FOV
% Define Debris 
h = 630;
e = 0;
i = 90; 
RAAN = 90;
w = 0;
f0 = 0; 
tspan = 0:dt:t_total;
Debris_kep = [Re+h, e, deg2rad(i), deg2rad(RAAN), deg2rad(w), deg2rad(f0)];
[SC_rv, SC_t] = Debris_Propagtor(SC_kep, tspan);
[Debris_rv, Debris_t] = Debris_Propagtor(Debris_kep, tspan);

% Define the Debris position relative to the satellite
Debris_r = Debris_rv(:, 1:3);
Debris_v = Debris_rv(:, 4:6);
SC_r = SC_rv(:, 1:3);                
SC_v = SC_rv(:, 4:6);                

% Convert ECI to LVLH and gate datapoints 
for i = 1:length(Debris_r(:,1))
    [Debris_dr_LVLH, Debris_dv_LVLH] = rotate_ECI2LVLH(Debris_r(i,:)', Debris_v(i,:)', SC_r(i,:)', SC_v(i,:)');
    Debris_state_LVLH(i, :) = [Debris_dr_LVLH; Debris_dv_LVLH];
    % Convert to spherical LVLH 
    Debris_state_LVLH_spherical(i, :) = CartLVLH2SphLVLH(Debris_state_LVLH(i, :));

    % Gate by Range (<50 km) and FOV for Azimuth and Elevation
    if Debris_state_LVLH_spherical(i,1) > 50
        Debris_state_LVLH_spherical(i, :) = NaN;
    end   
    if (deg2rad(-1.26) > Debris_state_LVLH_spherical(i,2)) || (Debris_state_LVLH_spherical(i,2) > deg2rad(1.26))
        Debris_state_LVLH_spherical(i, :) = NaN;
    end
    if (deg2rad(-15) > Debris_state_LVLH_spherical(i,3)) || (Debris_state_LVLH_spherical(i,3) > deg2rad(15))
        Debris_state_LVLH_spherical(i, :) = NaN;
    end 
end

% Keep only the measurable points (remove azimuth and elevation rates for radar)
Debris_state_LVLH_spherical_Measurable = Debris_state_LVLH_spherical;
Debris_state_LVLH_spherical_Measurable(:, 5) = NaN;
Debris_state_LVLH_spherical_Measurable(:, 6) = NaN;
Az = Debris_state_LVLH_spherical_Measurable(:, 2);
El = Debris_state_LVLH_spherical_Measurable(:, 3);

% Compute finite differences for azimuth and elevation rates
Az_rate = numerical_derivative(Az, dt);
El_rate = numerical_derivative(El, dt);
Debris_state_LVLH_spherical_Measurable(:, 5) = Az_rate;
Debris_state_LVLH_spherical_Measurable(:, 6) = El_rate;

Debris_error_LVLH_spherical = [erange eangle eangle erange_rate eaz_rate eel_rate];

% Convert the measurable spherical LVLH measurements to Cartesian LVLH for EKF
for i = 1:length(Debris_state_LVLH_spherical_Measurable(:,1))
    [Debris_state_LVLH_est(i, :), Debris_error_LVLH(:, :, i)] = SphLVLH2CartLVLH(Debris_state_LVLH_spherical_Measurable(i, :), Debris_error_LVLH_spherical(1, :));
end

% Define initial state and covariance for the EKF (in Cartesian LVLH)
x0 = Debris_state_LVLH(1, :);            % True initial state [position and velocity]
P0 = Debris_error_LVLH(:, :, 1);         % Initial covariance (Cartesian LVLH)

% Measurements for the EKF (using only range, az, el, and range rate from spherical LVLH)
Debris_state_LVLH_spherical_Measurable_EKF = Debris_state_LVLH_spherical_Measurable(:, 1:6);
measurements = Debris_state_LVLH_spherical_Measurable_EKF;  
t = Debris_t;
a = Re + h_SC;

%% Grid Search for Optimal Process Noise Matrix Q using EKF_IODV3
% Use the "true" state history as Debris_state_LVLH (in Cartesian LVLH)
[bestQ, bestResults] = gridSearch_EKF_IODV3(x0, P0, measurements, t, PRI, a, Debris_state_LVLH);

% Display the best Q and corresponding error metrics
fprintf('Best Q matrix found via grid search:\n');
disp(bestQ);
fprintf('Best scaling factors: alpha = %.2f, beta = %.2f\n', bestResults.alpha, bestResults.beta);
fprintf('RMSE per state dimension: %s\n', mat2str(bestResults.RMSE));
fprintf('Final covariance trace: %e\n', bestResults.covMetric);

%% (Optional) Use the Best Q to Run the EKF Once More
[X_hist_opt, P_hist_opt] = EKF_IODV3(x0, P0, measurements, t, PRI, a, bestQ);

%% Plots to Determine "Optimal" Number of Observations
% (Plotting the evolution of the filter covariance)
N_obs = size(P_hist_opt, 3);
threeSigma = zeros(6, N_obs);
for k = 1:N_obs
    sigma = sqrt(diag(P_hist_opt(:, :, k)));
    threeSigma(:, k) = 3 * sigma;
end

state_labels = {'x', 'y', 'z', 'v_x', 'v_y', 'v_z'};
figure;
for state = 1:6
    subplot(3, 2, state);
    plot(1:N_obs, threeSigma(state, :), 'LineWidth', 1.5);
    xlabel('Observation Number');
    ylabel('3\sigma Uncertainty (km)');
    title(sprintf('State %s', state_labels{state}));
    grid on;
end

%% The rest of the script follows...
% (Monte Carlo Simulation and additional propagation/plotting as in your original script)
% ...

%% End of Main Script
