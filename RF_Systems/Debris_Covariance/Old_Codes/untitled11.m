%% Main IOD script: Take Measurements a Propagate forward for a revisit time
clear; clc;
%Purpose and Methodology: 
% 1) Input SAR sizing/uncertanity parameters
% 2) Conducts coordinate transformations from spherical LVLH to Cartesian
% ECI
% 3b) If not Herricks Gibbs; Assume the angle rate can be extracted by
% finite differences (verified)
% 4) Pass Measurement Uncertanity through EKF to get Final measurement
% uncertanity/covariance (not verified)
% 5) Estimate the revisit time required based on the Projected Area of the
% Optical Instrument and covariance propagation (verified)

%% Constant 
mu = 3.986e5; % Earth's gravitational parameter (km^3/s^2)
Re = 6371;      %[km]
% SC parameters TBC 
h = 590;%km
e = 0;
i = 90;
Om = 90;
om = 0;
theta = 0;
c = 3e8;
f = 5e9;
lambda = c/f;
SC_kep= [Re + h, e, deg2rad(i), deg2rad(Om), deg2rad(om), deg2rad(theta)]; % [a e i Om om theta]

%% SAR sizing/uncertanity parameters [Herrick-Gibbs 1 satellite] (inputs)

v_rel = 7e3;            % Maximum Relative Velocity of the Debris for this design 
FOV = 11;               % Half cone FOV [deg]
Range_max = 50;         % [km]
Optical_SW = 100;       % [km] The optical Swath width TBD

%Angle Errors at design point
erange = 30.7e-3;        % Range Error [km]
eangle = deg2rad(0.1);           % Angle error [Degrees]
erange_rate = 0.005e-3;
ediff = 0.01;
eaz_rate = sqrt(2*(eangle^2));  %Finite difference error (taylor, introduction to error analysis)
eel_rate = eaz_rate;
%% Precompute Range, Az, El for verification

SW_max = 2*Range_max*tand(FOV);   %Half cone swath width
t_swath = 1e-3;                 %The time between pulses limits the number of data points

t_total = SW_max/v_rel;        % Timespan in view
N = floor(t_total/t_swath);      % Number of data points 
dt = t_swath;                    % The difference between points is roughly the integration time of a 'scene'
time = linspace(0,t_total,N);
%% Generate orbit data that passess through FOV

% Define Debris 
h = 630;
e = 0;
i = 90; 
RAAN = 90;
w = 0;
f0 = 0; 
tspan = 0:1e-3:2;
Debris_kep = [Re+h,e,deg2rad(i),deg2rad(RAAN), deg2rad(w),deg2rad(f0)];
[SC_rv,SC_t] = Debris_Propagtor(SC_kep,tspan);
[Debris_rv, Debris_t]  = Debris_Propagtor(Debris_kep,tspan);


%Define the Debris position relative to the satellite
Debris_r = Debris_rv(:,1:3);
Debris_v = Debris_rv(:,4:6);
SC_r = SC_rv(:,1:3);                
SC_v = SC_rv(:,4:6);                

%Convert ECI to LVLH and gate datapoints 

for i = 1:length(Debris_r(:,1))
    [Debris_dr_LVLH, Debris_dv_LVLH] = rotate_ECI2LVLH(Debris_r(i,:)', Debris_v(i,:)', SC_r(i,:)', SC_v(i,:)');
    Debris_state_LVLH(i,:) = [Debris_dr_LVLH;Debris_dv_LVLH];
    %Convert to spherical LVLH 
    [Debris_state_LVLH_spherical(i,:)] = CartLVLH2SphLVLH(Debris_state_LVLH(i,:));

    %Assess if the Range is <50km and the az<11 deg and the el>15 deg
    if Debris_state_LVLH_spherical(i,1) > 50
        Debris_state_LVLH_spherical(i,:) = NaN;
    end   
    %Assess if Debris is in Azimuth FOV
    if (deg2rad(-11) > Debris_state_LVLH_spherical(i,2)) || (Debris_state_LVLH_spherical(i,2) > deg2rad(11))
        Debris_state_LVLH_spherical(i,:) = NaN;
    end
    %Assess if the Debris is in Elevation FOV
    if (deg2rad(-15) > Debris_state_LVLH_spherical(i,3)) || (Debris_state_LVLH_spherical(i,3) > deg2rad(15))
        Debris_state_LVLH_spherical(i,:) = NaN;
    end 
end


%The non NaNs have entered the sensor view, simulate a streak observable by the radar system 
%A set of points [Range, Az, El] and the range rate through doppler shift
%Therefore remove az rate and el rate as that cant be measured.

Debris_state_LVLH_spherical_Measurable = Debris_state_LVLH_spherical;
Debris_state_LVLH_spherical_Measurable(:,5) = NaN;
Debris_state_LVLH_spherical_Measurable(:,6) = NaN;
Az = Debris_state_LVLH_spherical_Measurable(:,2);
El = Debris_state_LVLH_spherical_Measurable(:,3);

%Conduct finite differences between percieved points within 'image'
Az_rate = numerical_derivative(Az, dt);
El_rate = numerical_derivative(El, dt);
Debris_state_LVLH_spherical_Measurable(:,5) = Az_rate;
Debris_state_LVLH_spherical_Measurable(:,6) = El_rate;


Debris_error_LVLH_spherical = [erange eangle eangle erange_rate eaz_rate eel_rate];

%Kalman Filter requires cartesian measurements so convert from spherical to
%cartesian -> keep LVLH
Debris_state_LVLH_est = zeros(1,6);
for i = 1:length(Debris_state_LVLH_spherical_Measurable(:,1))
     [Debris_state_LVLH_est(i,:), Debris_error_LVLH(:,:,i)] = SphLVLH2CartLVLH(Debris_state_LVLH_spherical_Measurable(i,:), Debris_error_LVLH_spherical(1,:));
end

x0 = Debris_state_LVLH(1,:);                                %Require Cartesian LVLH
P0 = Debris_error_LVLH(:,:,1);                              %Require Cartesian LVLH
Debris_state_LVLH_spherical_Measurable_EKF = Debris_state_LVLH_spherical_Measurable(:,1:4); 
measurements = Debris_state_LVLH_spherical_Measurable_EKF;      %Require Spherical LVLH
t = Debris_t;

[X_hist, P_hist] = EKF_IOD(x0, P0, measurements, t);

%% Plots to determine "Optimal" number of observations

%We want the minimum number of observations as that means the debris doesnt
%need to linger in the FOV for as long
%But We want the maximum revisit time so mission don't kill me

%so we want to pick the point before deminishing returns set in

%Plot the number of observations vs covariance Shrinkage, select number of observations.
%Min Time in FOV = scene time * number of observations 

%Extract the 
% Assume covHistory is your 6x6xN covariance history variable
N = size(P_hist, 3);  % Number of observations

% Preallocate an array to store 3*sigma for each state over time
threeSigma = zeros(6, N);

% Compute 3-sigma values for each observation
for k = 1:N
    % Extract the 6x6 covariance matrix at iteration k
    P = P_hist(:,:,k);
    
    % Compute the standard deviation (square root of diagonal variances)
    sigma = sqrt(diag(P));
    
    % Multiply by 3 to get 3-sigma values
    threeSigma(:, k) = 3 * sigma;
end

% Create a figure with 6 subplots, one for each state
% Figure 1: Plot 3-sigma uncertainties for each state
state_labels = {'x', 'y', 'v_x', 'v_y', 'v_z', 'v_z'};
figure;
for state = 1:6
    subplot(3, 2, state); % 3 rows and 2 columns layout
    plot(1:N, threeSigma(state,:), 'LineWidth', 1.5);
    xlabel('Observation Number');
    ylabel('Uncertainty (km)');
    title(sprintf('3 sigma uncertainty of state %s', state_labels{state}));
    xlim([0 100]);
    grid on;
end

% Figure 2: Plot the rate of change of 3-sigma uncertainties for each state
figure;
for state = 1:6
    subplot(3, 2, state); % 3 rows and 2 columns layout
    % Compute the derivative using the gradient function
    d_threeSigma = gradient(threeSigma(state,:));
    plot(1:N, d_threeSigma, 'LineWidth', 1.5);
    xlabel('Observation Number');
    ylabel('Rate of Change (km/obs)');
    title(sprintf('Rate of change of 3 sigma uncertainty for state %s', state_labels{state}));
    xlim([0 100]);
    grid on;
end

%Min Time in FOV = scene time * number of observations 


%% Propagate Covariance through time till error limit
% Define propagation settings and initial conditions (as in your code)
mu = 3.986e5;
tspan = [0 10*24*3600];        % Propagation time: 10 days in seconds
P0 = P_hist(:,:,120);           % Covariance matrix when debris leaves the FOV
x0 = X_hist(end,:);            % State when leaving FOV

%Convert x0 LVLH to ECI 


% Define Optical Swath Width values (error limit in km)
Optical_SW = [0.5, 0.8, 1, 5, 10, 100];

% Prepare a figure for the 3-sigma uncertainty curves.
figure;
hold on;
colors = lines(length(Optical_SW));

% Loop through each Optical_SW value
for i = 1:length(Optical_SW)
    % Call the propagation function.
    % (Make sure your Covariance_Prop function is modified to output Time_out and Cov_out.)
    [final_state, final_cov, t_limit, Time_out, Cov_out] = Covariance_Prop(x0, P0, tspan, Optical_SW(i));
    
    % Compute uncertainties from the covariance history.
    % Here we extract the position uncertainties for x, y, and z.
    sigma_x = sqrt(Cov_out(:,7));   % Column 7 corresponds to variance in x
    sigma_y = sqrt(Cov_out(:,14));  % Column 14 corresponds to variance in y
    sigma_z = sqrt(Cov_out(:,21));  % Column 21 corresponds to variance in z
    
    % For this plot, we choose the x-position uncertainty as representative.
    threeSigma = 3 * sigma_x;
    
    % Plot the 3-sigma uncertainty versus time (converted to hours).
    plot(Time_out/3600, threeSigma, 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Optical SW = %.2f km', Optical_SW(i)));
     
    % Mark the revisit time on the curve with a marker.
    idx = find(Time_out >= t_limit, 1);
    plot(t_limit/3600, threeSigma(idx), 'o', 'Color', colors(i,:), ...
         'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
end

xlabel('Time [hrs]');
ylabel('3\sigma Uncertainty (km)');
title('3\sigma Uncertainty vs Time for Different Optical SW');
legend('Location', 'best');
grid on;
%% Propagate Covariance to get revisit time
%Plot revisit time vs number of observations