%% Main IOD script: Take Measurements a Propagate forward for a revisit time
clear; clc; close all
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

% SIDE LOOKING GEOMETRY

%% Constant 
mu = 3.986e5; % Earth's gravitational parameter (km^3/s^2)
Re = 6371;      %[km]
% SC parameters TBC 
h_SC = 590;%km
e = 0;
i = 90;
Om = 90;
om = 0;
theta = 0;
c = 3e8;
f = 5e9;
lambda = c/f;
SC_kep= [Re + h_SC, e, deg2rad(i), deg2rad(Om), deg2rad(om), deg2rad(theta)]; % [a e i Om om theta]

%% SAR sizing/uncertanity parameters [Herrick-Gibbs 1 satellite] (inputs)
Nsub = 8;
v_rel = 7.5;            % Maximum Relative Velocity of the Debris for this design [km]
FOV = 5.873;           % ScanSAR Half cone FOV [deg]
look_angle = 20;
FOV_strip = FOV/Nsub;   % FoV per strip [deg]
Range_max = 68.24;         % [km]
Optical_SW = 100;       % [km] The optical Swath width TBD
prf = 9375;
%Angle Errors at design point
erange = 45.3e-3;        % Range Error [km]
eangle = deg2rad(0.01);           % Angle error [Degrees]
erange_rate = 0.005e-3;
ediff = 0.01;
eaz_rate = sqrt(2*(eangle^2));  %Finite difference error (taylor, introduction to error analysis)
eel_rate = eaz_rate;
%% Precompute Range, Az, El for verification

SW_max = 2*Range_max*(tand(FOV_strip + look_angle) - tand(look_angle));   %strip swath width [km]
t_total = SW_max/(8*v_rel);        % Timespan in view



N = floor(prf*t_total);      % Number of data points 
PRI = 1/prf;
dt = PRI;                    % The difference between points is roughly the integration time of a 'scene'
time = linspace(0,t_total,N);
%% Generate orbit data that passess through FOV

% Define Debris 
h = 610;
e = 0;
i = 90; 
RAAN = 90;
w = 0;
f0 = 0; 
tspan = 0:dt:t_total;
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
    [Debris_state_LVLH_spherical(i,:)] = CartLVLH2SphLVLH(Debris_state_LVLH(i,:));  %[km] and [rad]

    %Assess if the Range is <50km and the az<FOV_strip deg and the el>15 deg
    if Debris_state_LVLH_spherical(i,1) > 50
        Debris_state_LVLH_spherical(i,:) = NaN;
    end   
    %Assess if Debris is in Azimuth FOV
    if (deg2rad(-1.26) > Debris_state_LVLH_spherical(i,2)) || (Debris_state_LVLH_spherical(i,2) > deg2rad(1.26))
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
Debris_state_LVLH_spherical_Measurable(:,5) = Az_rate; %[deg]
Debris_state_LVLH_spherical_Measurable(:,6) = El_rate; %[deg]


Debris_error_LVLH_spherical = [erange eangle eangle erange_rate eaz_rate eel_rate]; % [km] and [rad]

%Kalman Filter requires cartesian measurements so convert from spherical to
%cartesian -> keep LVLH
Debris_state_LVLH_est = zeros(1,6);
for i = 1:length(Debris_state_LVLH_spherical_Measurable(:,1))
     [Debris_state_LVLH_est(i,:), Debris_error_LVLH(:,:,i)] = SphLVLH2CartLVLH(Debris_state_LVLH_spherical_Measurable(i,:), Debris_error_LVLH_spherical(1,:));
end

x0 = Debris_state_LVLH(1,:);                                %Require Cartesian LVLH
P0 = Debris_error_LVLH(:,:,1);                                  %Require Cartesian LVLH
Debris_state_LVLH_spherical_Measurable_EKF = Debris_state_LVLH_spherical_Measurable(:,1:6); 
measurements = Debris_state_LVLH_spherical_Measurable_EKF;      %Require Spherical LVLH
t = Debris_t;
a = Re+h_SC;
[X_hist, P_hist] = EKF_IODV2(x0, P0, measurements, t, PRI, a);

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
state_labels = {'x', 'y', 'z', 'v_x', 'v_y', 'v_z'};
figure;
for state = 1:6
    subplot(3, 2, state); % 3 rows and 2 columns layout
    plot(1:N, threeSigma(state,:), 'LineWidth', 1.5);
    xlabel('Observation Number');
    ylabel('Uncertainty (km)');
    title(sprintf('3 sigma uncertainty of state %s', state_labels{state}));
    xlim([0 290]);
    ylim([0 2]);
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
    xlim([0 290]);
    grid on;
end

%Min Time in FOV = scene time * number of observations 


%% Monte Carlo Simulation 
% This program propagates the uncertainty distribution of a 
% satellite in Earth orbit using a Monte-Carlo simulation. In addition,
% it records the covariance matrix at every time step.


%% Covariance Propagator

%Propagate the Shrunken Covariance using the State transition Matrix 
t0 = 0;                 % Initial time [s]
t_final = 24*3600;      % End time [s]
num_steps = 10000;      % Number of time points to record
tspan = linspace(t0, t_final, num_steps);  % Time vector for output

X0 = Debris_rv(1,:);         % Debris Position and velocity ECI
P0 = P_hist(:,:,end);        % Shrunken Covariance Matrix

%% NOW OPEN MC_Simulation.m in the 'Github codes' folder










% %% Inputs
% t0 = 0;                 % Initial time [s]
% t_final = 2*3600;      % End time [s]
% num_steps = 10000;        % Number of time points to record
% t_vector = linspace(t0, t_final, num_steps);  % Time vector for output
% 
% % Initial mean state [x, y, z, xdot, ydot, zdot]
% SC_r = SC_rv(1,1:3);
% SC_v = SC_rv(1,4:6);
% Debris_r = x0(1:3);
% Debris_v = x0(4:6);
% 
% 
% [r_B, v_B] = rotate_LVLH2ECI(x0(1:3)', x0(4:6)', SC_r', SC_v');
% x0 = Debris_rv(1,:);                %Assume  
% % Initial covariance matrix (6x6) [km and km/s]
% p0 = P_hist(:,:,end);    % Provided covariance at convergence point
% 
% n = length(x0);         % State vector dimension


% %% Preallocate arrays for state trajectories and covariance history
% 
% [p_out, x_traj, STM] = phi(t_vector, x0, p0);
% x_all = x_traj;
% 
% 
% P_history = p_out;
% 
% %% Plotting: Evolution of Standard Deviations in x, y, and z over Time
% % Preallocate arrays for sigma_x, sigma_y, and sigma_z
% sigma_x = zeros(1, num_steps);
% sigma_y = zeros(1, num_steps);
% sigma_z = zeros(1, num_steps);
% 
% for j = 1:num_steps
%     sigma_x(j) = 3*sqrt(P_history(1,1,j));
%     sigma_y(j) = 3*sqrt(P_history(2,2,j));
%     sigma_z(j) = 3*sqrt(P_history(3,3,j));
% end
% 
% figure;
% plot(t_vector/3600, sigma_x, 'LineWidth', 1.5, 'DisplayName', '\sigma_x');
% hold on;
% plot(t_vector/3600, sigma_y, 'LineWidth', 1.5, 'DisplayName', '\sigma_y');
% plot(t_vector/3600, sigma_z, 'LineWidth', 1.5, 'DisplayName', '\sigma_z');
% xlabel('Time (hrs)');
% ylabel('Uncertainty (km)');
% title('Evolution of Uncertainty (\sigma) in x, y, and z');
% legend('show');
% ylim([0 1000])
% grid on;
% 
% % %% Plotting the Final Distribution using Scatter Histograms
% % 
% % % Plot y versus x
% % figure('Name', 'f1')
% % s1 = scatterhist(squeeze(x_all(:, end, 1)), squeeze(x_all(:, end, 2)), ...
% %     'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1);
% % xlabel('x (km)')
% % ylabel('y (km)')
% % 
% % % Plot z versus x
% % figure('Name', 'f2')
% % s2 = scatterhist(squeeze(x_all(:, end, 1)), squeeze(x_all(:, end, 3)), ...
% %     'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1);
% % xlabel('x (km)')
% % ylabel('z (km)')
% % 
% % % Plot ydot versus xdot
% % figure('Name', 'f3')
% % s3 = scatterhist(squeeze(x_all(:, end, 4)), squeeze(x_all(:, end, 5)), ...
% %     'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1);
% % xlabel('$\dot{x}$ (km/s)', 'Interpreter', 'latex')
% % ylabel('$\dot{y}$ (km/s)', 'Interpreter', 'latex')
% % 
% % % Plot zdot versus xdot
% % figure('Name', 'f4')
% % s4 = scatterhist(squeeze(x_all(:, end, 4)), squeeze(x_all(:, end, 6)), ...
% %     'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1);
% % xlabel('$\dot{x}$ (km/s)', 'Interpreter', 'latex')
% % ylabel('$\dot{z}$ (km/s)', 'Interpreter', 'latex')
% % 
% % % Combine the four scatter plots into one figure
% % figure('Name', 'Combined Plots')
% % u1 = uipanel('Position', [0 .5 .5 .5]);
% % set(s1, 'Parent', u1)
% % u2 = uipanel('Position', [.5 .5 .5 .5]);
% % set(s2, 'Parent', u2)
% % u3 = uipanel('Position', [0 0 .5 .5]);
% % set(s3, 'Parent', u3)
% % u4 = uipanel('Position', [.5 0 .5 .5]);
% % set(s4, 'Parent', u4)
% % close f1 f2 f3 f4
% 
% toc
% 
% %% Modified Propagation Function: phi
% 
% function [p_out, x_out, STM] = phi(t_vector, x0, p0)
%     % Propagates the two-body dynamics and the STM
%     % t_vector: vector of times to record the solution
%     % x0: initial state [r; v] (6x1)
%     % Returns:
%     %   t_out: time vector from ode45
%     %   x_out: state trajectory [num_steps x 6]
%     %   STM: state transition matrices at each time step [6 x 6 x num_steps]
%     dt = mean(diff(t_vector));
%     x_out = x0';
%     p_out = p0;
% 
%     for i = 2:(length(t_vector))
%         A = sm(x_out(:,i-1));
%         STM = expm(A*dt);
% 
%         x_out(:,i) = STM*x_out(:,i-1);
%         p_out(:,:,i) = STM*p_out(:,:,i-1) * STM';
%     end
%     x_out = x_out';
% 
% 
% end
% 
% 
% 
% % function [t_out, x_out, STM] = phi(t_vector, x0)
% %     % Propagates the two-body dynamics and the STM
% %     % t_vector: vector of times to record the solution
% %     % x0: initial state [r; v] (6x1)
% %     % Returns:
% %     %   t_out: time vector from ode45
% %     %   x_out: state trajectory [num_steps x 6]
% %     %   STM: state transition matrices at each time step [6 x 6 x num_steps]
% % 
% %     % Initial STM is the 6x6 identity matrix.
% %     phi0 = eye(6);
% %     % Flatten the STM to a vector.
% %     phi0_vec = reshape(phi0, [],1);
% % 
% %     % Augment the initial state with the STM vector.
% %     y0 = [x0'; phi0_vec];
% % 
% %     % Set ODE options (optional: adjust tolerances as needed)
% %     options = odeset('RelTol',1e-10,'AbsTol',1e-12);
% % 
% %     % Solve the augmented differential equation using ode45.
% %     [t_out, y_out] = ode45(@(t,y) two_body_aug(t,y), t_vector, y0, options);
% % 
% %     % Extract the state trajectory.
% %     x_out = y_out(:,1:6);
% % 
% %     % Number of time steps.
% %     num_steps = length(t_out);
% %     % Preallocate the 3D array for the STM history.
% %     STM = zeros(6,6,num_steps);
% % 
% %     % Extract the STM at each time step.
% %     for k = 1:num_steps
% %         phi_vec = y_out(k, 7:end);
% %         STM(:,:,k) = reshape(phi_vec, 6, 6);
% %     end
% % end
% % 
% % function dydt = two_body_aug(~, y)
% %     % Computes the derivatives for the augmented state vector [x; phi_vec]
% %     % where x is the state (6x1) and phi is the 6x6 STM (flattened to 36x1).
% %     % Gravitational parameter (e.g., for Earth: mu = 398600 km^3/s^2)
% %     mu = 398600;  % adjust as needed
% %     % Extract state vector.
% %     state = y(1:6);
% %     % Reshape the remaining 36 elements into the 6x6 STM.
% %     phi = reshape(y(7:end), 6, 6);
% % 
% %             % Component
% %     x = state(1);
% %     y = state(2);
% %     z = state(3);
% %     % Vector
% %     r = [x y z]';
% %     % Magnitude
% %     R = norm(r);
% % 
% % 
% %     % Component
% %     dx = state(4);
% %     dy = state(5);
% %     dz = state(6);
% %     % Vector
% %     v = [dx dy dz]';
% %     % Magnitude
% %     V = norm(v);
% % 
% % 
% %     % Component
% %     ddx = -mu*x/R^3;
% %     ddy = -mu*y/R^3;
% %     ddz = -mu*z/R^3;
% %     % Vector
% %     a = [ddx ddy ddz]';
% % 
% %     xdot = [v;a]; % Return state vector for next step
% % 
% % 
% % 
% %      %Extract the state
% %     rx = state(1);
% %     ry = state(2);
% %     rz = state(3);
% % 
% %     r0_norm = norm(state(1:3));
% %     coef1 = 3*mu/(2*r0_norm^5);
% %     coef2 = mu/(2*r0_norm^3);
% % 
% %     %Form 2BP dynamics matrix
% %     A_21 = [coef1*rx^2 - coef2, coef1*rx*ry, coef1*rx*rz;
% %             coef1*rx*ry, coef1*ry^2 - coef2, coef1*ry*rz;
% %             coef1*rx*rz, coef1*ry*rz, coef1*rz^2-coef1];
% % 
% %     A = [zeros(3), eye(3); A_21, zeros(3)];
% % 
% %     % Compute the derivative of the STM.
% %     phi_dot = A * phi;
% % 
% %     % Flatten phi_dot to a column vector.
% %     phi_dot_vec = reshape(phi_dot, [], 1);
% % 
% %     % Combine the derivatives.
% %     dydt = [xdot; phi_dot_vec];
% % end
% % 
