tic; close all; clc

%% Monte-Carlo Propagation of Satellite Uncertainty with Covariance History Recording
% This program propagates the uncertainty distribution of a 
% satellite in Earth orbit using a Monte-Carlo simulation. In addition,
% it records the covariance matrix at every time step.

%% Inputs
N = 100;                % Number of random samples
t0 = 0;                 % Initial time [s]
t_final = 24*3600;      % End time [s]
num_steps = 100;        % Number of time points to record
t_vector = linspace(t0, t_final, num_steps);  % Time vector for output

% Initial mean state [x, y, z, xdot, ydot, zdot]
m0 = [3.9583e4, -1.4667e4, 0.1035e4, 1.0583, 2.8815, 0.0842];

% Initial covariance matrix (6x6) [km and km/s]
p0 = (P_hist(:,:,125) + P_hist(:,:,125)') / 2;    % Provided covariance at index 66

n = length(m0);         % State vector dimension

%% Generate random samples from the initial Gaussian distribution
x0_samples = mvnrnd(m0, p0, N);

%% Preallocate arrays for state trajectories and covariance history
% x_all: [N samples x num_steps x n state elements]
x_all = zeros(N, num_steps, n);

% Propagate each Monte Carlo sample and store its full trajectory
for i = 1:N
    [t_out, x_traj] = phi(t_vector, x0_samples(i, :));
    x_all(i, :, :) = x_traj;
end

% Compute the covariance history at each time step
% P_history: [n x n x num_steps]
P_history = zeros(n, n, num_steps);
for j = 1:num_steps
    % Extract the states at time t_vector(j) for all samples (an N x n matrix)
    sample_states = squeeze(x_all(:, j, :));  
    P_history(:, :, j) = cov(sample_states);
end

%% Compute the final covariance matrix at t_final and output the ensemble mean
P_final = cov(squeeze(x_all(:, end, :)));
disp('Final Covariance Matrix:');
disp(P_final);

final_mean = mean(squeeze(x_all(:, end, :)), 1);
fprintf('Monte Carlo Simulation mean: \n x = %.3e, y = %.3e, z = %.3e, xdot = %.3e, ydot = %.3e, zdot = %.3e\n', final_mean);

%% Plotting: Evolution of Standard Deviations in x, y, and z over Time
% Preallocate arrays for sigma_x, sigma_y, and sigma_z
sigma_x = zeros(1, num_steps);
sigma_y = zeros(1, num_steps);
sigma_z = zeros(1, num_steps);

for j = 1:num_steps
    sigma_x(j) = 3*sqrt(P_history(1,1,j));
    sigma_y(j) = 3*sqrt(P_history(2,2,j));
    sigma_z(j) = 3*sqrt(P_history(3,3,j));
end

figure;
plot(t_vector/3600, sigma_x, 'LineWidth', 1.5, 'DisplayName', '\sigma_x');
hold on;
plot(t_vector/3600, sigma_y, 'LineWidth', 1.5, 'DisplayName', '\sigma_y');
plot(t_vector/3600, sigma_z, 'LineWidth', 1.5, 'DisplayName', '\sigma_z');
xlabel('Time (hrs)');
ylabel('Uncertainty (km)');
title('Evolution of Uncertainty (\sigma) in x, y, and z');
legend('show');
ylim([0 1000])
grid on;

%% Plotting the Final Distribution using Scatter Histograms

% Plot y versus x
figure('Name', 'f1')
s1 = scatterhist(squeeze(x_all(:, end, 1)), squeeze(x_all(:, end, 2)), ...
    'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1);
xlabel('x (km)')
ylabel('y (km)')

% Plot z versus x
figure('Name', 'f2')
s2 = scatterhist(squeeze(x_all(:, end, 1)), squeeze(x_all(:, end, 3)), ...
    'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1);
xlabel('x (km)')
ylabel('z (km)')

% Plot ydot versus xdot
figure('Name', 'f3')
s3 = scatterhist(squeeze(x_all(:, end, 4)), squeeze(x_all(:, end, 5)), ...
    'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1);
xlabel('$\dot{x}$ (km/s)', 'Interpreter', 'latex')
ylabel('$\dot{y}$ (km/s)', 'Interpreter', 'latex')

% Plot zdot versus xdot
figure('Name', 'f4')
s4 = scatterhist(squeeze(x_all(:, end, 4)), squeeze(x_all(:, end, 6)), ...
    'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1);
xlabel('$\dot{x}$ (km/s)', 'Interpreter', 'latex')
ylabel('$\dot{z}$ (km/s)', 'Interpreter', 'latex')

% Combine the four scatter plots into one figure
figure('Name', 'Combined Plots')
u1 = uipanel('Position', [0 .5 .5 .5]);
set(s1, 'Parent', u1)
u2 = uipanel('Position', [.5 .5 .5 .5]);
set(s2, 'Parent', u2)
u3 = uipanel('Position', [0 0 .5 .5]);
set(s3, 'Parent', u3)
u4 = uipanel('Position', [.5 0 .5 .5]);
set(s4, 'Parent', u4)
close f1 f2 f3 f4

toc

%% Modified Propagation Function: phi
function [t_out, x_out] = phi(t_vector, x0)
    % Solves the differential equation for the satellite's motion
    % over the specified time vector.
    %   t_vector: vector of times at which to record the solution
    %   x0: initial state vector [x, y, z, xdot, ydot, zdot]
    mu = 398600; % Earth's gravitational parameter (km^3/s^2)
    
    % Function to compute the magnitude of the position vector
    r = @(x) sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    
    % Use ode45 to solve the equations of motion at the times specified by t_vector
    [t_out, x_out] = ode45(@(ti, yi) [yi(4); yi(5); yi(6); ...
                            -mu * yi(1) / r(yi)^3; -mu * yi(2) / r(yi)^3; -mu * yi(3) / r(yi)^3], ...
                           t_vector, x0);
end
