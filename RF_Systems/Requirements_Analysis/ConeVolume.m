
clc; clear; close all
% FrustumVolume2DPlot_HalfAngle_Varying_h0.m
% This script visualizes the volume of a conical frustum as a function of 
% the frustum height (h) and the half cone angle (θ). Here, the lower 
% height h0 is defined as a function of θ:
%
%     h0(θ) = 875/tan(θ)
%
% The volume of the full cone at any height z is:
%     V(z) = (1/3)*pi*z^3*tan^2(θ)
%
% The frustum volume between z = h0(θ) and z = h0(θ)+h is:
%     V = (1/3)*pi*tan^2(θ)*[(h0(θ)+h)^3 - h0(θ)^3]
%
% Partial derivatives:
%   dV/dh = pi*tan^2(θ)*(h0(θ)+h)^2
%
%   dV/dθ = (1/3)*pi*sec^2(θ)*[ 2*tan(θ)*((h0(θ)+h)^3 - h0(θ)^3) - 3*875*((h0(θ)+h)^2 - h0(θ)^2) ]
%

%% Define ranges for frustum height h and half cone angle θ
h = linspace(0, 200, 2000);            % Frustum height from 0 to 8 units
theta = linspace(0.1, deg2rad(10), 2000);     % Half cone angle from ~0.1 rad to 45° (avoid theta=0 to prevent division by zero)

% Create a grid of (h, θ) values
[H, T] = meshgrid(h, theta);

%% Define h0 as a function of theta (for each point in T)
h_min = 1.75 ./ tan(T);  % h0 now is an array with the same dimensions as T

%% Compute the Frustum Volume 
V = (1/3) * pi * ((h.^3.*tan(T).^2) - ((0.875)^3./tan(theta)) );

%% Define discrete levels for contour plots
numBands = 200;
levelsV = linspace(min(V(:)), max(V(:)), numBands+1);

figure;

% Frustum Volume V plot
contourf(H, rad2deg(T), V, levelsV, 'LineColor', 'none');
shading interp;
xlabel('Frustum Height h');
ylabel('Half Cone Angle θ (deg)');
title('Frustum Volume V');
colorbar;

