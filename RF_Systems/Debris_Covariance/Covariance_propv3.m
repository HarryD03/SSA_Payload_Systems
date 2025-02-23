%% Observation uncertanity to covariance Propagation
% Used to define revisit time requirement
% By Harry Dawson

%Assumptionss:
%Observation is taken along radial vector (in LVLH frame), therefore
%spherical to LVLH frame is
%2BP dynamics only taken into account - polar orbit so less RAAN drift due
%to J2 and Low CSA of debris means low drag. Therefore 2BP dominates 

clear; clc; close all

%% Constants
mu = 3.986e5; % Earth's gravitational parameter (km^3/s^2)
Re = 6371;      %[km]
%% Inputs

%SC parameters TBC 
h = 500;%km
e = 0;
i = 80;
Om = 80;
om = 0;
theta = 0;

%Choose maximim values for Best case - FOV boundary
%Choose Minimum values for Worst Case - See "Robust initial orbit
%determination for short-arc Doppler radar observations"

%TBC: Choice dependant on the FOV from the payload
r = 100; % km 
az = 1;   %deg
elevation = 1; % deg

%Choice based on design choice
sigma_r = 10e-3; %km
sigma_ang = 0.1; %deg

%Spacecraft Cartesian Coordinates ECI
%There will be some error assocaited with this TO ADD

% Get Keplerian Orbit from mission 
SC_kep= [Re + h, e, deg2rad(i), deg2rad(Om), deg2rad(om), deg2rad(theta)]; % [a e i Om om theta]
SC_X = kep2car(SC_kep, mu);


%% Observation vectors into ECI Cartesian vectors 
% convert range, azimuth and elevation from one s/c and debris to ECI
% cartesian coordinates 

%Error rotation
sigma_az = deg2rad(sigma_ang);
sigma_el = deg2rad(sigma_ang);

P_sphere = diag([sigma_r^2 sigma_az^2 sigma_el^2]);                 %Covariance Matrix spherical coordinates
R = [sind(elevation)*cosd(az), r*cosd(elevation)*cosd(az), -r*cosd(elevation)*sind(az);
     sind(elevation)*sind(az), r*cosd(elevation)*sind(az), -r*sind(elevation)*cosd(az);
     cosd(elevation), -r*sin(elevation), 0];                        %Spherical to cartesian rotation matrix
P_cart = R*P_sphere*R';                                             %LVLH measurement Covariance matrix

%Coordinate rotation from spherical to LVLH
x = r*sind(elevation)*sind(az);
y = r*sind(elevation)*sind(az);
z = r*cosd(elevation);
x_rth = [x y z];

%LVLH to ECI 
Debris_X = rth2car(x_rth,SC_X);

%% Debris Covariance and State Propogation
% THis point onwards is validated with MC_Simulation from Github. 
X0 = [3.9583*10^4, -1.4667*10^4, 0.1035*10^4, 1.0583, 2.8815 0.0842];       %Initial state - cartesian ECI

% Initial Cartesian Covariance Matrix (P0)
P0 = [ 24287918.0736715,  5665435.69165342,  894967.886653782, -260.261998968652,  1843.15218222622,  25.0611461380351;
       5665435.69165342,  257826685.618099,  4584696.24234150, -19879.9603291687,  247.507838264477, -643.710040075805;
       894967.886653782,  4584696.24234150,  150514.823858886, -347.533001229772,  64.0785106860803, -7.14358443006258;
      -260.261998968652, -19879.9603291687, -347.533001229772,  1.53503734807762, -0.00580176692199,  0.04990688410132;
       1843.15218222622,  247.507838264477,  64.0785106860803, -0.00580176692199,  0.14029757209040,  0.00226834102064;
       25.0611461380351, -643.710040075805, -7.14358443006258,  0.04990688410132,  0.00226834102064,  0.00244767655339]*10^-3; % Initial covariance matrix of state vector distribution
% Position (km) and velocity (km/s) uncertainties
tspan = [0 10*24*3600];                                %Time propagation scale
tf = tspan(2);
t = tspan(1);
initial_STM = eye(6);
% Covariance propogation 
% PK+1 = Phi*PK*Phi';

initial_conditions = [X0, reshape(initial_STM,1,36)];
initial_conditions_cov = [X0, reshape(P0,1,36)];

%% Simulator Loop 

P_flat = reshape(P0,[],1);

options_ode45 = odeset('AbsTol',1e-6, 'RelTol',1e-9);
%State Propogation + STM prop
[Time_out,X_out] = ode45(@(t,x)eom_2BP_with_STM(t,x,mu),tspan,initial_conditions,options_ode45);
%Covariance Prop + State Prop
[Time_out,Cov_out] = ode45(@(t,x)cov_2BP_with_STM(t,x,mu),tspan,initial_conditions_cov,options_ode45);

sigma_x = sqrt(Cov_out(:,7));
sigma_y = sqrt(Cov_out(:,14));
sigma_z = sqrt(Cov_out(:,21));

final_cov_flat = Cov_out(end, 7:42);
final_cov = reshape(final_cov_flat, 6, 6);
disp('Final Covariance Matrix:');
disp(final_cov);

%% Plot Covariance Growth Over Time
Time_out_hrs = Time_out/3600; 
figure;
plot(Time_out_hrs, 3*sigma_x, 'r', 'LineWidth', 1.5); 
hold on;
plot(Time_out_hrs, 3*sigma_y, 'g', 'LineWidth', 1.5);
plot(Time_out_hrs, 3*sigma_z, 'b', 'LineWidth', 1.5);
xlabel('Time [hrs]');
ylabel('Position Uncertainty (km)');
legend('\sigma_x', '\sigma_y', '\sigma_z');
title('Covariance Growth Over Time');
grid on;
xlim([0 50])

%% Plot Position vs. Time with ±3σ Uncertainty Envelopes    
figure;

% Plot x-coordinate
subplot(3,1,1);
plot(Time_out_hrs, Cov_out(:,1), 'b', 'LineWidth', 1.5);
hold on;
plot(Time_out_hrs, Cov_out(:,1) + 3*sigma_x, 'r--', 'LineWidth', 1.0);
plot(Time_out_hrs, Cov_out(:,1) - 3*sigma_x, 'r--', 'LineWidth', 1.0);
xlabel('Time (hrs)');
ylabel('X Position (km)');
title('X Position with ±3σ Uncertainty');
grid on;
ylim([-1e5 1e5])
xlim([0 50])
% Plot y-coordinate
subplot(3,1,2);
plot(Time_out_hrs, Cov_out(:,2), 'b', 'LineWidth', 1.5);
hold on;
plot(Time_out_hrs, Cov_out(:,2) + 3*sigma_y, 'r--', 'LineWidth', 1.0);
plot(Time_out_hrs, Cov_out(:,2) - 3*sigma_y, 'r--', 'LineWidth', 1.0);
xlabel('Time (hrs)');
ylabel('Y Position (km)');
title('Y Position with ±3σ Uncertainty');
grid on;
ylim([-1e5 1e5])
xlim([0 50])
% Plot z-coordinate
subplot(3,1,3);
plot(Time_out_hrs, Cov_out(:,3), 'b', 'LineWidth', 1.5);
hold on;
plot(Time_out_hrs, Cov_out(:,3) + 3*sigma_z, 'r--', 'LineWidth', 1.0);
plot(Time_out_hrs, Cov_out(:,3) - 3*sigma_z, 'r--', 'LineWidth', 1.0);
xlabel('Time (hrs)');
ylabel('Z Position (km)');
title('Z Position with ±3σ Uncertainty');
grid on;

ylim([-1e5 1e5])
xlim([0 50])


function [dx_dt] = eom_2BP_with_STM(t,x,mu)
    
    dr_dt_and_dv_dt = [x(4:6);
                       -mu/(norm(x(1:3))^3) * x(1:3)];

    rx = x(1);
    ry = x(2);
    rz = x(3);

    r0_norm = norm(x(1:3));
    coef1 = 3*mu/(2*r0_norm^5);
    coef2 = mu/(2*r0_norm^3);

    A_21 = [coef1*rx^2 - coef2, coef1*rx*ry, coef1*rx*rz;
            coef1*rx*ry, coef1*ry^2 - coef2, coef1*ry*rz;
            coef1*rx*rz, coef1*ry*rz, coef1*rz^2-coef1];

    A = [zeros(3), eye(3); A_21, zeros(3)];

    STM = reshape(x(7:end), 6,6);
    dSTM_dt = A*STM;
    dx_dt = [dr_dt_and_dv_dt; reshape(dSTM_dt,36,1)];

end

function dXdt = cov_2BP_with_STM(t, X, mu)
    % X is assumed to contain [state (6); covariance (flattened 36 elements)]
    % Extract the state and covariance matrix
    state = X(1:6);
    P = reshape(X(7:42), 6, 6);
    
    % Compute the state derivative (position and velocity)
    r = state(1:3);
    rx = r(1);
    ry = r(2);
    rz = r(3);
    v = state(4:6);
    r_norm = norm(r);
    drdt = v;
    dvdt = -mu * r / r_norm^3;
    
    coef1 = 3*mu/(2*r_norm^5);
    coef2 = mu/(2*r_norm^3);
    
    % Build the system matrix A (Jacobian)
    A_21 = [coef1*rx^2 - coef2, coef1*rx*ry, coef1*rx*rz;
            coef1*rx*ry, coef1*ry^2 - coef2, coef1*ry*rz;
            coef1*rx*rz, coef1*ry*rz, coef1*rz^2-coef1];

    A = [zeros(3), eye(3); A_21, zeros(3)];
    
    % Compute the derivative of the covariance matrix
    dPdt = A * P + P * A';
    
    % Flatten dPdt to a vector
    dPdt_flat = reshape(dPdt, 36, 1);
    
    % Construct the derivative of the complete state
    dXdt = [drdt; dvdt; dPdt_flat];
end

function x_car = rth2car(x_rth,s_car)

% rth2car.m - Vector reference frame transformation.
%   Radial-trasversal-h to Cartesian reference frame.
%
% PROTOTYPE:
%	x_car = rth2car(x_rth, s_car)
%
% DESCRIPTION:
%   Transformation from radial-trasversal-h to Cartesian reference frame.
%   Radial-transversal-h (rth) reference frame: {r,t,h}
%       r-axis: direction of the orbit radius.
%       h-axis: direction of the angular momentum.
%       t-axis: in the orbit plane, completes the reference frame. If the
%               orbit is circular it is in the direction of the velocity
%               vector.
%   Cartesian (car) reference frame: {x,y,z}
%       inertial reference frame.
%
% INPUT:
%	x_rth[3]    Vector to be transformed, expressed in {r,t,h}.
%  	s_car[6]    State vector (position [L], velocity [L/T]) of the orbiting
%               body, expressed in {x,y,z}.
%
% OUTPUT:
%  	x_car[3,1]  Vector transformed into {x,y,z}.
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {r,t,h};
%       - we want the thrust vector in {x,y,z}.
%   In this case:
%       x_rth = Thrust vector in {r,t,h};
%       s_car = [position, velocity] of the spacecraft in {x,y,z};
%       x_car = Thrust vector, transformed in {x,y,z}.
%
% CALLED FUNCTIONS:
%   crossFast
%
% AUTHOR:
%   Camilla Colombo, 03/03/2006, MATLAB, rth2car.m
%
% PREVIOUS VERSION:
%   Camilla Colombo, 03/03/2006, MATLAB, rth_carT.m
%       - Header and function name in accordance with guidlines.
%   
% CHANGELOG:
%   10/01/2007, REVISION: Matteo Ceriotti
%   11/02/2008, Matteo Ceriotti: Help improved.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   29/03/2010, Camilla Colombo: crossFast used.
%
% -------------------------------------------------------------------------

x_rth = x_rth(:);
s = s_car(:);
r = s(1:3);
v = s(4:6);
r_ = r/norm(r);
h = crossFast(r,v);
h_ = h/norm(h);
t_ = crossFast(h_,r_);

A = [r_ t_ h_];

x_car = A*x_rth;
end

function out = kep2car(kep,mu,p)

% kep2car.m - Convertion from Keplerian orbital elements to Cartesian
%   position and velocity.
%
% PROTOTYPE:
%   out = kep2car(kep, mu, p)
%
% DESCRIPTION:
%   Converts from Keplerian orbital elements to Cartesian position and
%   velocity. All units to be consistent each other. Angles in radians.
%   Note: In the case of hyperbola, theta must be such that the point is on
%       the physical leg of the hyperbola (the leg around the attracting
%       mass).
%  
% INPUT:
%	kep[6]      Vector of keplerian elements: [a e i Om om theta], where
%               theta is the true anomaly. a in [L], angles in [rad].
%               In case of hyperbola (e>1), it must be: kep(1)<0.
%   mu          Planetary gravity constant [L^3/(M*T^2)].
%   p           Semi-latus rectum [L]. Only used for parabola case.
%
% OUTPUT:
% 	out[1,6]    State vector in cartesian coordinates (position [L],
%               velocity [L/T]).
%
% CALLED FUCTIONS:
%   (none)
%
% REFERENCES:
%	ASCL, "ToolboxManual", cap 2, 2010 - for setting on threshold on the
%       eccentricity for considering the orbit to be parabolic.
%
% AUTHOR:
%   Massimiliano Vasile, 2002, MATLAB, kep2car.m
%
% PREVIOUS VERSION:
%   Massimiliano Vasile, 2002, MATLAB, kep2cart.m
%       - Header and function name in accordance with guidlines.
%
% CHANGELOG:
%   12/02/2007, REVISION: Matteo Ceriotti
%   17/04/2007, Camilla Colombo, Matteo Ceriotti: checked case of the
%   	non-physical leg of the hyperbola and added note.
%   27/02/2008, Matteo Ceriotti, Daniel Novak: added the parabolic case, in
%     	which kep(1)=rp.
%   30/09/2009, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%   12/11/2009, Matteo Ceriotti, Jeannette Heiligers: changed the criterion
%       to fall into the parabolic case from e==1 to e>=0.99999999999999999
%       to prevent singularities with highly eccentric orbits. Added the
%       semi-latus rectum to the input to eliminate the assumption
%       kep(1)=rp.
%   07/10/2010, Matteo Ceriotti:
%       - Added condition on parabolic case for eccentricity to be smaller
%         than (1+elimitpar). This corrects a bug that made function not
%         working on all hyperbolic orbits.
%       - Removed "if nargin < 3" that made function use input p value even
%         if not parabola.
%
% -------------------------------------------------------------------------

% Threshold on the eccentricity for considering the orbit to be parabolic
%   Reference to report ToolboxTest cap 1.

%elimitpar = 0.99999999999999999;
% ---> CL: 07/10/2010, Matteo Ceriotti: changed way of expressing elimitpar
elimitpar = 1e-17;

e   = kep(2);
i   = kep(3);
Om  = kep(4);
om  = kep(5);
tho = kep(6);

% Rotation matrix
R(1,1) = cos(om)*cos(Om)-sin(om)*cos(i)*sin(Om);
R(2,1) = cos(om)*sin(Om)+sin(om)*cos(i)*cos(Om);
R(3,1) = sin(om)*sin(i);

R(1,2) = -sin(om)*cos(Om)-cos(om)*cos(i)*sin(Om);
R(2,2) = -sin(om)*sin(Om)+cos(om)*cos(i)*cos(Om);
R(3,2) = cos(om)*sin(i);

R(1,3) = sin(i)*sin(Om);
R(2,3) = -sin(i)*cos(Om);
R(3,3) = cos(i);

% In plane Parameters
% ---> CL: 07/10/2010, Matteo Ceriotti: added condition e <= (1+elimitpar)
if e >= (1-elimitpar) && e <= (1+elimitpar) % Parabola
    if nargin < 3
        error('Parabolic case: the semi-latus rectum needs to be provided')
    end
else
    % ---> CL: 07/10/2010, Matteo Ceriotti: Removed if nargin < 3. p shall
    % be computed even if it is given, if it is not a parabola.
    p = kep(1)*(1-e^2);     % Value of p in the input is not considered and overwritten with this one
end

r = p/(1+e*cos(tho));
xp = r*cos(tho);
yp = r*sin(tho);
wom_dot = sqrt(mu*p)/r^2;
r_dot   = sqrt(mu/p)*e*sin(tho);
vxp = r_dot*cos(tho)-r*sin(tho)*wom_dot;
vyp = r_dot*sin(tho)+r*cos(tho)*wom_dot;

% 3D cartesian vector
out(1) = R(1,1)*xp+R(1,2)*yp;
out(2) = R(2,1)*xp+R(2,2)*yp;
out(3) = R(3,1)*xp+R(3,2)*yp;

out(4) = R(1,1)*vxp+R(1,2)*vyp;
out(5) = R(2,1)*vxp+R(2,2)*vyp;
out(6) = R(3,1)*vxp+R(3,2)*vyp;

end
