%% Two Body Problem
function [dstatedt] = two_body(t,state)
%Purpose: Define the equations of motion for the 2BP orbit.

%t : used so be put into ODE45
%state : [6x1] Position and velocity of object in orbit

mue = 3.986e5;
% Component
x = state(1);
y = state(2);
z = state(3);
% Vector
r = [x y z]';
% Magnitude
R = norm(r);


% Component
dx = state(4);
dy = state(5);
dz = state(6);
% Vector
v = [dx dy dz]';
% Magnitude
V = norm(v);


% Component
ddx = -mue*x/R^3;
ddy = -mue*y/R^3;
ddz = -mue*z/R^3;
% Vector
a = [ddx ddy ddz]';

dstatedt = [v;a]; % Return state vector for next step
end

