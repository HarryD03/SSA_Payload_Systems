%% Herrick Gibbs
function [v2] = herrickgibbs(r1, r2, r3, t1, t2, t3,mue)
%Purpose: complete herrick gibbs IOD using 3 different orbit position
%vectors, and 3 different times associated

%[r1,r2,r3] are in geocentric co-ordinates where geocentric angle between
%them is less the 1 degree

%[t1, t2, t3] are the time associated with the position vector [seconds]

% mue is Standard gravitional constant of Earth

dt31 = (t3 - t1);
dt32 = (t3 - t2);
dt21 = (t2 - t1);
r1mag = norm(r1);
r2mag = norm(r2);
r3mag = norm(r3);

%Check if CoPlanar
Z_23 = cross(r2,r3);
ang_coplanar = 90 - acosd((dot(Z_23,r1))/(norm(Z_23)*norm(r1)));
ang_12 = acosd(dot(r1,r2)/(norm(r1)*norm(r2)));
ang_23 = acosd(dot(r2,r3)/(norm(r2)*norm(r3)));

v2 = -dt32*((1/(dt21*dt31))+(mue/(12*r1mag^3)))*r1 +...
    (dt32-dt21)*((1/(dt21*dt32))+(mue/(12*r2mag^3)))*r2 +...
    dt21*((1/(dt32*dt31))+(mue/(12*r3mag^3)))*r3; % Velocity around the 2nd point

end

%This is valid for IOD. If the norm of the acceleration between the points
%is ~0, then the points must correspond to a similar debris.

%Verified 