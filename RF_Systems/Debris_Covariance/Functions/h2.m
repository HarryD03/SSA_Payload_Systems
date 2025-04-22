function obvs = h2(cart_state)
%Purpose: convert the cartesian estimated state into the observation model
% cart_state : [6x1] cartesian vector 

%Unpack the vector
    x = cart_state(1);
    y = cart_state(2);
    z = cart_state(3);

    vx = cart_state(4);
    vy = cart_state(5);
    vz = cart_state(6);

    rho = sqrt(x^2 + y^2 + z^2);
    az = atan2(y,x);
    el = atan2(z, sqrt(x^2 + y^2));
    
    rho_rate = ((x*vx) + (y*vy) + (z*vz))/(rho);
    az_rate = ((x*vy) - (y*vx)) / (x^2 +y^2);
    el_rate = ( (vz*(x^2+y^2)) - (z*((x*vx) + (y*vy))) ) / (sqrt(x^2 + y^2) * (x^2 + y^2 +z^2));


    obvs = [rho; az; el; rho_rate; az_rate; el_rate];

end