function [cart_state,cart_error] = SphLVLH2CartLVLH(sph_state,sph_error)
%Purpose: to transform the spherical LVLH to the Cartesian LVLH.
%Sph_state : [6x1]
%Sph_error : [1x6]

%unpack variables
    rho = sph_state(1); az = sph_state(2); el = sph_state(3);
    rho_rate = sph_state(4); az_rate = sph_state(5); el_rate = sph_state(6);

%get sph_error as a 6x6
    sph_error = diag(sph_error);
    %Get Cartesian position 
    [x,y,z] = sph2cart(az,el,rho);

    %Get Cartesian velocity
    vx = (rho_rate*cos(az)*cos(el)) - (rho*az_rate*sin(az)*cos(el)) - (rho*cos(az)*el_rate*sin(el));
    vy = rho_rate*sin(az)*cos(el) + rho*((az_rate * cos(az)*cos(el)) - (sin(az)*el_rate*sin(el)));
    vz = rho_rate*sin(el) + rho*el_rate*cos(el);

    cart_state = [x,y,z,vx,vy,vz];

    %Covariance Conversion
    J_11 = [cos(el)*cos(az), -rho*cos(el), -rho*sin(el)*cos(az);
            cos(el)*sin(az), rho*cos(el)*cos(az), -rho*sin(el)*sin(az);
            sin(el), 0, rho*cos(el)];
    
    J_21 = [(-sin(el)*cos(az)*el_rate)-(cos(el)*sin(az)*(az_rate)), (-rho_rate*cos(el)*sin(az)) + (rho*sin(el)*sin(az)*el_rate) - (rho*(cos(el)*cos(az)*az_rate)), (-rho_rate*sin(el)*cos(az))-(rho*cos(el)*cos(az)*el_rate) + (rho*sin(el)*sin(az)*az_rate);
            (-sin(el)*sin(az)*el_rate) + (cos(el)*cos(az)*az_rate), (rho_rate*cos(el)*cos(az)) - (rho*sin(el)*cos(az)*el_rate) - (rho*cos(el)*sin(az)*az_rate), (-rho_rate*sin(el)*sin(az)) - (rho*cos(el)*sin(az)*el_rate) - (rho*sin(el)*cos(az)*az_rate);
            (cos(el)), 0, (rho_rate*cos(el)) - (rho*sin(el)*el_rate)];
    
    J_22 = [cos(el)*cos(az), -rho*cos(el)*sin(az), -rho*sin(el)*cos(az);
            cos(el)*sin(az), rho*cos(el)*cos(az), -rho*sin(el)*sin(az);
            sin(el), 0, rho*cos(el)];

    J = [J_11 zeros(3,3);
         J_21 J_22];

    cart_error = J* sph_error *J';
end

