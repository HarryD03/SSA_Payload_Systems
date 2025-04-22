function [sph_state] =CartLVLH2SphLVLH(cart)
%Purpose: Transform cartesian LVLH to the Spherical LVLH 
%Cart : [6x1] The cartesian state vector [L]
%Sph : [6x1] The spherical state vector [rad]
%Extract states 
    x = cart(1); y = cart(2); z = cart(3);
    vx = cart(4); vy = cart(5); vz = cart(6);

    %Get the spherical position
    [azimuth, elevation, range] = cart2sph(x,y,z); %Use MATLAB provided function
    
    %Get the spherical velocity
    %Differentiate the range, az, el functions wrt time
    range_rate = (x*vx + y*vy + z*vz)/range;
    azimuth_rate = ((x*vy) - (y*vx)) / (x^2 + y^2);
    elevation_rate = ((x^2 + y^2)*vz - (z*((x*vx) + (y*vy))))/(sqrt(x^2 + y^2) * (x^2 + y^2 +z^2));

    sph_state = [range, azimuth, elevation, range_rate, azimuth_rate, elevation_rate];
end