function [cart_state] = SphLVLH2CartLVLH(sph_state)
%Purpose: to transform the spherical LVLH to the Cartesian LVLH

%unpack variables
    rho = sph_state(1); az = sph_state(2); el = sph_state(3);
    rho_rate = sph_state(4); az_rate = sph_state(5); el_rate = sph_state(6);

    %Get Cartesian position 
    [x,y,z] = sph2cart(az,el,rho);
    
    %Get Cartesian velocity 

end

