function H = om(cart_state,meas_state)
%Purpose : Obtain the Jacobian of the measurement function constant
%refrence frame.
%cart_state : [6x1] This is the estimated state cartesian [r,v]'
%meas_state : [4x1] This is the estimated measured state in spherical [rho,az,el,rho_rate]

%unpack vectors:
    x = cart_state(1);
    y = cart_state(2);
    z = cart_state(3);
    vx = cart_state(4);
    vy = cart_state(5);
    vz = cart_state(6);

    rho = meas_state(1);
    az = meas_state(2);
    el = meas_state(3);
    rho_rate = meas_state(4);

%Conduct partial derivatives for conversion
    
    %Range partial Derivatives
    drho_dx = x/rho;
    drho_dy = y/rho;
    drho_dz = z/rho;
    drho_dvx = 0;
    drho_dvy = 0;
    drho_dvz = 0;
    
    %Azimuth Partial Derivatives
    daz_dx = -y/rho;
    daz_dy = x/rho^2;
    daz_dz = 0;
    daz_dvx = 0;
    daz_dvy = 0;
    daz_dvz = 0;

    %Elevation Partial Derivatives
    del_dx = (z*x)/(sqrt(x^2+y^2)*rho^2);
    del_dy = -(z*y)/(sqrt(x^2+y^2)*rho^2);
    del_dz = sqrt(x^2+y^2)/(rho^2);
    del_dvx= 0;
    del_dvy = 0;
    del_dvz = 0;

    %Range Rate Partial Derivatives 
    drangerate_dx = ((vx*rho^2) - (x*((x*vx) + (y*vy) + (z*vz))))/(rho^3);
    drangerate_dy = ((vy*rho^2) - (y*((x*vx) + (y*vy) + (z*vz))))/(rho^3);
    drangerate_dz = ((vz*rho^2) - (z*((x*vx) + (y*vy) + (z*vz))))/(rho^3);
    drangerate_dvx = x/rho;
    drangerate_dvy = y/rho;
    drangerate_dvz = z/rho;

    %Construct H matrix 

    H = [drho_dx drho_dy drho_dz drho_dvx drho_dvy drho_dvz;
         daz_dx  daz_dy  daz_dz  daz_dvx  daz_dvy  daz_dvz;
         del_dx  del_dy  del_dz  del_dvx  del_dvy  del_dvz;
         drangerate_dx drangerate_dy drangerate_dz drangerate_dvx drangerate_dvy drangerate_dvz];
end