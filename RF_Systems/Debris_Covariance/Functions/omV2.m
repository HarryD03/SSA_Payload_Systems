function H = omV2(cart_state,meas_state)
%Purpose : Obtain the Jacobian of the measurement function constant
%refrence frame.
%cart_state : [6x1] This is the estimated state cartesian [r,v]'
%meas_state : [6x1] This is the estimated measured state in spherical [rho, az, el, rho_rate, az_rate, el_rate]

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
    az_rate = meas_state(5);
    el_rate = meas_state(6);

%Conduct partial derivatives for conversion
    rho_h = sqrt(x^2+y^2);
    %Range partial Derivatives
    drho_dx = x/rho;
    drho_dy = y/rho;
    drho_dz = z/rho;
    drho_dvx = 0;
    drho_dvy = 0;
    drho_dvz = 0;
    
    %Azimuth Partial Derivatives
    daz_dx = -y/rho_h^2;
    daz_dy = x/rho_h^2;
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

    %Az Rate Partial Derivative
    rho_h = sqrt(x^2+y^2);

    daz_rate_dx = (((y^2 - x^2)*vy) + (2*x*y*vx))/(rho_h^4);
    daz_rate_dy = ((-vx*(x^2-y^2)) + (2*x*y*vy))/(rho_h^4);
    daz_rate_dz = 0;
    daz_rate_dvx = -y/rho_h^2;
    daz_rate_dvy = x/rho_h^2;
    daz_rate_dvz = 0;
    
    %El rate partial derivative 
    R = sqrt(x^2+y^2+z^2);

    del_rate_dx = ((((2*x*vz) - (z*vx))*(rho_h*R^2)) - (vz*(x^2+y^2) - (z * ((x*vx) + (y*vy)))) * ((x*R^2)/rho_h + (2*x*rho_h))) / (rho_h^2 * (R^4));
    del_rate_dy = ((((2*y*vz) - (z*vy)) *(rho_h*R^2)) - (vz*(x^2 + y^2) - (z * ((x*vx) + (y*vy)))) * ((y*R^2)/rho_h + (2*y*rho_h))) / (rho_h^2*(R^4));
    del_rate_dz = ( ((x*vx) + (y*vy))*(z^2 - (x^2+y^2)) - (2*rho_h^2*z*vz) ) / (rho_h*R^4);
    del_rate_dvx = -(z*x)/(rho_h*R^2);
    del_rate_dvy = -(z*y)/(rho_h*R^2);
    del_rate_dvz = rho_h/(R^2);
    %Construct H matrix 

    H = [drho_dx drho_dy drho_dz drho_dvx drho_dvy drho_dvz;
         daz_dx  daz_dy  daz_dz  daz_dvx  daz_dvy  daz_dvz;
         del_dx  del_dy  del_dz  del_dvx  del_dvy  del_dvz;
         drangerate_dx drangerate_dy drangerate_dz drangerate_dvx drangerate_dvy drangerate_dvz;
         daz_rate_dx daz_rate_dy daz_rate_dz daz_rate_dvx daz_rate_dvy daz_rate_dvz;
         del_rate_dx del_rate_dy del_rate_dz del_rate_dvx del_rate_dvy del_rate_dvz];
end