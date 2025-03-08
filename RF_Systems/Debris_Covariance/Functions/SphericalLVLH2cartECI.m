function [Debris_state, P_cart] = SphericalLVLH2cartECI_RTH(state_spherical, error_spherical, s_car, velFlag)

%Inputs : 
%State_spherical : [3x1]/[6x1] The position vector or the state vector in
%                 spherical coordinates [km]
%Error_spherical : [3x3/[6x6] The Observation errors associated with the
%                  debris
%s_car           : [6x1] The Cartesian ECI vector of the observation
%                  platform [km]
%VelFlag         : Boolean value. True means that the debris measurements
%                  are full or partial state observations

% Converts (r,az,el) in the R–T–H frame into ECI,
% and propagates the 3x3 position covariance.

    % Unpack spherical coords
    r      = state_spherical(1);
    az_deg = state_spherical(2);
    el_deg = state_spherical(3);

    sigma_r    = error_spherical(1);
    sigma_az   = error_spherical(2);  % deg
    sigma_el   = error_spherical(3);  % deg

    az_rad = deg2rad(az_deg);
    el_rad = deg2rad(el_deg);

    % Create 3x3 spherical covariance (position only)
    P_sph_pos = diag([sigma_r^2, sigma_az^2, sigma_el^2]);

    % 1) Build x_rth
    x_r = r*cos(el_rad)*cos(az_rad);
    x_t = r*cos(el_rad)*sin(az_rad);
    x_h = r*sin(el_rad);
    x_rth = [x_r; x_t; x_h];    % in km

    % 2) LVLH to ECI
    r_A = s_car(1:3);           %[km] The position of the Radar Satellite
    v_A = s_car(4:6);           %[km] The velocity of the Radar Satellite
    dr_LVLH = x_rth;
    dv_LVLH = 0;                %[kms-1] This depends on the doppler processing
    [Debris_pos_ECI] = rotate_LVLH2ECI(dr_LVLH, dv_LVLH, r_A, v_A);        

    % 3) Build the partial-derivative matrix from (r,az_deg,el_deg)->(x_rth)
    c = pi/180;
    dxr_dr   = cos(el_rad)*cos(az_rad);
    dxr_daz  = -r*cos(el_rad)*sin(az_rad)*c;
    dxr_del  = -r*sin(el_rad)*cos(az_rad)*c;

    dxt_dr   = cos(el_rad)*sin(az_rad);
    dxt_daz  =  r*cos(el_rad)*cos(az_rad)*c;
    dxt_del  = -r*sin(el_rad)*sin(az_rad)*c;

    dxh_dr   = sin(el_rad);
    dxh_daz  = 0;
    dxh_del  = r*cos(el_rad)*c;

    J_rth_sph = [dxr_dr, dxr_daz, dxr_del;
                 dxt_dr, dxt_daz, dxt_del;
                 dxh_dr, dxh_daz, dxh_del];

    % 4) The RTH->ECI rotation matrix from s_car
    rVec = s_car(1:3);
    vVec = s_car(4:6);
    r_  = rVec / norm(rVec);
    h_  = cross(rVec, vVec);   
    h_ = h_/norm(h_);
    t_  = cross(h_, r_);
    A   = [r_', t_', h_'];        % 3x3

    % 5) Full Jacobian from spherical -> ECI
    J_ECI_sph = A * J_rth_sph;

    % 6) Map the covariance
    P_cart_pos = J_ECI_sph * P_sph_pos * J_ECI_sph.';

    % If velocity is included
    if velFlag == 1
        % ... do velocity partials similarly and produce a 6x6 ...
        % For now, just store zeros or do separate steps.
        P_cart = zeros(6);
        P_cart(1:3, 1:3) = P_cart_pos;
        % Also do the actual velocity transform ...
        drdt    = state_spherical(4);
        dazdt_d = state_spherical(5);
        deldt_d = state_spherical(6);
        ...
        % form the velocity in rth, then call rth2car, etc.

        % Combine final state
        Debris_state = [Debris_pos_ECI;  SomeVelocity_ECI];
    else
        % Only position measured
        Debris_state = Debris_pos_ECI;
        P_cart       = P_cart_pos;  
    end

end
