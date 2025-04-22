function [r_B, v_B] = rotate_LVLH2ECI(dr_LVLH, dv_LVLH, r_A, v_A)
%ROTATE_LVLH2ECI  
% Converts a relative position/velocity from the rotating spacecraft LVLH frame 
% to the Earth-Centered Inertial (ECI) frame.
%
% Inputs:
%   dr_LVLH : (3×1) relative position in LVLH [km]
%   dv_LVLH : (3×1) relative velocity in LVLH [km/s]
%   r_A, v_A: (3×1) reference position & velocity in ECI 
%             that define the LVLH frame.
%
% Outputs:
%   r_B, v_B: (3×1) final position & velocity in ECI.
%
% LVLH axes:
%   x̂_LVLH = r_A / norm(r_A)
%   ẑ_LVLH = (r_A × v_A) / norm(r_A × v_A)
%   ŷ_LVLH = cross(ẑ_LVLH, x̂_LVLH)
%
% The LVLH frame rotates with angular velocity:
%   Ω_A = (r_A × v_A) / norm(r_A)^2   (in ECI coordinates)
%
% Transformation:
%   dr_ECI = R * dr_LVLH
%   dv_ECI = R * dv_LVLH + cross(Ω_A, dr_ECI)
%
%   r_B = r_A + dr_ECI
%   v_B = v_A + dv_ECI
%
% where R = [x̂_LVLH, ŷ_LVLH, ẑ_LVLH] is the rotation matrix from LVLH to ECI.

    %-- 1) Build LVLH unit vectors in ECI coordinates ---------------------
    % x̂_LVLH: along r_A
    Lx = r_A / norm(r_A);

    % ẑ_LVLH: orbit-normal direction
    Lz = cross(r_A, v_A);
    Lz = Lz / norm(Lz);

    % ŷ_LVLH: completes the right-handed coordinate system
    Ly = cross(Lz, Lx);

    %-- 2) Build the rotation matrix from LVLH to ECI ---------------------
    % The columns of R_l2i are the unit vectors in ECI.
    R_l2i = [Lx, Ly, Lz];

    %-- 3) Compute the LVLH frame's angular velocity in ECI --------------
    Omega_A = cross(r_A, v_A) / (norm(r_A)^2);

    %-- 4) Transform the relative position to ECI --------------------------
    dr_ECI = R_l2i * dr_LVLH;

    %-- 5) Transform the relative velocity to ECI --------------------------
    dv_ECI = R_l2i * dv_LVLH + cross(Omega_A, dr_ECI);

    %-- 6) Compute the final ECI coordinates -------------------------------
    r_B = r_A + dr_ECI;
    v_B = v_A + dv_ECI;
end

