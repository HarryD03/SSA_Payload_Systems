function [r_B] = rotate_LVLH2ECI(dr_LVLH, dv_LVLH, r_A, v_A)
%ROTATE_LVLH2ECI  Inverts rotate_ECI2LVLH, converting from LVLH coords back to ECI.
%
% Inputs:
%   dr_LVLH : relative position in LVLH   (3x1)  if input [km], output [km]
%   dv_LVLH : relative velocity in LVLH   (3x1) [km]
%   r_A, v_A: "reference" position & velocity in ECI (3x1 each).
%             These define the LVLH axes at that instant.
%             This is the spacecraft position
%
% Outputs:
%   r_B, v_B: final position & velocity in ECI coordinates (3x1 each).
%
% The axes are the same as in rotate_ECI2LVLH:
%   L_x = r_A / |r_A| 
%   L_z = cross(r_A, v_A) / |r_A x v_A| 
%   L_y = L_x x L_z
%   R_i = [L_x L_y L_z]^T 
% so that:
%   dr_LVLH = R_i * (r_B - r_A),  and
%   dv_LVLH = [v_B - v_A] - cross(Om_A, (r_B - r_A)),
% with Om_A = cross(r_A, v_A) / |r_A|^2.
%

    % 1) Build the local LVLH unit axes (same as in rotate_ECI2LVLH)
    L_x = r_A / norm(r_A);
    L_z = cross(r_A, v_A);
    L_z = L_z / norm(L_z);
    L_y = cross(L_x, L_z);

    % 2) Rotation matrix from ECI to LVLH was R_i = [L_x L_y L_z]^T
    %    => from LVLH to ECI is R_i' = [L_x L_y L_z].
    R_i = [L_x', L_y', L_z'];   % ECI→LVLH, so its transpose is LVLH→ECI

    % 3) Reconstruct the ECI relative position from dr_LVLH
    %    dr_LVLH = R_i * dr_i => dr_i = R_i' * dr_LVLH
    dr_i = R_i * dr_LVLH;    

    % 4) Similarly for velocity:
    %    dv_LVLH = dv_i - cross(Om_A, dr_i)
    %    => dv_i = dv_LVLH + cross(Om_A, dr_i)
    Om_A = cross(r_A, v_A) / (norm(r_A)^2);
    % dv_i = dv_LVLH + cross(Om_A, dr_i);

    % 5) Now convert dv_i from LVLH axes back to ECI orientation:
    %    Actually, note that dv_LVLH is an LVLH *inertial* measure, 
    %    but the direction is the same as if it were in R_i coords. 
    %    We must also rotate dv_LVLH by R_i' if it was in r_A-based orientation.
    %    However, from the original rotate_ECI2LVLH code, you see dv_LVLH is 
    %    effectively "in the same basis" as dr_LVLH. 
    %
    %    If we used exactly the same definitions as rotate_ECI2LVLH, 
    %    then the final step is:
    % dv_i = R_i' * dv_i;   % rotate the local inertial vector back to ECI frame

    % 6) Finally, add back the reference to get the absolute ECI state
    r_B = r_A + dr_i';
    % v_B = v_A + dv_i;

end
