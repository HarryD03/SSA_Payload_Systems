function P_ECI = LVLH2ECI_cov(P_LVLH, r_A, v_A)
% ROTATECOV_LVLH2ECI Transforms a 6x6 covariance matrix from the LVLH frame 
%                   into the ECI frame.
%
%   P_ECI = rotateCov_LVLH2ECI(P_LVLH, r_A, v_A)
%
%   Inputs:
%       P_LVLH - (6×6) covariance matrix in LVLH coordinates (position and velocity)
%       r_A    - (3×1) reference position in ECI [km]
%       v_A    - (3×1) reference velocity in ECI [km/s]
%
%   Outputs:
%       P_ECI  - (6×6) covariance matrix in ECI coordinates
%
%   The LVLH frame is defined by:
%       x̂_LVLH = r_A / norm(r_A)
%       ẑ_LVLH = (r_A × v_A) / norm(r_A × v_A)
%       ŷ_LVLH = cross(ẑ_LVLH, x̂_LVLH)
%
%   The transformation from LVLH to ECI for a relative state [dr_LVLH; dv_LVLH]
%   is:
%
%       dr_ECI = R_l2i * dr_LVLH
%       dv_ECI = R_l2i * dv_LVLH + cross(Ω_A, dr_ECI)
%
%   where:
%       R_l2i = [x̂_LVLH, ŷ_LVLH, ẑ_LVLH]  (3×3 rotation matrix)
%       Ω_A   = (r_A × v_A) / norm(r_A)^2     (angular velocity of LVLH in ECI)
%
%   The Jacobian of this transformation with respect to [dr_LVLH; dv_LVLH]
%   is given by:
%
%       T = [ R_l2i,              0;
%             S(Ω_A)*R_l2i,       R_l2i ]
%
%   where S(Ω_A) is the 3×3 skew-symmetric matrix of Ω_A.
%
%   The covariance transformation is then:
%       P_ECI = T * P_LVLH * T'

    % Ensure inputs are column vectors.
    r_A = r_A(:);
    v_A = v_A(:);
    
    % 1) Build LVLH unit vectors in ECI coordinates
    % x̂_LVLH: along r_A
    Lx = r_A / norm(r_A);
    
    % ẑ_LVLH: orbit-normal direction
    Lz = cross(r_A, v_A);
    Lz = Lz / norm(Lz);
    
    % ŷ_LVLH: completes the right-handed coordinate system
    Ly = cross(Lz, Lx);
    
    % 2) Build the rotation matrix from LVLH to ECI.
    % The columns of R_l2i are the LVLH unit vectors expressed in ECI.
    R_l2i = [Lx, Ly, Lz];
    
    % 3) Compute the LVLH frame's angular velocity in ECI.
    Omega_A = cross(r_A, v_A) / (norm(r_A)^2);
    
    % 4) Construct the skew-symmetric matrix S(Omega_A).
    S_Omega = [   0       -Omega_A(3)  Omega_A(2);
                Omega_A(3)     0       -Omega_A(1);
               -Omega_A(2) Omega_A(1)      0      ];
    
    % 5) Construct the full 6x6 transformation matrix T.
    % This matrix maps a state in LVLH coordinates [dr_LVLH; dv_LVLH] into ECI.
    T = [R_l2i,            zeros(3,3);
         S_Omega*R_l2i,    R_l2i];
     
    % 6) Transform the covariance matrix.
    P_ECI = T * P_LVLH * T';
end