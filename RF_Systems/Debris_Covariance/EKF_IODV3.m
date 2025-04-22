function [X_final, P_hist] = EKF_IODV3(X_0, P_0, Z, t, PRI, a, Q)
% EKF_IODV3: Extended Kalman Filter for IOD with process noise Q.
%
% Inputs:
%   X_0 : [6x1] Initial state (position and velocity in LVLH) [m]
%   P_0 : [6x6] Initial covariance matrix
%   Z   : [4xN] Measurements [range, az, el, range_rate]
%   t   : [1xN] Timestamps [s]
%   PRI : Pulse Repetition Interval (time step) [s]
%   a   : Parameter required for the dynamics (used in sm_LVLH)
%   Q   : [6x6] Process noise covariance matrix
%
% Outputs:
%   X_final : [Nx6] Filtered state history (Cartesian LVLH coordinates)
%   P_hist  : [6x6xN] Covariance history

    % Set ODE options (if needed for dynamics integration)
    options = odeset('AbsTol',1e-6, 'RelTol',1e-9);
    
    % Initialise variables
    dt = mean(diff(t));  % Average time step from timestamps
    N = length(t);       % Number of measurement epochs
    P = zeros(6,6,N+1);  % Pre-allocate covariance history
    X = X_0';            % Initial state as row vector (1x6)
    P(:,:,1) = P_0;      % Initial covariance matrix
    X_hist = zeros(N,6); % Pre-allocate state history
    
    % Measurement noise covariance matrix (assumed constant)
    R = diag([46.1e-3, deg2rad(0.1), deg2rad(0.1), 0.005e-3, sqrt(2*(deg2rad(0.1)^2)), sqrt(2*(deg2rad(0.1)^2))]);
    
    % Kalman Filter Loop
    for i = 1:N
        % Prediction Step
        A = sm_LVLH(a);          % Get state dynamics matrix in LVLH frame
        phi = expm(A * dt);      % State transition matrix (using matrix exponential)
        X_pred = phi * X;       % Predicted state (column vector)
        P_pred = phi * P(:,:,i) * phi' + Q;  % Covariance prediction
        
        % Measurement Update Step
        Z_pred = h(X_pred);      % Compute predicted measurement (convert state to measurement space)
        H = omV2(X_pred, Z_pred);  % Measurement matrix (assumes spherical LVLH measurements)
        y = Z(i,:)' - Z_pred;    % Innovation (measurement residual)
        
        % Compute Kalman Gain
        K = P_pred * H' / (H * P_pred * H' + R);
        
        % Update state and covariance using the residual
        X_updated = X_pred + K * y;
        P_updated = (eye(6) - K * H) * P_pred;
        
        % Save the updated state estimate (transpose to row vector for history)
        X_hist(i,:) = X_updated';
        
        % Save updated covariance for the next iteration
        P(:,:,i+1) = P_updated;
        
        % Field-Of-View (FOV) Gating
        if Z(i,2) > deg2rad(11) || Z(i,2) < deg2rad(-11)
            disp('Outside FOV (azimuth)...');
            break;
        end
        if Z(i,3) > deg2rad(15) || Z(i,3) < deg2rad(-15)
            disp('Outside FOV (elevation)...');
            break;
        end
        if Z(i,1) > 50
            disp('Outside FOV (range)...');
            break;
        end
        
        % Set state for next iteration
        X = X_updated;
    end

    % Return the state and covariance history (only up to valid estimates)
    X_final = X_hist(1:i,:);
    P_hist = P(:,:,1:i+1);
end

