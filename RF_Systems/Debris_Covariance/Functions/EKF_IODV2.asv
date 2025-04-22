%% EKF Test Function for IOD
% Assume pure EKF code

function [X_final, P_hist] = EKF_IODV2(X_0, P_0,Z, t, PRI,a)
%Purpose: Reduce the Uncertanity of State through EKF so can revisit later

%X_0 : [6x1] Initial position and Velocity in the LVLH Frame [m]
%P_0 : [6x6] The intial covariance matrix for the group of measurements LVLH [m]
%Z   : [4xN] Debris State measurements [range,az,el,range_rate]
%t   : [1xN] timestamps of the measurement epochs [s]

options = odeset('AbsTol',1e-6, 'RelTol',1e-9);
dt = PRI;     % Find the average difference between 
Q = 0*eye(6,6);           %Process Noise, How different are the dynamics to the 2BP?
P = zeros(6,6,length(t));
%Initialise state and covariance with inputs
X = X_0';               % initial state
P(:,:,1) = P_0;         % covariance Matrix [6x6]
N = length(t);          % Number of effective measurements post herrick gibbs

R = diag([46.1e-3, deg2rad(0.0155), deg2rad(0.0155), 0.005e-4, sqrt(2*(deg2rad(0.0155)^2)), sqrt(2*(deg2rad(0.0155)^2))]);                        % Measurement errors
t = dt;                                     % Measurement weights
%Kalman Loop


for i = 1:N
    % Prediction step 
    A = sm_LVLH(a);            %Get State dynamics Matrix
    phi = expm(A*dt);
    X_pred = phi*X;                    %Predict state
    P_pred = phi*P(:,:,i)*phi' + Q;    %Covariance prediction

    %Measurement Update Step
    %Convert Measurement into Cartesian through measurement Matrix H then
    %find the residual between the predicted step and the measurements
    
    Z_pred = h2(X_pred);                    %Find the Spherical Measurements
    %estimate the predicted obsverations

    H = omV2(X_pred,Z_pred);             %Assume the measurements Z are in spherical LVLH
    y = Z(i,:) - Z_pred';                   %residual measurement error

    
    %Compute Kalman Gain
    K = P_pred * H' / (H * P_pred * H' + R);

    %Update Step with residual
    X = X_pred + K*y';
    X_hist(i,:) = X';                   %Cartesian LVLH

    P(:,:,i+1) = (eye(6) - K *H)* P_pred;
    t = t + dt;
    
    %FOV gating
    if Z(i,2) > deg2rad(11) || Z(i,2) < deg2rad(-11)
        disp('Outside Fov...')
        break
    end
    if Z(i,3) > deg2rad(15) || Z(i,3) < deg2rad(-15)
        disp('Outside Fov...')
        break
    end
    if Z(i,1) > 50
        disp('Outside Fov...')
        break
    end
end

%Return Final state and covariance after Kalman filter
%Output should be in cartesian LVLH.

%Eliminate the terms outside FOV

X_final = X_hist;            %Position at the end of field of view
P_hist = P(:,:,1:length(X_hist));             %Kalman filter history
end