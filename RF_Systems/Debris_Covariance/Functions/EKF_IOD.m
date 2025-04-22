%% EKF Test Function for IOD
% Assume pure EKF code

function [X_final, P_hist] = EKF_IOD(X_0, P_0,Z,t,PRI)
%Purpose: Reduce the Uncertanity of State through EKF so can revisit later

%X_0 : [6x1] Initial position and Velocity in the LVLH Frame [m]
%P_0 : [6x6] The intial covariance matrix for the group of measurements LVLH [m]
%Z   : [4xN] Debris State measurements [range,az,el,range_rate]
%t   : [1xN] timestamps of the measurement epochs [s]

options = odeset('AbsTol',1e-6, 'RelTol',1e-9);
dt = mean(diff(t));     % Find the average difference between 
Q = 0*diag([dt^3 dt^3 dt^3 dt dt dt]);           %Process Noise, How different are the dynamics to the 2BP?
P = zeros(6,6,length(t));
%Initialise state and covariance with inputs
X = X_0;                % initial state
P(:,:,1) = P_0;         % covariance Matrix [6x6]
N = length(t);          % Number of effective measurements post herrick gibbs
dt = mean(diff(t));     % Find the average difference between 

R = diag([46.1e-3, deg2rad(0.1), deg2rad(0.1), 0.005e-3]);                        % Measurement errors
t = dt;                                     % Measurement weights
%Kalman Loop
for i = 1:N
    % Prediction step 
    tspan = 0:PRI:t;
    [~, X_out] = ode45(@two_body, tspan, X, options);   %Predict state to next measurement  
    X_pred = X_out(end,:)';

    %Propagate Covriance Matrix
    phi0 = eye(6);                          % initial STM is the identity matrix
    phi0_flat = [reshape(phi0, 36, 1)'];    % flatten the matrix for ODE integration
    F = sm(X_pred(1:3));
    [~,phi_out] = ode45(@STM, tspan, phi0_flat, options,F); %predict covariance for next timestep
    phi = reshape(phi_out(end,:), 6, 6);  % reshape back to 6x6

    P_pred = phi*P(:,:,i)*phi' + Q;    %Covariance prediction

    %Measurement Update Step
    %Convert Measurement into Cartesian through measurement Matrix H then
    %find the residual between the predicted step and the measurements
    Z_pred = h(X_pred);
    %estimate the predicted obsverations

    H = om(X_pred,Z_pred);             %Assume the measurements Z are in spherical LVLH
    y = Z(i,:) - Z_pred';               %residual

    
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