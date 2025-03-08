function dXdt = cov_2BP_with_STM(t, X, mu)
    % X is assumed to contain [state (6); covariance (flattened 36 elements)]
    % Extract the state and covariance matrix
    state = X(1:6);
    P = reshape(X(7:42), 6, 6);
    
    % Compute the state derivative (position and velocity)
    r = state(1:3);
    rx = r(1);
    ry = r(2);
    rz = r(3);
    v = state(4:6);
    r_norm = norm(r);
    drdt = v;
    dvdt = -mu * r / r_norm^3;
    
    coef1 = 3*mu/(2*r_norm^5);
    coef2 = mu/(2*r_norm^3);
    
    % Build the system matrix A (Jacobian)
    A_21 = [coef1*rx^2 - coef2, coef1*rx*ry, coef1*rx*rz;
            coef1*rx*ry, coef1*ry^2 - coef2, coef1*ry*rz;
            coef1*rx*rz, coef1*ry*rz, coef1*rz^2-coef1];

    A = [zeros(3), eye(3); A_21, zeros(3)];
    
    % Compute the derivative of the covariance matrix
    dPdt = A * P + P * A';
    
    % Flatten dPdt to a vector
    dPdt_flat = reshape(dPdt, 36, 1);
    
    % Construct the derivative of the complete state
    dXdt = [drdt; dvdt; dPdt_flat];
end

%Verified through the Covariance_Prop function and MC