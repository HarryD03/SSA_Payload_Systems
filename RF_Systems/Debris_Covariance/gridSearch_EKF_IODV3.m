function [bestQ, bestResults] = gridSearch_EKF_IODV3(X_0, P_0, Z, t, PRI, a, Y_hist)
% gridSearch_EKF_IODV3: Performs a grid search over candidate process noise
% matrices Q to find the best one for the EKF.
%
% Inputs:
%   X_0   : [6x1] Initial state in LVLH coordinates.
%   P_0   : [6x6] Initial covariance matrix.
%   Z     : [4xN] Measurement data.
%   t     : [1xN] Timestamps.
%   PRI   : Pulse Repetition Interval (time step).
%   a     : Parameter for the dynamics (used in sm_LVLH).
%   Y_hist: [Nx6] "True" state history for error evaluation.
%
% Outputs:
%   bestQ     : The process noise matrix Q that minimizes the cost function.
%   bestResults: Structure containing the EKF outputs and error metrics for the best Q.
    
    % Define scaling factors for the candidate Q matrices
    scales = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1, 2, 5, 10];
    bestCost = inf;
    bestQ = [];
    
    % Loop over candidate scaling factors for the position and velocity noise
    for alpha = scales
        for beta = scales
            % Construct candidate Q matrix:
            Q_candidate = diag([alpha * 1e-7, alpha * 1e-7, alpha * 1e-7, ...
                                beta * 1e-10, beta * 1e-10, beta * 1e-10]);
            
            % Run the EKF with the candidate Q
            [X_est, P_hist] = EKF_IODV3(X_0, P_0, Z, t, PRI, a, Q_candidate);
            
            % Compute the RMSE per state dimension
            errorVec = X_est - Y_hist(1:size(X_est,1),:);  % Align sizes if filter stops early
            RMSE = sqrt(mean(errorVec.^2, 1));  % RMSE for each state dimension
            totalRMSE = sum(RMSE);              % Sum RMSE (or use max(RMSE) if desired)
            
            % Compute the covariance metric using the final covariance matrix
            finalP = P_hist(:,:,end);
            covMetric = trace(finalP);
            
            % Combined cost function: weights can be adjusted if desired
            cost = 0.5 * totalRMSE + 0.5 * covMetric;
            
            % Save the candidate if it has the lowest cost so far
            if cost < bestCost
                bestCost = cost;
                bestQ = Q_candidate;
                bestResults.X_est = X_est;
                bestResults.P_hist = P_hist;
                bestResults.RMSE = RMSE;
                bestResults.covMetric = covMetric;
                bestResults.alpha = alpha;
                bestResults.beta = beta;
            end
        end
    end
    
    % Display the best candidate's information
    fprintf('Best Q matrix found:\n');
    disp(bestQ);
    fprintf('Scaling factors: alpha = %.2f, beta = %.2f\n', bestResults.alpha, bestResults.beta);
    fprintf('RMSE per state dimension: %s\n', mat2str(bestResults.RMSE));
    fprintf('Final covariance trace: %e\n', bestResults.covMetric);
end
