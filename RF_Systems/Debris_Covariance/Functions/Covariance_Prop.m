function [final_state, final_cov, t_limit, Time_out, Cov_out] = Covariance_Prop(X0, P0, tspan, error_limit)
    mu = 3.986e5;

    % Prepare the initial condition vector:
    % X0 is the state (6 elements) and P0 is the 6x6 covariance matrix.
    initial_conditions_cov = [X0, reshape(P0,1,36)];
    
    % Set ODE options.
    options_ode45 = odeset('AbsTol',1e-6, 'RelTol',1e-9);

    % Propagate the combined state and covariance.
    [Time_out, Cov_out] = ode45(@(t,x)cov_2BP_with_STM(t,x,mu), tspan, initial_conditions_cov, options_ode45);

    % Extract the standard deviations (square roots of the variances)
    sigma_x = sqrt(Cov_out(:,7));    % Covariance element for x variance
    sigma_y = sqrt(Cov_out(:,14));   % Covariance element for y variance
    sigma_z = sqrt(Cov_out(:,21));   % Covariance element for z variance
    
    % Find the time when any position uncertainty exceeds the error limit.
    for i = 1:length(sigma_x)
        if (sigma_x(i) >= error_limit) || (sigma_y(i) >= error_limit) || (sigma_z(i) >= error_limit)
            t_limit = Time_out(i);
            % Extract the covariance matrix at this time:
            final_cov_flat = Cov_out(i, 7:42);
            final_cov = reshape(final_cov_flat, 6, 6);
            % Extract the state at this time:
            final_state_flat = Cov_out(i, 1:6);
            final_state = reshape(final_state_flat, 6, 1);
            break;
        end
    end

    % Optional: The function can still generate plots if needed.
    % However, it is often better to separate plotting from computation.

    % For example, to plot covariance growth, you might include:
    Time_out_hrs = Time_out / 3600;
    figure;
    plot(Time_out_hrs, 3*sigma_x, 'r', 'LineWidth', 1.5);
    hold on;
    plot(Time_out_hrs, 3*sigma_y, 'g', 'LineWidth', 1.5);
    plot(Time_out_hrs, 3*sigma_z, 'b', 'LineWidth', 1.5);
    xlabel('Time [hrs]');
    ylabel('Position Uncertainty (km)');
    legend('\sigma_x', '\sigma_y', '\sigma_z');
    title('Covariance Growth Over Time');
    grid on;
    xlim([0 50]);
end
