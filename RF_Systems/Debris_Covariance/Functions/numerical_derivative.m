function dydx = numerical_derivative(y, dx)
    % NUMERICAL_DERIVATIVE
    % Approximate derivative of y with respect to x,
    % where the entries of y are uniformly spaced by dx.
    %
    % dydx = numerical_derivative(y, dx)
    %
    % Inputs:
    %   y  : vector of function values
    %   dx : uniform spacing between x-values
    %
    % Output:
    %   dydx : approximate derivative (same size as y)
    
    n = length(y);
    dydx = zeros(size(y));
    
    % Forward difference for the first point
    dydx(1) = (y(2) - y(1)) / dx;
    
    % Central difference for interior points
    for i = 2:n-1
        dydx(i) = (y(i+1) - y(i-1)) / (2*dx);
    end
    
    % Backward difference for the last point
    dydx(n) = (y(n) - y(n-1)) / dx;
end