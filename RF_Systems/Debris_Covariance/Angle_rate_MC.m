%% GOAL: Estimate the Error in angle rate estimates through a Monte Carlo Simulation
% Monte Carlo gives 'true values' 
clc; clear; close all
orbit = 100;                %Number of different Orbits/Debris
MeasurementNoise = 0;       %Measurment noise.
tspan = [0 90*60];          %Propogation time - one orbit
dt = 10;                    %Time difference for Simulation - pulse width for chosen design? 
FOV = 11;                   %Half cone field of view
monteCarloFiniteDiffDemo(orbit,MeasurementNoise,tspan,dt)            %Call Monte Carlo Simulation




function monteCarloFiniteDiffDemo(numOrbits, measurementNoiseDeg, tSpan, dt,FOV)
    % monteCarloFiniteDiffDemo
    % ------------------------------------------------
    % Monte Carlo to estimate the finite-difference error
    % in azimuth/elevation rates for debris orbits in
    % the 600--800 km altitude range, with inclination
    % in [80,100] degrees.  A radar "station" is fixed
    % at 600 km altitude, range limit 50 km -> TBD.
    %
    % By Harry Dawson, 25/02/2025

    % Earth/radar constants
    RE = 6378.0;                    % Earth radius [km]
    muEarth = 398600.4418;         % Earth gravitational parameter [km^3/s^2]
    radarAltitude = 600;           % radar-sat altitude [km]
    radarRangeLimit = 50;          % radar max range [km]
    
    % Radar "station" fixed in ECI at some point on x-axis for simplicity
    %    radius = Earth radius + 600 km

    % SC parameters TBC 
    h = 590;%km
    e = 0;
    i = 80;
    Om = 80;
    om = 0;
    theta = 0;

   
    radarPosECI = coe2rv(RE + h, e, deg2rad(i), deg2rad(Om), deg2rad(om), deg2rad(theta),muEarth);  % [x, y, z] km
    
    % We will collect FD errors over ALL random orbits/time steps
    azRateErrors = [];
    elRateErrors = [];
    
    % For repeatability
    rng(123);
    
    % 2) Monte Carlo over multiple random orbits
    for n = 1:numOrbits
        % Generate random orbit within constraints
        %  Altitude ~ 600--800 km => a ~ (RE+600) -- (RE+800)
        a  = randUniform(RE+600, RE+650);    % [km]
        e  = randUniform(0.0, 0.1);          % eccentricity
        iDeg = 80;         % inclination in degrees
        iRad = deg2rad(iDeg);
        RAAN = deg2rad(randUniform(0, 360));
        ArgP = deg2rad(randUniform(0, 360));
        nu   = deg2rad(randUniform(0, 360)); % true anomaly
        
        % Convert COEs to initial position/velocity in ECI
        [r0, v0] = coe2rv(a, e, iRad, RAAN, ArgP, nu, muEarth);
        
        % Propagate with two-body ODE
        [tArray, stateArray] = propagateOrbit(r0, v0, muEarth, tSpan, dt);
        % stateArray: Nx6 => [x, y, z, vx, vy, vz] each row
        
        % Convert to AZ/EL from radar's perspective
        Npts = size(stateArray,1);
        azTrue = zeros(Npts,1);
        elTrue = zeros(Npts,1);
        for k = 1:Npts                     %Cycle through the number of points to define the time when the debris within FOV
                                           %%Create IF loop which ensures
                                           %%the Debris is within the FOV
                                           %%then execute the rest.


            rDebris = stateArray(k,1:3)';  % ECI position of debris
            rVec = rDebris - radarPosECI;  % station->debris vector
            dist = norm(rVec);             % Line of Sight distance
            
            if dist <= radarRangeLimit    
                [azTrue(k), elTrue(k)] = eciToAzEl(rVec);
                if (azTrue(k) < FOV_az) && (elTrue(k) < FOV_el)
                    azTrue(k) = azTrue(k);
                    elTrue(k) = elTrue(k);
                else
                    elTrue(k) = NaN;
                    azTrue(k) = NaN;
                end
            else
                azTrue(k) = NaN;
                elTrue(k) = NaN;
            end
        end
        
        % "True" angle rates (central difference) ignoring measurement
        % noise Within the Field of View. Currently this is wrong

        azRateMeas = centralDiff(azTrue, dt);
        elRateMeas = centralDiff(elTrue, dt);
        
        azRateTrue = 
        elRateTrue =
        % Now add measurement noise (in degrees => radians)
        noiseAz = deg2rad(measurementNoiseDeg)*randn(Npts,1);
        noiseEl = deg2rad(measurementNoiseDeg)*randn(Npts,1);
        
        azMeas = azTrue + noiseAz;
        elMeas = elTrue + noiseEl;
        
        % If the object is out of range, we have NaNs.
        % We'll only do FD where we have valid data around k-1, k, k+1.
        
        azRateFD = centralDiff(azMeas, dt);
        elRateFD = centralDiff(elMeas, dt);
        
        % Compute error = (FD - True) only where both are valid
        validMask = ~isnan(azRateFD) & ~isnan(azRateTrue);
        azDiff = azRateFD(validMask) - azRateTrue(validMask);
        elDiff = elRateFD(validMask) - elRateTrue(validMask);
        
        % Accumulate differences for statistics
        azRateErrors = [azRateErrors; azDiff(:)];
        elRateErrors = [elRateErrors; elDiff(:)];
    end
    
    % 3) Final statistics
    sigmaAzRate = std(azRateErrors, 'omitnan');
    sigmaElRate = std(elRateErrors, 'omitnan');
    
    fprintf('\n--------------------------------------------------\n');
    fprintf('Monte Carlo Finite-Difference Results:\n');
    fprintf('  Number of orbits: %d\n', numOrbits);
    fprintf('  Altitude range:   600--800 km\n');
    fprintf('  Inclination range: 80--100 deg\n');
    fprintf('  Radar altitude:   600 km, range limit: 50 km\n');
    fprintf('  dt = %.1f s, total time = %.1f s\n', dt, tSpan(2));
    fprintf('  Measurement noise = %.2f deg (1-sigma)\n', measurementNoiseDeg);
    fprintf('\nResults (rad/s):\n');
    fprintf('  STD of Az rate error: %.6g rad/s\n', sigmaAzRate);
    fprintf('  STD of El rate error: %.6g rad/s\n', sigmaElRate);
    fprintf('--------------------------------------------------\n');
end

% =====================================================
% =============== Supporting Functions ================
% =====================================================

function val = randUniform(a, b)
    % Generate a uniform random number in [a, b].
    val = a + (b - a)*rand();
end

function [rECI, vECI] = coe2rv(a, e, i, RAAN, argp, nu, mu)
    % coe2rv: Convert classical orbital elements to ECI state
    % Inputs:
    %   a     - semi-major axis [km]
    %   e     - eccentricity
    %   i     - inclination [rad]
    %   RAAN  - right ascension of ascending node [rad]
    %   argp  - argument of perigee [rad]
    %   nu    - true anomaly [rad]
    %   mu    - gravitational parameter [km^3/s^2]
    %
    % Outputs:
    %   rECI, vECI in km, km/s
    
    

    % 1) Orbital parameters
    p = a*(1 - e^2);         % semi-latus rectum
    rOrbit = p/(1 + e*cos(nu));
    
    % 2) Position in "perifocal" frame
    xOrb = rOrbit*cos(nu);
    yOrb = rOrbit*sin(nu);
    zOrb = 0;
    
    % 3) Velocity in "perifocal" frame
    h = sqrt(mu*p);
    vxOrb = -mu/h * sin(nu);
    vyOrb =  mu/h * (e + cos(nu));
    vzOrb = 0;
    
    % 4) Rotation from perifocal to ECI (3-1-3) with RAAN, i, argp
    cO = cos(RAAN); sO = sin(RAAN);
    ci = cos(i);    si = sin(i);
    cw = cos(argp); sw = sin(argp);
    
    % Combined rotation matrix R = Rz(RAAN)*Rx(i)*Rz(argp)
    RzO = [ cO  -sO  0
            sO   cO  0
            0    0   1];
    Rxi = [ 1  0   0
            0  ci  -si
            0  si   ci];
    Rzw = [ cw  -sw  0
            sw   cw  0
            0    0   1];
    
    ROT = RzO*Rxi*Rzw;
    
    % 5) Rotate position, velocity
    rECI = ROT*[xOrb; yOrb; zOrb];
    vECI = ROT*[vxOrb; vyOrb; vzOrb];
end

function [tOut, stateOut] = propagateOrbit(r0, v0, mu, tSpan, dt)
    % propagateOrbit
    % Simple two-body propagation using ODE45.
    % r0, v0 in km, km/s
    % tSpan = [t0, tEnd]
    % dt = output time step
    
    y0 = [r0; v0];  % 6x1 initial state
    % time vector
    tVec = (tSpan(1):dt:tSpan(2))';
    
    opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
    sol = ode45(@(t, y) twoBodyOde(t, y, mu), tSpan, y0, opts);
    
    % Interpolate solution on uniform time grid
    stateOut = zeros(length(tVec), 6);
    for k = 1:length(tVec)
        stateOut(k,:) = deval(sol, tVec(k));
    end
    tOut = tVec;
end

function dydt = twoBodyOde(~, y, mu)
    % twoBodyOde
    % y = [x, y, z, vx, vy, vz]
    rx = y(1); ry = y(2); rz = y(3);
    vx = y(4); vy = y(5); vz = y(6);
    
    r = sqrt(rx^2 + ry^2 + rz^2);
    ax = -mu*rx/(r^3);
    ay = -mu*ry/(r^3);
    az = -mu*rz/(r^3);
    
    dydt = [vx; vy; vz; ax; ay; az];
end

function [az, el] = eciToAzEl(rVec)
    % eciToAzEl
    % Given station->debris vector (rVec) in ECI or local frame,
    % return Az, El in [radians].
    %
    % For simplicity, we interpret rVec as if the station is the origin
    % in a local-horizon sense. This is purely illustrative.
    
    % Distance in horizontal plane
    rx = rVec(1);
    ry = rVec(2);
    rz = rVec(3);
    
    rhoXY = sqrt(rx^2 + ry^2);
    el = atan2(rz, rhoXY);
    az = atan2(ry, rx);
    
    % Force az in [0, 2*pi) if you want, but not required
    % az = mod(az, 2*pi);
end

function rate = centralDiff(values, dt)
    % centralDiff
    % Simple central difference for approximate derivative:
    %   rate(k) = [values(k+1) - values(k-1)] / (2*dt)
    %
    % rate(1) and rate(end) set to NaN because we can't do a 2-sided difference.
    
    N = length(values);
    rate = NaN*ones(N,1);
    for k = 2:N-1
        if ~isnan(values(k-1)) && ~isnan(values(k+1)) && ~isnan(values(k))
            rate(k) = (values(k+1) - values(k-1)) / (2*dt);
        end
    end
end
