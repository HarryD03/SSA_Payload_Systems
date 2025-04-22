%%% Modified Code for scanSAR Mode (FOV vs. Uncertainty Pareto Front), now including FOV vs. Range and Frequency Dependence
clear; clc; close all;

%% Choose operational mode: 'stripSAR' or 'scanSAR'
operational_mode = 'scanSAR';
if strcmpi(operational_mode, 'scanSAR')
    Nsub = 8;  % number of subswaths in scanSAR mode (must be >= 1)
else
    Nsub = 1;
end

%% Define simulation parameters that are independent of frequency
antenna_width_vec = linspace(0.01, 10, 50);     % Antenna widths [m]
range_vec         = linspace(40e3, 200e3, 50);  % Slant range values [m]
res_along_vec     = 10e-2;                      % Azimuth resolution [m]
FOV_limit = 5.5;       % FOV limit for Herrick Gibbs
Power_limit = 50*0.55; % Peak Power limit
Power_target = 20;     
bandwidth_vec = [1, 10, 50, 80]*1e6;   % Bandwidth vector in Hz
p_peak_vec = 25; % Peak Power vector [W]
Lp = length(p_peak_vec);

% Fixed system parameters (independent of frequency, except for wavelength)
c           = 7e8;              % Speed of light [m/s]
v_rel       = 7e3;              % Relative velocity [m/s]
boltz_const = 1.380649e-23;     % Boltzmann constant [J/K]
T_sys       = 300;              % System noise temperature [K]
receiver_noise_db = 2;          % [dB]
receiver_noise    = 10^(receiver_noise_db/10);
L_sys_db    = 5;                % System losses [dB]
L_sys       = 10^(L_sys_db/10);
quantisation_bits = 2;          % Quantisation bits for data rate calculation

% Desired range resolution and pulse compression
desired_range_resolution = 0.1;         % 10 cm resolution
B_eff_required = c/(2*desired_range_resolution);  % Effective bandwidth [Hz]
T_min = 1/B_eff_required;               % Minimum pulse width for 10 cm resolution

% Define frequency vector (Hz) and preallocate structure for design points
freq_vec = [5e9, 45.5e9, 69e9];
numFreq = length(freq_vec);
designPoints = struct('freq',[],'found_range',[],'found_FOV',[],...
    'found_ant_width',[],'found_res_along',[],'found_ant_length',[],...
    'found_SNR_dB',[],'found_data_rate',[],'found_mass',[],...
    'found_bw',[],'found_sigma',[],'found_p_peak',[],'found_CR',[],...
    'found_sw',[]);

%% Loop over frequencies
for fi = 1:numFreq
    f = freq_vec(fi);
    lambda = c/f;  % Wavelength [m]
    
    % Grid dimensions based on independent vectors
    N  = length(antenna_width_vec);  
    Lr = length(range_vec);          
    M  = 1;  % Only one azimuth resolution value (scalar)
    K  = length(bandwidth_vec);
    
    %% Preallocate Arrays (5-D arrays: [antenna_width, range, res_along, bandwidth, power])
    SW_all         = zeros(N, Lr, M, K, Lp);
    NESZ_all       = zeros(N, Lr, M, K, Lp);
    data_rate_all  = nan(N, Lr, M, K, Lp);  
    mass_all       = nan(N, Lr, M, K, Lp);  
    SNR_dB_all     = nan(N, Lr, M, K, Lp);  
    CR_all         = nan(N, Lr, M, K, Lp);  
    
    % Preallocate storage for physical antenna area, bandwidth used, and peak power
    A_phys_all   = zeros(N, M);
    bw_used_all  = zeros(N, Lr, M, K, Lp);
    p_peak_store = zeros(N, Lr, M, K, Lp);
    
    %% Compute Antenna Area (dependent on lambda via graz_ang)
    for j = 1:Lr
        for i = 1:N
            for m = 1:M
                D_AT = 2 * res_along_vec;  % Along-track dimension [m]
                A_phys_all(i,m) = antenna_width_vec(i) * D_AT;  % [m^2]
                R_val = range_vec(j);
                graz_ang = lambda/2 * antenna_width_vec(i);
                A_min = (4*v_rel*lambda*R_val)/c * tand(graz_ang);
                if A_phys_all(i,m) < A_min
                    A_phys_all(i,m) = NaN;
                end
            end
        end
    end
    
    %% Main Loop: Loop Over Power and Bandwidth Values
    for l = 1:Lp
        p_peak_current = p_peak_vec(l);
        for k = 1:K
            bw = bandwidth_vec(k);
            t_pulse = 1/bw;           % Transmitted pulse duration [s]
            Res_range = c*t_pulse/2;  % Uncompressed range resolution [m]
            CR_required = t_pulse / T_min;  % Pulse Compression Ratio required
            
            for i = 1:N
                for m = 1:M
                    for j = 1:Lr
                        p_peak_store(i,j,m,k,l) = p_peak_current;
                    end
                    
                    for j = 1:Lr
                        R_val = range_vec(j);
                        % Compute swath width: for stripmap mode then expand for scanSAR
                        SW_strip = lambda * R_val / antenna_width_vec(i);
                        if Nsub > 1
                            SW = Nsub * SW_strip;
                        else
                            SW = SW_strip;
                        end
                        SW_all(i,j,m,k,l) = SW;
                        
                        graz_ang = lambda/2 * antenna_width_vec(i);
                        % Use instantaneous swath for PRF bounds in scanSAR
                        if Nsub > 1
                            [PRF_min_local, PRF_max] = sarprfbounds(v_rel, res_along_vec, SW_strip, graz_ang);
                        else
                            [PRF_min_local, PRF_max] = sarprfbounds(v_rel, res_along_vec, SW, graz_ang);
                        end
                        t_swath = 2*SW/c;  % swath dwell time
                        bw_used_all(i,j,m,k,l) = bw;
                        
                        if PRF_max < PRF_min_local
                            NESZ_all(i,j,m,k,l) = NaN;
                            data_rate_all(i,j,m,k,l) = NaN;
                            mass_all(i,j,m,k,l) = NaN;
                            SNR_dB_all(i,j,m,k,l) = NaN;
                            CR_all(i,j,m,k,l) = NaN;
                        else
                            prf = PRF_min_local;
                            if Nsub > 1
                                duty = t_pulse * (prf / Nsub);
                            else
                                duty = t_pulse * prf;
                            end
                            P_avg = duty * p_peak_current;
                            
                            azGain = sarazgain(R_val, lambda, v_rel, res_along_vec, prf);
                            if Nsub > 1
                                azGain = azGain / sqrt(Nsub);
                            end
                            FilterGain = t_pulse*bw;
                            NESZ_all(i,j,m,k,l) = 10*log10( ...
                                (2*v_rel * (4*pi*R_val)^3 * boltz_const * T_sys * receiver_noise * L_sys) ...
                                / (P_avg * FilterGain * azGain * ...
                                (4*pi * 0.7 * A_phys_all(i,m)/(lambda^2))^2 * lambda^3 * Res_range) );
                            
                            SNR_dB_all(i,j,m,k,l) = 10*log10(pi*(10e-2)^2) - NESZ_all(i,j,m,k,l);
                            if Nsub > 1
                                SNR_dB_all(i,j,m,k,l) = SNR_dB_all(i,j,m,k,l) - 10*log10(Nsub);
                            end
                            
                            data_rate_all(i,j,m,k,l) = 2 * bw * quantisation_bits * t_swath * prf * duty;
                            
                            mass_all(i,j,m,k,l) = 20.5044 * antenna_width_vec(i) * (2*res_along_vec);
                            
                            CR_all(i,j,m,k,l) = CR_required;
                        end
                    end
                end
            end
        end
    end
    
    %% Flatten 5D Arrays to 1D Vectors (for post-processing)
    numDesigns = N * Lr * M * K * Lp;
    SW_vec       = reshape(SW_all,       numDesigns, 1);
    SNR_vec      = reshape(SNR_dB_all,    numDesigns, 1);
    dataRate_vec = reshape(data_rate_all, numDesigns, 1);
    mass_vec     = reshape(mass_all,     numDesigns, 1);
    Range_vec    = reshape(repmat(reshape(range_vec, [1,Lr]), [N,1,M,K,Lp]), numDesigns, 1);
    bw_vec_store = reshape(bw_used_all,  numDesigns, 1);
    p_peak_vec_flat = reshape(p_peak_store, numDesigns, 1);
    CR_vec       = reshape(CR_all, numDesigns, 1);
    
    antennaWidth_vecStore = reshape(repmat(reshape(antenna_width_vec, [N,1]), [1,Lr,M,K,Lp]), numDesigns, 1);
    resAlong_vecStore     = reshape(repmat(reshape(res_along_vec, [1,M]), [N,Lr,1,K,Lp]), numDesigns, 1);
    antennaLength_vecStore = 2 * resAlong_vecStore;
    
    %% Compute SNR in Linear and the Uncertainty (sigma)
    SNR_linear = 10.^(SNR_vec/10);
    sigma_r = c ./ (2 .* bw_vec_store .* sqrt(2 .* SNR_linear));
    % Use the antenna length (2*res_along) for beamwidth calculation
    beamwidth = lambda ./ (antennaLength_vecStore);
    sigma_ang = beamwidth ./ sqrt(2 .* SNR_linear);
    sigma = (sigma_r + sigma_ang) / 2;
    
    %% Compute Field-of-View (FOV) and update Range
    FOV = asind((SW_vec/2)./Range_vec);  % in degrees
    Range_vec = Range_vec .* cosd(FOV);
    
    %% Define Feasible Designs (SNR > 5 dB)
    feasibleIdx = ~isnan(SNR_vec) & (SNR_vec > 5);
    
    %% USER-QUERIED DESIGN POINT: Target FOV and Range
    target_FOV   = 11;     % Target FOV angle in degrees
    target_range = 50;     % Target range in km
    
    % Calculate errors between the design and the targets:
    FOV_error   = abs(FOV(feasibleIdx) - target_FOV);
    range_error = abs((Range_vec(feasibleIdx)/1000) - target_range);
    total_error = FOV_error/target_FOV + range_error/target_range;
    [~, idx_target_feasible] = min(total_error);
    feasible_full_idx = find(feasibleIdx);
    idx_target = feasible_full_idx(idx_target_feasible);
    
    % Retrieve the design parameters at the target index:
    found_range      = Range_vec(idx_target);       % [m]
    found_FOV        = FOV(idx_target);               % [deg]
    found_ant_width  = antennaWidth_vecStore(idx_target);
    found_res_along  = resAlong_vecStore(idx_target);
    found_ant_length = antennaLength_vecStore(idx_target);
    found_SNR_dB     = SNR_vec(idx_target);
    found_data_rate  = dataRate_vec(idx_target);
    found_mass       = mass_vec(idx_target);
    found_bw         = bw_vec_store(idx_target);
    found_sigma      = sigma(idx_target);
    found_p_peak     = p_peak_vec_flat(idx_target);
    found_CR         = CR_vec(idx_target);
    found_sw         = SW_vec(idx_target);
    
    % Display the design data for this frequency
    fprintf('\n=== DESIGN POINT AT %.2f GHz ===\n', f/1e9);
    fprintf('Found Range            = %.2f km\n', found_range/1000);
    fprintf('Found FOV Angle        = %.2f deg\n', found_FOV);
    fprintf('Antenna Width          = %.2f m\n', found_ant_width);
    fprintf('Along-track Resolution = %.4f m\n', found_res_along);
    fprintf('Antenna Length         = %.4f m\n', found_ant_length);
    fprintf('SNR (dB)               = %.2f dB\n', found_SNR_dB);
    fprintf('Data Rate              = %.2f bps\n', found_data_rate);
    fprintf('Mass                   = %.2f kg\n', found_mass);
    fprintf('Bandwidth              = %.2f MHz\n', found_bw/1e6);
    fprintf('Uncertainty (\x03C3)         = %.2f m\n', found_sigma);
    fprintf('Peak Power             = %.2f W\n', found_p_peak);
    fprintf('Required PCR           = %.1f\n', found_CR);
    fprintf('Swath Width            = %.1f m\n', found_sw);
    
    % Store design point data for later comparison
    designPoints(fi).freq = f;
    designPoints(fi).found_range = found_range;
    designPoints(fi).found_FOV = found_FOV;
    designPoints(fi).found_ant_width = found_ant_width;
    designPoints(fi).found_res_along = found_res_along;
    designPoints(fi).found_ant_length = found_ant_length;
    designPoints(fi).found_SNR_dB = found_SNR_dB;
    designPoints(fi).found_data_rate = found_data_rate;
    designPoints(fi).found_mass = found_mass;
    designPoints(fi).found_bw = found_bw;
    designPoints(fi).found_sigma = found_sigma;
    designPoints(fi).found_p_peak = found_p_peak;
    designPoints(fi).found_CR = found_CR;
    designPoints(fi).found_sw = found_sw;
    
end

%% Plot: Compare Design Points (FOV vs. Range) across Frequencies
figure;
hold on;
colors = lines(numFreq);
for fi = 1:numFreq
    scatter(designPoints(fi).found_range/1000, designPoints(fi).found_FOV, 100, colors(fi,:), 'filled');
end
xlabel('Range [km]');
ylabel('FOV [deg]');
title('Design Points for Different Frequencies');
legend(arrayfun(@(x) sprintf('%.2f GHz', x.freq/1e9), designPoints, 'UniformOutput', false));
grid on;
hold off;
