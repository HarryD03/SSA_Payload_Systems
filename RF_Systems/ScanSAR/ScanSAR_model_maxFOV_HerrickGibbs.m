%%% Modified Code for scanSAR Mode (FOV vs. Uncertainty Pareto Front), now including FOV vs. Range
clear; clc; close all;
% Todo:
% Write explaination for FOV, range and uncertanity as FoM 
%% Choose operational mode: 'stripSAR' or 'scanSAR'
operational_mode = 'scanSAR';
if strcmpi(operational_mode, 'scanSAR')
    Nsub = 10;  % number of subswaths in scanSAR mode (must be >= 1)
else
    Nsub = 1;
end

%% Design Variables
antenna_width_vec = linspace(0.01, 5, 20);     % Antenna widths [m]
range_vec         = linspace(50e3, 200e3, 30);  % Slant range values [m]
res_along_vec     = 10e-2;                      % Azimuth resolution [m]
FOV_limit = 5.5; %FOV limit for Herrick Gibbs
Power_limit = 50*0.55; %Peak Power limt
Power_target = 20; % 
% Bandwidth vector (MHz converted to Hz)
bandwidth_vec = [0.5, 1, 10, 50, 80]*1e6;   % [Hz]

% Peak Power vector [W]
p_peak_vec =  25  % watts
Lp = length(p_peak_vec);

%% Fixed system parameters (other than power, which is now varied)
f = 5e9;  % Frequency [Hz] based on ITU 

%% Constants and Other Parameters
c           = 7e8;              % Speed of light [m/s]
v_rel       = 7e3;              % Relative velocity [m/s]
boltz_const = 1.380649e-23;     % Boltzmann constant [J/K]
T_sys       = 300;              % System noise temperature [K] define (look at phased array paper) 
receiver_noise_db = 2;          % [dB] define
receiver_noise    = 10^(receiver_noise_db/10);
L_sys_db    = 5;                % System losses [dB] define
L_sys       = 10^(L_sys_db/10);
eff_trans   = 1;                % Transmitter efficiency (assumed)
lambda      = c/f;              % Wavelength [m]

% Quantization bits for data rate calculation
quantisation_bits = 2;          % Where? define

%% Desired Range Resolution and Pulse Compression Calculation
desired_range_resolution = 0.1;         % 10 cm resolution
B_eff_required = c/(2*desired_range_resolution);  % Required effective bandwidth [Hz]
T_min = 1/B_eff_required;               % Minimum pulse width for 10 cm resolution

%% Create Design Grid Dimensions
N  = length(antenna_width_vec);  % Antenna width index
M  = length(res_along_vec);      % Azimuth resolution index (only one value here)
Lr = length(range_vec);          % Range index
K  = length(bandwidth_vec);      % Bandwidth index

%% Preallocate Arrays (5-D arrays: [antenna_width, range, res_along, bandwidth, power])
SW_all         = zeros(N, Lr, M, K, Lp);
NESZ_all       = zeros(N, Lr, M, K, Lp);
data_rate_all  = nan(N, Lr, M, K, Lp);  
mass_all       = nan(N, Lr, M, K, Lp);  
SNR_dB_all     = nan(N, Lr, M, K, Lp);  
CR_all         = nan(N, Lr, M, K, Lp);  

% Also preallocate storage for physical antenna area, bandwidth used, and peak power
A_phys_all   = zeros(N, M);
bw_used_all  = zeros(N, Lr, M, K, Lp);
p_peak_store = zeros(N, Lr, M, K, Lp);

%% Compute Antenna Area (does not depend on range, bandwidth, or power)
for j = 1:Lr
    for i = 1:N
        for m = 1:M
            D_AT = 2 * res_along_vec(m);  % Azimuth dimension [m]
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
        
        % Pulse Compression Ratio required:
        CR_required = t_pulse / T_min;
        
        for i = 1:N
            for m = 1:M
                % For each range value, store current p_peak
                for j = 1:Lr
                    p_peak_store(i,j,m,k,l) = p_peak_current;
                end
                
                for j = 1:Lr
                    R_val = range_vec(j);
                    
                    % In stripmap mode, swath width:
                    SW_strip = lambda * R_val / antenna_width_vec(i);
                    % In scanSAR mode, total swath is expanded over Nsub subswaths:
                    if Nsub > 1
                        SW = Nsub * SW_strip;
                    else
                        SW = SW_strip;
                    end
                    SW_all(i,j,m,k,l) = SW;
                    
                    graz_ang = lambda/2 * antenna_width_vec(i);
                    % For PRF bounds in scanSAR, use the instantaneous (per-subswath) swath:
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
                        % In scanSAR, each subswath receives pulses only every Nsub-th pulse:
                        if Nsub > 1
                            duty = t_pulse * (prf / Nsub);
                        else
                            duty = t_pulse * prf;
                        end
                        P_avg = duty * p_peak_current;
                        
                        % Compute azimuth gain (reduced by sqrt(Nsub) for scanSAR):
                        azGain = sarazgain(R_val, lambda, v_rel, res_along_vec, prf);
                        if Nsub > 1
                            azGain = azGain / sqrt(Nsub);
                        end
                        FilterGain = t_pulse*bw;
                        % Compute NESZ (Noise Equivalent Sigma Zero)
                        NESZ_all(i,j,m,k,l) = 10*log10( ...
                            (2*v_rel * (4*pi*R_val)^3 * boltz_const * T_sys * receiver_noise * L_sys) ...
                            / (P_avg * FilterGain * azGain * ...
                            (eff_trans * 4*pi * 0.7 * A_phys_all(i,m)/(lambda^2))^2 * lambda^3 * Res_range) );
                        
                        % Compute SNR (in dB) from NESZ
                        SNR_dB_all(i,j,m,k,l) = 10*log10(pi*(10e-2)^2) - NESZ_all(i,j,m,k,l);
                        % Further reduce SNR in scanSAR mode (shorter dwell)
                        if Nsub > 1
                            SNR_dB_all(i,j,m,k,l) = SNR_dB_all(i,j,m,k,l) - 10*log10(Nsub);
                        end
                        
                        % Data rate calculation
                        data_rate_all(i,j,m,k,l) = 2 * bw * quantisation_bits * t_swath * prf * duty;
                        
                        % Mass estimation (example formula)
                        mass_all(i,j,m,k,l) = 20.5044 * antenna_width_vec(i) * (2*res_along_vec(m));
                        
                        % Pulse compression ratio
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
beamwidth = lambda./(D_AT);  % Synthetic aperture beamwidth
sigma_ang = beamwidth./sqrt(2.*SNR_linear);
sigma = (sigma_r + sigma_ang)/2;          % average uncertainty

%% Compute Field-of-View (FOV)
FOV = asind((SW_vec/2)./Range_vec);  % in degrees
Range_vec = Range_vec .* cosd(FOV);

%% Define Feasible Designs
feasibleIdx = ~isnan(SNR_vec) & (SNR_vec > 5); %define SNR limit 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global Pareto Front Calculation (FOV vs. Uncertainty)
% Objectives:
%   - Maximize FOV -> minimize -FOV
%   - Minimize uncertainty (sigma)
feasible_FOV = FOV(feasibleIdx);
feasible_sigma = sigma(feasibleIdx);
obj1_all = -feasible_FOV;  % To maximize FOV, we minimize -FOV
obj2_all = feasible_sigma;
feasiblePoints_all = [obj1_all, obj2_all];
numFeasible_all = size(feasiblePoints_all, 1);
isPareto_all = true(numFeasible_all, 1);
for a = 1:numFeasible_all
    for b = 1:numFeasible_all
        if a ~= b
            if (feasiblePoints_all(b,1) <= feasiblePoints_all(a,1)) && ...
               (feasiblePoints_all(b,2) <= feasiblePoints_all(a,2)) && ...
               ((feasiblePoints_all(b,1) < feasiblePoints_all(a,1)) || (feasiblePoints_all(b,2) < feasiblePoints_all(a,2)))
                isPareto_all(a) = false;
                break;
            end
        end
    end
end
paretoIdx_feasible = find(isPareto_all);

%% Plotting: Subplots for All p_peak Levels (FOV vs. Uncertainty)
feasible_FOV_all = feasible_FOV;
feasible_sigma_all = feasible_sigma;
feasible_pPeak = p_peak_vec_flat(feasibleIdx);
unique_power = unique(feasible_pPeak);
nLevels = length(unique_power);
colors = lines(nLevels);
nrows = ceil(sqrt(nLevels));
ncols = ceil(nLevels / nrows);

figure;
for ip = 1:nLevels
    idx = (feasible_pPeak == unique_power(ip));
    group_sigma = feasible_sigma_all(idx);
    group_FOV = feasible_FOV_all(idx);
    
    % Group-specific Pareto front for FOV vs. Uncertainty
    group_obj1 = -group_FOV;
    group_obj2 = group_sigma;
    numGroup = length(group_obj1);
    isParetoGroup = true(numGroup,1);
    for a = 1:numGroup
        for b = 1:numGroup
            if a ~= b
                if (group_obj1(b) <= group_obj1(a)) && (group_obj2(b) <= group_obj2(a)) && ...
                   ((group_obj1(b) < group_obj1(a)) || (group_obj2(b) < group_obj2(a)))
                    isParetoGroup(a) = false;
                    break;
                end
            end
        end
    end
    groupParetoIdx = find(isParetoGroup);
    
    subplot(nrows, ncols, ip);
    scatter(group_sigma, group_FOV, 20, colors(ip,:), 'filled');
    hold on;
    scatter(group_sigma(groupParetoIdx), group_FOV(groupParetoIdx), 30, 'r', 'o');
    yline(FOV_limit)
    xlabel('Uncertainty (\sigma) [m]');
    ylabel('FOV [deg]');
    title(sprintf('p_{peak} = %g W', unique_power(ip)));
    grid on;
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure: Group-Specific Pareto Points for Each p_peak Level (Peak Power vs. FOV)
figure;
hold on;
colors = lines(length(unique_power));

maxFOV = zeros(length(unique_power), 1);
minFOV = zeros(length(unique_power), 1);

for ip = 1:length(unique_power)
    group_idx = (feasible_pPeak == unique_power(ip));
    group_FOV = feasible_FOV_all(group_idx);
    group_sigma = feasible_sigma_all(group_idx);
    
    group_obj1 = -group_FOV;
    group_obj2 = group_sigma;
    numGroup = length(group_obj1);
    isParetoGroup = true(numGroup, 1);
    for a = 1:numGroup
        for b = 1:numGroup
            if a ~= b
                if (group_obj1(b) <= group_obj1(a)) && (group_obj2(b) <= group_obj2(a)) && ...
                   ((group_obj1(b) < group_obj1(a)) || (group_obj2(b) < group_obj2(a)))
                    isParetoGroup(a) = false;
                    break;
                end
            end
        end
    end
    groupPareto = group_FOV(isParetoGroup);
    
    maxFOV(ip) = max(groupPareto);
    minFOV(ip) = min(groupPareto);
    
    x_values = unique_power(ip) * ones(size(groupPareto));
    scatter(x_values, groupPareto, 20, 'filled', 'MarkerFaceColor', colors(ip,:));
end

plot(unique_power, maxFOV, ':', 'LineWidth', 1);
plot(unique_power, minFOV, ':', 'LineWidth', 1);
hold off;
xlabel('Peak Power (W)');
ylabel('FOV [deg]');
title('Pareto Points across All p_{peak} Levels (Peak Power vs. FOV)');
grid on;

%% Swath Width vs. Data Rate Plot (unchanged)
feasible_SW_km = SW_vec(feasibleIdx) / 1000;  % Convert from m to km
feasible_dataRate = dataRate_vec(feasibleIdx);

figure;
scatter(feasible_SW_km, feasible_dataRate, 20, 'filled');
xlabel('Swath Width (km)');
ylabel('Data Rate (bps)');
title('Swath Width vs. Data Rate');
grid on;

%% Maximum FOV Design Data
[max_FOV, idx_local] = max(FOV(feasibleIdx));
feasible_full_idx = find(feasibleIdx);
idx_max = feasible_full_idx(idx_local);

found_FOV_max          = FOV(idx_max);
found_range_max       = Range_vec(idx_max);
found_ant_width_max   = antennaWidth_vecStore(idx_max);
found_res_along_max   = resAlong_vecStore(idx_max);
found_ant_length_max  = antennaLength_vecStore(idx_max);
found_SNR_dB_max      = SNR_vec(idx_max);
found_data_rate_max   = dataRate_vec(idx_max);
found_mass_max        = mass_vec(idx_max);
found_bw_max          = bw_vec_store(idx_max);
found_sigma_max       = sigma(idx_max);
found_p_peak_max      = p_peak_vec_flat(idx_max);
found_CR_max          = CR_vec(idx_max);
found_sw_max          = SW_vec(idx_max);

fprintf('\n=== MAXIMUM FOV DESIGN DATA ===\n');
fprintf('Maximum FOV     = %.2f deg\n', found_FOV_max);
fprintf('Found Range              = %.2f km\n', found_range_max/1000);
fprintf('Antenna Width            = %.2f m\n', found_ant_width_max);
fprintf('Along-track Resolution   = %.4f m\n', found_res_along_max);
fprintf('Antenna Length (2*res_along) = %.4f m\n', found_ant_length_max);
fprintf('SNR (dB)                 = %.2f dB\n', found_SNR_dB_max);
fprintf('Data Rate                = %.2f bps\n', found_data_rate_max);
fprintf('Mass                     = %.2f kg\n', found_mass_max);
fprintf('Bandwidth                = %.2f MHz\n', found_bw_max/1e6);
fprintf('Uncertainty (\x03C3)         = %.2f m\n', found_sigma_max);
fprintf('Peak Power               = %.2f W\n', found_p_peak_max);
fprintf('Required Pulse Compression Ratio = %.1f\n', found_CR_max);
fprintf('Swath width                   = %.1f\n', found_sw_max);

%% User Query: Find Design Data for a Given FOV Angle and Range
% Define target design parameters:
target_FOV   = 11;     % Target FOV angle in degrees
target_range = 50;  % Target range in km

% Calculate errors between the design and the targets:
% (Note: Convert Range_vec from m to km for the error calculation)
FOV_error   = abs(FOV(feasibleIdx) - target_FOV);
range_error = abs((Range_vec(feasibleIdx)/1000) - target_range);
total_error = FOV_error/target_FOV + range_error/target_range;
[~, idx_target_feasible] = min(total_error);
feasible_full_idx = find(feasibleIdx);
idx_target = feasible_full_idx(idx_target_feasible);

% Retrieve the design parameters at the target index:
found_range      = Range_vec(idx_target);       % in meters
found_FOV        = FOV(idx_target);               % in degrees
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
found_sigma_r    = sigma_r(idx_target);
found_sigma_ang  = sigma_ang(idx_target);

% Print the design data:
fprintf('\n=== USER-QUERIED DESIGN DATA (FOV Angle and Range Query) ===\n');
fprintf('Target FOV Angle       = %.2f deg\n', target_FOV);
fprintf('Target Range           = %.2f km\n\n', target_range);
fprintf('Found Range            = %.2f km\n', found_range/1000);
fprintf('Found FOV Angle        = %.2f deg\n', found_FOV);
fprintf('Antenna Width          = %.2f m\n', found_ant_width);
fprintf('Along-track Resolution = %.4f m\n', found_res_along);
fprintf('Antenna Length (2*res_along) = %.4f m\n', found_ant_length);
fprintf('SNR (dB)               = %.2f dB\n', found_SNR_dB);
fprintf('Data Rate              = %.2f bps\n', found_data_rate);
fprintf('Mass                   = %.2f kg\n', found_mass);
fprintf('Bandwidth              = %.2f MHz\n', found_bw/1e6);
fprintf('Uncertainty (\x03C3)         = %.2f m\n', found_sigma);
fprintf('Peak Power             = %.2f W\n', found_p_peak);
fprintf('Required Pulse Compression Ratio = %.1f\n', found_CR);
fprintf('Range Uncertainty      = %.1f m\n', found_sigma_r);
fprintf('Angular Uncertainty    = %.1f m\n', found_sigma_ang);

%save design point 


%% NEW SECTION: PLOTTING FOV vs. RANGE
% We also do a Pareto analysis, since we want to maximize FOV and maximize Range.
% That means we'll minimize -FOV and minimize -Range.

feasible_FOV_all = FOV(feasibleIdx);
feasible_range_all = Range_vec(feasibleIdx)/1000;  % convert to km for plotting
feasible_pPeak_all = p_peak_vec_flat(feasibleIdx);

% Pareto front for each power level: objectives are -FOV, -Range
figure;
colors = lines(length(unique_power));
nrows = ceil(sqrt(nLevels));
ncols = ceil(nLevels / nrows);

for ip = 1:nLevels
    subplot(nrows, ncols, ip);
    groupIdx = find(feasible_pPeak_all == unique_power(ip));
    groupFOV = feasible_FOV_all(groupIdx);
    groupRange = feasible_range_all(groupIdx);
    
    % Build objectives to minimize: -FOV and -Range
    group_obj1 = -groupFOV; 
    group_obj2 = -groupRange; 
    numGroup = length(group_obj1);
    isParetoGroup = true(numGroup,1);
    
    for a = 1:numGroup
        for b = 1:numGroup
            if a ~= b
                % If design b is better or equal in both objectives, design a is not Pareto
                if (group_obj1(b) <= group_obj1(a)) && (group_obj2(b) <= group_obj2(a)) && ...
                   ((group_obj1(b) < group_obj1(a)) || (group_obj2(b) < group_obj2(a)))
                    isParetoGroup(a) = false;
                    break;
                end
            end
        end
    end
    paretoIdx_group = find(isParetoGroup);
    
    scatter(groupRange, groupFOV, 20, colors(ip,:), 'filled');
    hold on;
    scatter(groupRange(paretoIdx_group), groupFOV(paretoIdx_group), 30, 'r', 'o');
    xlabel('Range [km]');
    ylabel('FOV [deg]');
    title(sprintf('FOV vs. Range at p_{peak} = %g W', unique_power(ip)));
    grid on;
    hold off;
end

%% New Section: Uncertainty vs. Peak Power for Various FOV-Dependent Designs (p_peak >= 1W)
% Ensure the feasible arrays for uncertainty, range, FOV, and peak power are defined:
feasible_sigma_all = sigma(feasibleIdx);           % Uncertainty [m]
feasible_range_all = Range_vec(feasibleIdx);         % Range [m]
feasible_FOV_all   = FOV(feasibleIdx);               % FOV [deg]
feasible_pPeak_all = p_peak_vec_flat(feasibleIdx);   % Peak Power [W]

% Get unique peak power values from the feasible designs (only >= 1W)
unique_power_all = unique(feasible_pPeak_all);
unique_power     = unique_power_all(unique_power_all >= 1);
numLevels        = length(unique_power);

% Preallocate arrays to store the uncertainty values for each power level
uncertainty_maxFOV  = zeros(numLevels,1);  % Uncertainty for design with max FOV
uncertainty_100km   = nan(numLevels,1);    % Uncertainty for design with range closest to 100 km
uncertainty_150km   = nan(numLevels,1);    % Uncertainty for design with range closest to 150 km
uncertainty_200km   = nan(numLevels,1);    % Uncertainty for design with range closest to 200 km

% Define target ranges in meters
target_ranges = [100e3, 150e3, 200e3];

% Loop over each unique peak power level
for ip = 1:numLevels
    % Select designs with the current peak power value
    idx = (feasible_pPeak_all == unique_power(ip));
    group_sigma = feasible_sigma_all(idx);
    group_range = feasible_range_all(idx);
    group_FOV   = feasible_FOV_all(idx);
    
    % 1. For the maximum FOV design: record its uncertainty
    [~, maxIdx] = max(group_FOV);
    uncertainty_maxFOV(ip) = group_sigma(maxIdx);
    
    % 2. For the design with range closest to 100 km: record its uncertainty
    diff_100 = abs(group_range - target_ranges(1));
    [~, idx100] = min(diff_100);
    uncertainty_100km(ip) = group_sigma(idx100);
    
    % 3. For the design with range closest to 150 km: record its uncertainty
    diff_150 = abs(group_range - target_ranges(2));
    [~, idx150] = min(diff_150);
    uncertainty_150km(ip) = group_sigma(idx150);
    
    % 4. For the design with range closest to 200 km: record its uncertainty
    diff_200 = abs(group_range - target_ranges(3));
    [~, idx200] = min(diff_200);
    uncertainty_200km(ip) = group_sigma(idx200);
end

% Plot the results: Uncertainty (σ) versus Peak Power (p_peak)
figure;
hold on;
plot(unique_power, uncertainty_maxFOV, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Max FOV design');
plot(unique_power, uncertainty_100km,  'ro-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Design at 100 km');
plot(unique_power, uncertainty_150km,  'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Design at 150 km');
plot(unique_power, uncertainty_200km,  'go-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Design at 200 km');
xline(Power_limit)
hold off;
xlabel('Peak Power (W)');
ylabel('Uncertainty (\sigma) [m]');
title('Uncertainty vs. Peak Power for Various FOV-Dependent Designs (p_{peak} \geq 1W)');
legend('show');
grid on;
