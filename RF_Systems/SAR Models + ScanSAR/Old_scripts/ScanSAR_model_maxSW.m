%%% Modified Code for scanSAR Mode
clear; clc; close all;

%% Choose operational mode: 'stripSAR' or 'scanSAR'
operational_mode = 'scanSAR';
if strcmpi(operational_mode, 'scanSAR')
    Nsub = 12;  % number of subswaths in scanSAR mode (must be >= 1)
else
    Nsub = 1;
end

%% Design Variables
antenna_width_vec = linspace(0.01, 10, 20);     % Antenna widths [m]
range_vec         = linspace(50e3, 200e3, 30);    % Slant range values [m]
res_along_vec     = 10e-2;                      % Azimuth resolution [m]

% Bandwidth vector (MHz converted to Hz)
bandwidth_vec = [0.1, 0.5, 1, 10, 50, 80]*1e6;   % [Hz]

% Peak Power vector [W]
p_peak_vec = [0.1, 1, 5, 10, 20, 35, 50];  % watts
Lp = length(p_peak_vec);

%% Fixed system parameters (other than power, which is now varied)
f = 5e9;  % Frequency [Hz]

%% Constants and Other Parameters
c           = 7e8;              % Speed of light [m/s]
v_rel       = 3e3;              % Relative velocity [m/s]
boltz_const = 1.380649e-23;     % Boltzmann constant [J/K]
T_sys       = 300;              % System noise temperature [K]
receiver_noise_db = 2;          % [dB]
receiver_noise    = 10^(receiver_noise_db/10);
L_sys_db    = 5;                % System losses [dB]
L_sys       = 10^(L_sys_db/10);
eff_trans   = 1;                % Transmitter efficiency (assumed)
lambda      = c/f;              % Wavelength [m]

% Quantization bits for data rate calculation
quantisation_bits = 2;

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
        Res_range = c*t_pulse/2;    % Uncompressed range resolution [m]
        
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
                    % In scanSAR mode, total swath width is expanded over Nsub subswaths:
                    if Nsub > 1
                        SW = Nsub * SW_strip;
                    else
                        SW = SW_strip;
                    end
                    SW_all(i,j,m,k,l) = SW;
                    
                    graz_ang = lambda/2 * antenna_width_vec(i);
                    % For PRF bounds in scanSAR, use the instantaneous (per subswath) swath:
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
                        prf = PRF_max;
                        % In scanSAR, each subswath receives pulses only every Nsub-th pulse:
                        if Nsub > 1
                            duty = t_pulse * (prf / Nsub);
                        else
                            duty = t_pulse * prf;
                        end
                        P_avg = duty * p_peak_current;
                        
                        % Compute azimuth gain using a helper function.
                        % For scanSAR, reduce the gain by sqrt(Nsub) due to shorter integration:
                        azGain = sarazgain(R_val, lambda, v_rel, res_along_vec, prf);
                        if Nsub > 1
                            azGain = azGain / sqrt(Nsub);
                        end
                        FilterGain = 1;
                        
                        % Compute NESZ (Noise Equivalent Sigma Zero)
                        NESZ_all(i,j,m,k,l) = 10*log10(...
                            (2*v_rel * (4*pi*R_val)^3 * boltz_const * T_sys * receiver_noise * L_sys) ...
                            / (P_avg * FilterGain * azGain * (eff_trans * 4*pi * 0.7 * A_phys_all(i,m)/(lambda^2))^2 * lambda^3 * Res_range) );
                        
                        % Compute SNR (in dB) from NESZ.
                        SNR_dB_all(i,j,m,k,l) = 10*log10(pi*(10e-2)^2) - NESZ_all(i,j,m,k,l);
                        % Further reduce SNR in scanSAR mode due to reduced dwell time:
                        if Nsub > 1
                            SNR_dB_all(i,j,m,k,l) = SNR_dB_all(i,j,m,k,l) - 10*log10(Nsub);
                        end
                        
                        % Data rate calculation:
                        data_rate_all(i,j,m,k,l) = 2 * bw * quantisation_bits * t_swath * prf * duty;
                        
                        % Mass estimation (using an example empirical formula)
                        mass_all(i,j,m,k,l) = 20.5044 * antenna_width_vec(i) * (2*res_along_vec(m));
                        
                        % Store the pulse compression ratio for this design point
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
beamwidth = lambda./(antenna_width_vec);  % Synthetic aperture beamwidth
sigma_ang = beamwidth./sqrt(2.*SNR_linear);
sigma = (sigma_r + sigma_ang)/2;          % average uncertainty

FOV = asind((SW_vec/2)./Range_vec);         % Half cone FOV
Range_vec = Range_vec .* cosd(FOV);


%% Define Feasible Designs
feasibleIdx = ~isnan(SNR_vec) & (SNR_vec > 5);

%% Global Pareto Front Calculation (Swath Width vs. Uncertainty)
% Objectives:
%   Objective 1: maximize Swath Width -> minimize -Swath Width
%   Objective 2: minimize uncertainty (sigma)
obj1_all = -SW_vec(feasibleIdx);
obj2_all = sigma(feasibleIdx);
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

%% Plotting: Subplots for All p_peak Levels (Swath Width vs. Uncertainty)
feasible_sigma = sigma(feasibleIdx);
feasible_SW = SW_vec(feasibleIdx);
feasible_pPeak = p_peak_vec_flat(feasibleIdx);
unique_power = unique(feasible_pPeak);
nLevels = length(unique_power);
colors = lines(nLevels);
nrows = ceil(sqrt(nLevels));
ncols = ceil(nLevels / nrows);

figure;
for ip = 1:nLevels
    idx = (feasible_pPeak == unique_power(ip));
    group_sigma = feasible_sigma(idx);
    group_SW = feasible_SW(idx);
    
    % Compute group-specific Pareto front for swath width vs. uncertainty:
    group_obj1 = -group_SW;
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
    scatter(group_sigma, group_SW/1000, 20, colors(ip,:), 'filled');
    hold on;
    scatter(group_sigma(groupParetoIdx), group_SW(groupParetoIdx)/1000, 30, 'r', 'o');
    xlabel('Uncertainty (\sigma) [m]');
    ylabel('Swath Width [km]');
    title(sprintf('p_{peak} = %d W', unique_power(ip)));
    grid on;
    hold off;
end

%% --- Compute Union of Group-Specific Pareto Front Indices ---
% (These indices refer to the feasible set for swath width vs. uncertainty)
groupParetoIdx_all = [];
for ip = 1:nLevels
    group_idx = find(feasible_pPeak == unique_power(ip));
    group_SW = feasible_SW(group_idx);
    group_sigma = feasible_sigma(group_idx);
    
    group_obj1 = -group_SW;
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
    groupParetoIdx_all = [groupParetoIdx_all; group_idx(groupParetoIdx)];
end

%% Figure 3: Subplots for Range vs. Uncertainty by Power Level
% In each subplot, all feasible designs for a given power level (color coded)
% are plotted (x-axis: Uncertainty, y-axis: Range [km]). The Pareto-optimal
% points (from the swath width vs. uncertainty analysis) for that power level are highlighted.
feasible_range = Range_vec(feasibleIdx);         % [m]
feasible_sigma = sigma(feasibleIdx);             % [m]
feasible_pPeak = p_peak_vec_flat(feasibleIdx);     % [W]

figure;
nrows = ceil(sqrt(nLevels));
ncols = ceil(nLevels/nrows);
for ip = 1:nLevels
    subplot(nrows, ncols, ip);
    % Get indices in the feasible set for the current power level
    groupIndices = find(feasible_pPeak == unique_power(ip));
    group_range = feasible_range(groupIndices);
    group_sigma = feasible_sigma(groupIndices);
    
    % Plot all feasible points for this power level
    scatter(group_sigma, group_range/1000, 20, colors(ip,:), 'filled');
    hold on;
    % Determine Pareto-optimal points (from the swath width vs. uncertainty Pareto)
    % that belong to this power level
    groupParetoForThis = intersect(groupIndices, groupParetoIdx_all);
    if ~isempty(groupParetoForThis)
        scatter(feasible_sigma(groupParetoForThis), feasible_range(groupParetoForThis)/1000, 50, 'r', 'o');
    end
    xlabel('Uncertainty (\sigma) [m]');
    ylabel('Range [km]');
    title(sprintf('p_{peak} = %d W', unique_power(ip)));
    grid on;
    ylim([0 250])
    hold off;
end

%% Figure 3: Group-Specific Pareto Points for Each p_peak Level (Peak Power vs. Swath Width)
% This plot is similar to the original "Peak Power vs. Uncertainty" plot, but here
% we replace range with swath width. For each power level, the Pareto-optimal swath width
% values (from the swath width vs. uncertainty Pareto analysis) are plotted.
figure;
hold on;
unique_power = unique(feasible_pPeak);  % Unique power levels from feasible designs
colors = lines(length(unique_power));     % Generate distinct colors for each group

% Preallocate arrays to store the extreme Pareto swath width values for each power group.
maxSW = zeros(length(unique_power), 1);
minSW = zeros(length(unique_power), 1);

for ip = 1:length(unique_power)
    % Get indices for the current power level.
    group_idx = (feasible_pPeak == unique_power(ip));
    
    % Extract swath width and uncertainty values for this group.
    group_SW = feasible_SW(group_idx);   % [m]
    group_sigma = feasible_sigma(group_idx);     % [m]
    
    % Compute group-specific Pareto front.
    % Objectives:
    %   obj1: maximize swath width (we use -swath width to convert to minimization)
    %   obj2: minimize uncertainty (sigma)
    group_obj1 = -group_SW;
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
    % Extract the Pareto front for this power level.
    groupSWPareto = group_SW(isParetoGroup);
    
    % Store the extreme swath width values (converted to km).
    maxSW(ip) = max(groupSWPareto);
    minSW(ip) = min(groupSWPareto);
    
    % All points in the group have the same peak power value.
    x_values = unique_power(ip) * ones(size(groupSWPareto));
    
    % Plot these Pareto points using a distinct color and marker.
    scatter(x_values, groupSWPareto/1000, 20, 'filled', 'MarkerFaceColor', colors(ip,:));
end

% Overlay lines connecting the extreme Pareto swath width values across power levels.
plot(unique_power, maxSW/1000, ':', 'LineWidth', 1, 'DisplayName', 'Max SW (High)');
plot(unique_power, minSW/1000, ':', 'LineWidth', 1, 'DisplayName', 'Min SW (Low)');

hold off;
xlabel('Peak Power (W)');
ylabel('Swath Width (km)');
title('Pareto Points across All p_{peak} Levels (Peak Power vs. Swath Width)');
legend('Location', 'best');
grid on;

%% Swath Width vs. Data Rate Plot
% Only use feasible designs (where SNR is valid and above the threshold)
feasible_SW_km = SW_vec(feasibleIdx) / 1000;  % Convert from m to km
feasible_dataRate = dataRate_vec(feasibleIdx);

figure;
scatter(feasible_SW_km, feasible_dataRate, 20, 'filled');
xlabel('Swath Width (km)');
ylabel('Data Rate (bps)');
title('Swath Width vs. Data Rate');
grid on;
x = FOV(feasibleIdx);
%% Maximum Swath Width Design Data
% Find the design with the maximum swath width among the feasible designs.
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
found_CR_max          = CR_vec(idx_max);  % Pulse Compression Ratio
found_sw_max          = SW_vec(idx_max);
%% Display the Maximum Swath Width Design Data
fprintf('\n=== MAXIMUM SWATH WIDTH DESIGN DATA ===\n');
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
%% User Query: Find Design Data for a Given Swath Width and Uncertainty
target_sw = 3.44e3;    % e.g., 5000 m (5 km) target swath width
target_sigma = 12.917;    % e.g., 5 m uncertainty

sw_error = abs(SW_vec(feasibleIdx) - target_sw);
sigma_error = abs(sigma(feasibleIdx) - target_sigma);
total_error = sw_error + sigma_error;

[~, idx_target_feasible] = min(total_error);
feasible_full_idx = find(feasibleIdx);
idx_target = feasible_full_idx(idx_target_feasible);

found_sw          = SW_vec(idx_target);
found_range       = Range_vec(idx_target);
found_ant_width   = antennaWidth_vecStore(idx_target);
found_res_along   = resAlong_vecStore(idx_target);
found_ant_length  = antennaLength_vecStore(idx_target);
found_SNR_dB      = SNR_vec(idx_target);
found_data_rate   = dataRate_vec(idx_target);
found_mass        = mass_vec(idx_target);
found_bw          = bw_vec_store(idx_target);
found_sigma       = sigma(idx_target);
found_p_peak      = p_peak_vec_flat(idx_target);
found_CR          = CR_vec(idx_target);  % Pulse Compression Ratio
found_FOV         = FOV(idx_target);
found_sigma_r     = sigma_r(idx_target);
found_sigma_ang   = sigma_ang(idx_target);
%% Display the User-Queried Results (Swath Width Query)
fprintf('\n=== USER-QUERIED DESIGN DATA (Swath Width Query) ===\n');
fprintf('Target Swath Width       = %.2f m\n', target_sw);
fprintf('Target Uncertainty       = %.2f m\n\n', target_sigma);
fprintf('Found Range              = %.2f km\n', found_range/1000);
fprintf('Found Swath Width        = %.2f m\n', found_sw);
fprintf('Antenna Width            = %.2f m\n', found_ant_width);
fprintf('Along-track Resolution   = %.4f m\n', found_res_along);
fprintf('Antenna Length (2*res_along) = %.4f m\n', found_ant_length);
fprintf('SNR (dB)                 = %.2f dB\n', found_SNR_dB);
fprintf('Data Rate                = %.2f bps\n', found_data_rate);
fprintf('Mass                     = %.2f kg\n', found_mass);
fprintf('Bandwidth                = %.2f MHz\n', found_bw/1e6);
fprintf('Swath Width              = %.2f m\n', found_sw);
fprintf('Uncertainty (\x03C3)        = %.2f m\n', found_sigma);
fprintf('Peak Power               = %.2f W\n', found_p_peak);
fprintf('Required Pulse Compression Ratio = %.1f\n', found_CR);
fprintf('Required Half Cone Angle = %.1f\n', found_FOV);
fprintf('Required range uncertanity = %.1f\n', found_sigma_r);
fprintf('Required angular uncertanity = %.1f\n', found_sigma_ang);