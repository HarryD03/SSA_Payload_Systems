%% StripSAR Sizing: Pareto Front with Varying Power (Legend by p_peak)
% Goal is to maximise range and minimise uncertanity 
%TODO: 
% Add coherent pulse integration: sarazgain function - to get processing
% gain and changes antenna dimensions
clear; close all; clc

%% Design Variables
antenna_width_vec = linspace(0.1, 10, 10);     % Antenna widths [m] (assume SAR antenna is 0.5 cm thick)
range_vec         = linspace(50e3, 200e3, 10);   % Range values [m]
res_along_vec     = 10e-2;                     % Azimuth resolution [m] (fixed, e.g., 0.2 m)

% Bandwidth vector (MHz converted to Hz)
bandwidth_vec = [1, 10, 50, 80, 100, 200]*1e6;  % [Hz]

% Power vector [W]
p_peak_vec = [10, 20, 50, 70, 80, 100];  
Lp = length(p_peak_vec);

%% Fixed system parameters (other than power, which is now varied)
f = 5.5e9;  % Frequency [Hz]

%% Constants and Other Parameters
c           = 3e8;              % Speed of light [m/s]
v_rel       = 7.5e3;            % Relative velocity [m/s]
boltz_const = 1.380649e-23;     % Boltzmann constant [J/K]
T_sys       = 300;              % System noise temperature [K]
receiver_noise_db = 2;          % [dB]
receiver_noise    = 10^(receiver_noise_db/10);
L_sys_db    = 8;                % System losses [dB]
L_sys       = 10^(L_sys_db/10);
eff_trans   = 1;                % Transmitter efficiency (assumed)
lambda      = c/f;              % Wavelength [m]

% Quantization bits for data rate calculation
quantisation_bits = 2;

%% Desired Range Resolution and Pulse Compression Calculation
desired_range_resolution = 0.1;         % 10 cm resolution
B_eff_required = c/(2*desired_range_resolution);  % Required effective bandwidth [Hz] (1.5e9 Hz)
T_min = 1/B_eff_required;               % Minimum pulse width for 10 cm resolution (~6.67e-10 s)

%% Create Design Grid Dimensions
N  = length(antenna_width_vec);  % Antenna width index
M  = length(res_along_vec);      % Azimuth resolution index (only one value here)
Lr = length(range_vec);          % Range index
K  = length(bandwidth_vec);      % Bandwidth index
% Lp already defined

%% Preallocate Arrays (5-D arrays: [antenna_width, range, res_along, bandwidth, power])
SW_all         = zeros(N, Lr, M, K, Lp);
NESZ_all       = zeros(N, Lr, M, K, Lp);
data_rate_all  = nan(N, Lr, M, K, Lp);  % Data rate
mass_all       = nan(N, Lr, M, K, Lp);  % Mass (if using a formula)
SNR_dB_all     = nan(N, Lr, M, K, Lp);  % SNR in dB
CR_all         = nan(N, Lr, M, K, Lp);  % Pulse Compression Ratio

% Also preallocate storage for physical antenna area, bandwidth used, and peak power
A_phys_all   = zeros(N, M);              % Depends only on antenna width and resolution
bw_used_all  = zeros(N, Lr, M, K, Lp);
p_peak_store = zeros(N, Lr, M, K, Lp);     % To store the transmit power for each design point

%% Compute Antenna Area (does not depend on range, bandwidth, or power)
for i = 1:N
    for m = 1:M
        D_AT = 2 * res_along_vec(m);  % Azimuth dimension [m]
        A_phys_all(i,m) = antenna_width_vec(i) * D_AT;  % [m^2]
    end
end

%% Main Loop: Loop Over Power and Bandwidth Values as well
for l = 1:Lp
    p_peak_current = p_peak_vec(l);
    for k = 1:K
        bw = bandwidth_vec(k);
        t_pulse = 1/bw;           % Transmitted pulse duration [s] for this bandwidth
        Res_range = c*t_pulse/2;    % Uncompressed Range resolution [m]
        
        % Calculate the pulse compression ratio required:
        % CR = transmitted pulse width / minimum pulse width needed for 10 cm resolution
        CR_required = t_pulse / T_min;  
        
        for i = 1:N
            for m = 1:M
                % Compute average resolution for this (res_along, bw)
                avgRes = (Res_range + res_along_vec(m)) / 2;
                
                % Store current p_peak for all range values at (i,m,k,l)
                for j = 1:Lr
                    p_peak_store(i,j,m,k,l) = p_peak_current;
                end
                
                % We no longer mark as infeasible based on avgRes; instead, we compute the required CR.
                for j = 1:Lr
                    R_val = range_vec(j);
                    
                    % Compute Swath Width
                    SW = lambda * R_val / antenna_width_vec(i);
                    SW_all(i,j,m,k,l) = SW;
                    
                    % Compute PRF and swath time
                    PRF_max = 1 / (2*t_pulse + (2*SW)/c);
                    t_swath = 2*SW/c;
                    
                    bw_used_all(i,j,m,k,l) = bw;
                    
                    D_AT = 2 * res_along_vec(m);
                    PRF_min_local = 2*v_rel / D_AT;
                    
                    if PRF_max < PRF_min_local
                        NESZ_all(i,j,m,k,l) = NaN;
                        data_rate_all(i,j,m,k,l) = NaN;
                        mass_all(i,j,m,k,l) = NaN;
                        SNR_dB_all(i,j,m,k,l) = NaN;
                        CR_all(i,j,m,k,l) = NaN;
                    else
                        duty  = t_pulse * PRF_max;
                        P_avg = duty * p_peak_current;
                        
                        [azGain] = sarazgain(R_val,lambda,v_rel,res_along_vec,PRF_max);

                        NESZ_all(i,j,m,k,l) = 10*log10( ...
                            (2*v_rel * (4*pi*R_val)^3 * boltz_const * T_sys * receiver_noise * L_sys) ...
                            / (P_avg * azGain*(eff_trans * 4*pi * 0.7 * A_phys_all(i,m)/(lambda^2))^2 * lambda^3 * Res_range) );
                        
                        SNR_dB_all(i,j,m,k,l) = 10*log10(pi*(10e-2)^2) - NESZ_all(i,j,m,k,l);
                        
                        data_rate_all(i,j,m,k,l) = 2 * bw * quantisation_bits * t_swath * PRF_max * duty;
                        
                        mass_all(i,j,m,k,l) = 20.5044 * antenna_width_vec(i) * (2*res_along_vec(m));
                        
                        % Store the pulse compression ratio for this design point
                        CR_all(i,j,m,k,l) = CR_required;
                    end
                end
            end
        end
    end
end

%% Flatten 5D Arrays to 1D Vectors
numDesigns = N * Lr * M * K * Lp;
SW_vec       = reshape(SW_all,       numDesigns, 1);
SNR_vec      = reshape(SNR_dB_all,    numDesigns, 1);
dataRate_vec = reshape(data_rate_all, numDesigns, 1);
mass_vec     = reshape(mass_all,     numDesigns, 1);
Range_vec    = reshape(repmat(reshape(range_vec, [1,Lr]), [N,1,M,K,Lp]), numDesigns, 1);
bw_vec_store = reshape(bw_used_all,  numDesigns, 1);
p_peak_vec_flat = reshape(p_peak_store, numDesigns, 1);  % Flattened peak power
CR_vec       = reshape(CR_all, numDesigns, 1);  % Flattened compression ratio

antennaWidth_vecStore = reshape(repmat(reshape(antenna_width_vec, [N,1]), [1,Lr,M,K,Lp]), numDesigns, 1);
resAlong_vecStore     = reshape(repmat(reshape(res_along_vec, [1,M]), [N,Lr,1,K,Lp]), numDesigns, 1);
antennaLength_vecStore = 2 * resAlong_vecStore;

%% Compute SNR in Linear and the Uncertainty (sigma)
SNR_linear = 10.^(SNR_vec/10);
sigma = c ./ (2 .* bw_vec_store .* sqrt(2 .* SNR_linear));

%% Define Feasible Designs
feasibleIdx = ~isnan(SNR_vec) & (SNR_vec > 5);

%% Global Pareto Front Calculation
% Objectives:
%   Objective 1: maximize Range -> minimize -Range
%   Objective 2: minimize sigma
obj1_all = -Range_vec(feasibleIdx);
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

%% Plotting: Subplots for All p_peak Levels (Range vs. Uncertainty)
feasible_range = Range_vec(feasibleIdx);
feasible_sigma = sigma(feasibleIdx);
feasible_pPeak = p_peak_vec_flat(feasibleIdx);
feasible_SW = SW_vec(feasibleIdx);
unique_power = unique(feasible_pPeak);
nLevels = length(unique_power);
colors = lines(nLevels);
nrows = ceil(sqrt(nLevels));
ncols = ceil(nLevels / nrows);

figure;
for ip = 1:nLevels
    idx = (feasible_pPeak == unique_power(ip));
    group_range = feasible_range(idx);
    group_sigma = feasible_sigma(idx);
    
    % Compute group-specific Pareto front
    group_obj1 = -group_range;
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
    scatter(group_range/1000, group_sigma, 20, colors(ip,:), 'filled');
    hold on;
    scatter(group_range(groupParetoIdx)/1000, group_sigma(groupParetoIdx), 30, 'r', 'o');
    xlabel('Range [km]');
    ylabel('\sigma (Uncertainty) [m]');
    title(sprintf('p_{peak} = %d W', unique_power(ip)));
    grid on;
    hold off;
end

%% --- Compute Union of Group-Specific Pareto Front Indices ---
% For each unique p_peak level, compute its own Pareto front points (in terms of range and uncertainty),
% then combine all these indices into a single array.
groupParetoIdx_all = [];
for ip = 1:nLevels
    % Find indices for the current peak power group
    group_idx = find(feasible_pPeak == unique_power(ip));
    group_range = feasible_range(group_idx);
    group_sigma = feasible_sigma(group_idx);
    
    % Define objectives: maximize range (by using -range) and minimize uncertainty
    group_obj1 = -group_range;
    group_obj2 = group_sigma;
    numGroup = length(group_obj1);
    isParetoGroup = true(numGroup,1);
    
    % Compute the group-specific Pareto front indices
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
    % Map back to the overall feasible design indices and append
    groupParetoIdx = group_idx(isParetoGroup);
    groupParetoIdx_all = [groupParetoIdx_all; groupParetoIdx];
end

%% --- Figure 2: Additional Plots ---
figure;
% Subplot 1: Uncertainty vs. Swath Width
subplot(1,2,1);
hold on;
h_all = zeros(nLevels,1);
for ip = 1:nLevels
    idx = (feasible_pPeak == unique_power(ip));
    % Plot all feasible points for the current power level
    h_all(ip) = scatter(feasible_sigma(idx), feasible_SW(idx)/1000, 20, colors(ip,:), 'filled');
end
% Overlay union of all group-specific Pareto front points (red circles)
h_pareto = scatter(feasible_sigma(groupParetoIdx_all), feasible_SW(groupParetoIdx_all)/1000, 50, 'r', 'o');
legend_entries = cell(nLevels+1,1);
for ip = 1:nLevels
    legend_entries{ip} = sprintf('p_{peak} = %d W', unique_power(ip));
end
legend_entries{end} = 'Pareto Front';
legend([h_all; h_pareto], legend_entries, 'Location', 'Best');
xlabel('Uncertainty (\sigma) [m]');
ylabel('Swath Width [km]');
title('Uncertainty vs. Swath Width');
grid on;
hold off;

% Subplot 2: Range vs. Swath Width (existing code remains unchanged)
subplot(1,2,2);
hold on;
h_all = zeros(nLevels,1);
for ip = 1:nLevels
    idx = (feasible_pPeak == unique_power(ip));
    h_all(ip) = scatter(feasible_range(idx)/1000, feasible_SW(idx)/1000, 20, colors(ip,:), 'filled');
    hold on
end
h_pareto = scatter(feasible_range(groupParetoIdx_all)/1000, feasible_SW(groupParetoIdx_all)/1000, 50, 'r', 'o');
legend_entries = cell(nLevels+1,1);
for ip = 1:nLevels
    legend_entries{ip} = sprintf('p_{peak} = %d W', unique_power(ip));
end
legend_entries{end} = 'Pareto Front';
legend([h_all; h_pareto], legend_entries, 'Location', 'Best');
xlabel('Range [km]');
ylabel('Swath Width [km]');
title('Range vs. Swath Width');
grid on;
hold off;

%% Plot: Group-Specific Pareto Points for Each p_peak Level (Peak Power vs. Uncertainty)
figure;
hold on;
unique_power = unique(feasible_pPeak);  % Unique power levels from the feasible designs
colors = lines(length(unique_power));     % Generate distinct colors for each group

% Preallocate arrays to store the extreme Pareto uncertainty values for each power group.
maxSigma = zeros(length(unique_power), 1);
minSigma = zeros(length(unique_power), 1);

for ip = 1:length(unique_power)
    % Get indices for the current power level.
    group_idx = (feasible_pPeak == unique_power(ip));
    
    % Extract range and uncertainty values for this group.
    group_range = feasible_range(group_idx);   % [m]
    group_sigma = feasible_sigma(group_idx);     % [m]
    
    % Compute group-specific Pareto front.
    % Objectives:
    %   obj1: maximize range (we use -range to convert to minimization)
    %   obj2: minimize uncertainty (sigma)
    group_obj1 = -group_range;
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
    groupSigmaPareto = group_sigma(isParetoGroup);
    
    % Store the extreme values: maximum uncertainty (corresponding to max range)
    % and minimum uncertainty (corresponding to min range).
    maxSigma(ip) = max(groupSigmaPareto);
    minSigma(ip) = min(groupSigmaPareto);
    
    % All points in the group have the same peak power value.
    x_values = unique_power(ip) * ones(size(groupSigmaPareto));
    
    % Plot these Pareto points using a distinct color and marker.
    scatter(x_values, groupSigmaPareto, 20, 'filled', 'MarkerFaceColor', colors(ip,:));
end

% Overlay lines connecting the extreme Pareto points across all power levels.
plot(unique_power, maxSigma, ':', 'LineWidth', 1, 'DisplayName', 'Max Range (High \sigma)');
plot(unique_power, minSigma, ':', 'LineWidth', 1, 'DisplayName', 'Min Range (Low \sigma)');

hold off;
xlabel('Peak Power (W)');
ylabel('Uncertainty (\sigma) (m)');
title('Pareto Points across All p_{peak} Levels (Peak Power vs. Uncertainty)');
legend('Location', 'best');
grid on;

%% User Query: Find Design Data for a Given Range and Uncertainty
target_range = 100e3;    % e.g., 66.67 km [m]
target_sigma = 5;       % e.g., 0.43 m

range_error = abs(Range_vec(feasibleIdx) - target_range);
sigma_error = abs(sigma(feasibleIdx) - target_sigma);
total_error = range_error + sigma_error;

[~, idx_target_feasible] = min(total_error);
feasible_full_idx = find(feasibleIdx);
idx_target = feasible_full_idx(idx_target_feasible);

found_range       = Range_vec(idx_target);
found_sw          = SW_vec(idx_target);
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

%% Display the User-Queried Results
fprintf('\n=== USER-QUERIED DESIGN DATA ===\n');
fprintf('Target Range             = %.2f km\n', target_range/1000);
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
