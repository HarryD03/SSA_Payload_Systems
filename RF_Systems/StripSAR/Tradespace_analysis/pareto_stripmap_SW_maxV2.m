%% StripSAR Sizing: Pareto Front with Varying Power (Legend by p_peak)
% Goal is to maximise range and minimise uncertainty.
% TODO:
%   
clear; clc;

%% Design Variables
antenna_width_vec = linspace(0.01, 10, 50);     % Antenna widths [m] 
range_vec         = linspace(50e3, 200e3, 20);   % Range values [m]
res_along_vec     = 10e-2;                     % Azimuth resolution [m]

% Bandwidth vector (MHz converted to Hz)
bandwidth_vec = linspace(0.01,100,20)*1e6;  % [Hz]

% Power vector [W]
p_peak_vec = [5,10, 20, 50, 70,90];  
Lp = length(p_peak_vec);

% New: Relative velocity vector [m/s]
v_rel_vec = [7.5]*1e3;  % 400, 7500, 15000 m/s           % When I change the velocity within the vector 
Lv = length(v_rel_vec);

%% Fixed system parameters (other than power and v_rel)
f = 4e9;  % Frequency [Hz]

%% Constants and Other Parameters
c           = 3e8;              % Speed of light [m/s]
boltz_const = 1.380649e-23;     % Boltzmann constant [J/K]
T_sys       = 300;              % System noise temperature [K]
receiver_noise_db = 2;          % [dB]
receiver_noise    = 10^(receiver_noise_db/10);
L_sys_db    = 8;                % System losses [dB]
L_sys       = 10^(L_sys_db/10);
eff_trans   = 1;                % Transmitter efficiency (assumed)
lambda      = c/f;              % Wavelength [m]

% Quantisation bits for data rate calculation
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

%% Preallocate Arrays (6-D arrays: [antenna_width, range, res_along, bandwidth, power, v_rel])
SW_all         = zeros(N, Lr, M, K, Lp, Lv);
NESZ_all       = zeros(N, Lr, M, K, Lp, Lv);
data_rate_all  = nan(N, Lr, M, K, Lp, Lv);
mass_all       = nan(N, Lr, M, K, Lp, Lv);
SNR_dB_all     = nan(N, Lr, M, K, Lp, Lv);
CR_all         = nan(N, Lr, M, K, Lp, Lv);

% Also preallocate storage for physical antenna area, bandwidth used, and peak power
A_phys_all   = zeros(N, M);  % Depends only on antenna width and resolution
bw_used_all  = zeros(N, Lr, M, K, Lp, Lv);
p_peak_store = zeros(N, Lr, M, K, Lp, Lv);  % Transmit power for each design point

%% Compute Antenna Area (physical; independent of range, bandwidth, etc.)
%% Compute Antenna Area (physical) with Dimensions: [antenna_width, range, res_along, v_rel]
% This array will be used to check feasibility at each range and v_rel.
A_phys_all = nan(N, Lr, M, Lv);
for v = 1:Lv
    current_v_rel = v_rel_vec(v);
    for j = 1:Lr
        R_val = range_vec(j);
        for i = 1:N
            for m = 1:M
                D_AT = 2 * res_along_vec(m);
                A_phys = antenna_width_vec(i) * D_AT;
                graz_ang = lambda/2 * antenna_width_vec(i);
                A_min = (4 * current_v_rel * lambda * R_val) / c * tand(graz_ang);
                if A_phys < A_min
                    A_phys_all(i,j,m,v) = NaN;
                else
                    A_phys_all(i,j,m,v) = A_phys;
                end
            end
        end
    end
end

%% Main Loop: Loop Over v_rel, Power, Bandwidth, etc.
for v = 1:Lv
    current_v_rel = v_rel_vec(v);
    for l = 1:Lp
        p_peak_current = p_peak_vec(l);
        for k = 1:K
            bw = bandwidth_vec(k);
            t_pulse = 1/bw;           
            Res_range = c*t_pulse/2;    
            CR_required = t_pulse / T_min;  
            for i = 1:N
                for m = 1:M
                    avgRes = (Res_range + res_along_vec(m)) / 2;
                    for j = 1:Lr
                        p_peak_store(i,j,m,k,l,v) = p_peak_current;
                        R_val = range_vec(j);
                        
                        % Compute Swath Width (independent of v_rel)
                        SW = lambda * R_val / antenna_width_vec(i);
                        SW_all(i,j,m,k,l,v) = SW;
                        PRF_max = 1 / (2*t_pulse + (2*SW)/c);
                        
                        t_swath = 2*SW/c;
                        bw_used_all(i,j,m,k,l,v) = bw;
                        D_AT = 2 * res_along_vec(m);
                        PRF_min_local = 2*current_v_rel / D_AT;
                        if PRF_max < PRF_min_local
                            NESZ_all(i,j,m,k,l,v) = NaN;
                            data_rate_all(i,j,m,k,l,v) = NaN;
                            mass_all(i,j,m,k,l,v) = NaN;
                            SNR_dB_all(i,j,m,k,l,v) = NaN;
                            CR_all(i,j,m,k,l,v) = NaN;
                        else
                            duty  = t_pulse * PRF_max;
                            P_avg = duty * p_peak_current;
                            % Call to sarazgain (assumed to return azGain)
                            [azGain] = sarazgain(R_val, lambda, current_v_rel, res_along_vec, PRF_max);
                            
                            NESZ_all(i,j,m,k,l,v) = 10*log10( ...
                                (2*current_v_rel * (4*pi*R_val)^3 * boltz_const * T_sys * receiver_noise * L_sys) ...
                                / (P_avg * azGain * (eff_trans * 4*pi * 0.7 * A_phys_all(i,m)/(lambda^2))^2 * lambda^3 * Res_range) );

                            SNR_dB_all(i,j,m,k,l,v) = 10*log10(pi*(10e-2)^2) - NESZ_all(i,j,m,k,l,v);
                            data_rate_all(i,j,m,k,l,v) = 2 * bw * quantisation_bits * t_swath * PRF_max * duty;
                            mass_all(i,j,m,k,l,v) = 20.5044 * antenna_width_vec(i) * (2*res_along_vec(m));
                            CR_all(i,j,m,k,l,v) = CR_required;
                        end
                    end
                end
            end
        end
    end
end

%% Flatten 6D Arrays to 1D Vectors
numDesigns = N * Lr * M * K * Lp * Lv;
SW_vec       = reshape(SW_all,       numDesigns, 1);
SNR_vec      = reshape(SNR_dB_all,    numDesigns, 1);
dataRate_vec = reshape(data_rate_all, numDesigns, 1);
mass_vec     = reshape(mass_all,     numDesigns, 1);
Range_vec    = reshape(repmat(reshape(range_vec, [1,Lr]), [N,1,M,K,Lp,Lv]), numDesigns, 1);
bw_vec_store = reshape(bw_used_all,  numDesigns, 1);
p_peak_vec_flat = reshape(p_peak_store, numDesigns, 1);
CR_vec       = reshape(CR_all, numDesigns, 1);

antennaWidth_vecStore = reshape(repmat(reshape(antenna_width_vec, [N,1]), [1,Lr,M,K,Lp,Lv]), numDesigns, 1);
resAlong_vecStore     = reshape(repmat(reshape(res_along_vec, [1,M]), [N,Lr,1,K,Lp,Lv]), numDesigns, 1);
antennaLength_vecStore = 2 * resAlong_vecStore;

% Also flatten the relative velocity array:
v_rel_array = repmat(reshape(v_rel_vec, [1, Lv]), [N, Lr, M, K, Lp]);
v_rel_vec_flat = reshape(v_rel_array, numDesigns, 1);

%% Compute SNR in Linear and the Uncertainty (sigma)
SNR_linear = 10.^(SNR_vec/10);
sigma = c ./ (2 .* bw_vec_store .* sqrt(2 .* SNR_linear));

%% Define Feasible Designs
feasibleIdx = ~isnan(SNR_vec) & (SNR_vec > 5);

%% Global Pareto Front Calculation (Swath Width vs. Uncertainty)
% Objectives:
%   1. Maximize Swath Width (minimize -SW)
%   2. Minimize Uncertainty (sigma)
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

%% NEW: Separate Figures for Different Relative Velocity Values
% Each figure corresponds to one relative velocity value.
% Within each figure, create subplots for each peak power value.
% Plot: x-axis = Uncertainty (sigma) [m], y-axis = Swath Width [km].
% The Pareto front in each subplot (for that group) is highlighted with red circles.
unique_v = unique(v_rel_vec_flat);
nV = length(unique_v);
for iv = 1:nV
    % Filter feasible designs for this relative velocity
    feasibleIdx_v = feasibleIdx & (v_rel_vec_flat == unique_v(iv));
    
    % Create a new figure for this velocity value
    figure;
    % Get unique peak power values for this velocity subset:
    feasible_pPeak_v = p_peak_vec_flat(feasibleIdx_v);
    unique_power_v = unique(feasible_pPeak_v);
    nLevels_v = length(unique_power_v);
    nrows_v = ceil(sqrt(nLevels_v));
    ncols_v = ceil(nLevels_v / nrows_v);
    
    for ip = 1:nLevels_v
        subplot(nrows_v, ncols_v, ip);
        % Filter indices for current peak power in this velocity subset
        group_idx = feasibleIdx_v & (p_peak_vec_flat == unique_power_v(ip));
        group_sigma = sigma(group_idx);
        group_SW = SW_vec(group_idx);
        
        % Plot all design points for this group
        scatter(group_sigma, group_SW/1000, 20, 'b', 'filled');
        hold on;
        % Compute Pareto front for this group:
        numGroup = length(group_sigma);
        if numGroup > 0
            isPareto = true(numGroup, 1);
            for a = 1:numGroup
                for b = 1:numGroup
                    if a ~= b
                        % Objectives: maximize SW (i.e., minimize -SW) and minimize sigma.
                        if ( (-group_SW(b) <= -group_SW(a)) && (group_sigma(b) <= group_sigma(a)) && ...
                             ((-group_SW(b) < -group_SW(a)) || (group_sigma(b) < group_sigma(a)) ) )
                            isPareto(a) = false;
                            break;
                        end
                    end
                end
            end
            pareto_indices = find(isPareto);
            scatter(group_sigma(pareto_indices), group_SW(pareto_indices)/1000, 50, 'r', 'o');
        end
        xlabel('Uncertainty (\sigma) [m]');
        ylabel('Swath Width [km]');
        title(sprintf('p_{peak} = %d W', unique_power_v(ip)));
        grid on;
        hold off;
    end
    sgtitle(sprintf('Swath Width vs. Uncertainty for v_{rel} = %.1f km/s', unique_v(iv)/1000));
end

%% Maximum Swath Width Design Data
[max_sw, idx_local] = max(SW_vec(feasibleIdx));
feasible_full_idx = find(feasibleIdx);
idx_max = feasible_full_idx(idx_local);

found_sw_max          = SW_vec(idx_max)/1000;
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

%% Display the Maximum Swath Width Design Data
fprintf('\n=== MAXIMUM SWATH WIDTH DESIGN DATA ===\n');
fprintf('Maximum Swath Width      = %.2f km\n', found_sw_max);
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

%% User Query: Find Design Data for a Given Swath Width and Uncertainty
target_sw = 20e3;    % e.g., 5000 m (5 km) target swath width
target_sigma = 10;    % e.g., 5 m uncertainty

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
found_CR          = CR_vec(idx_target);

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
