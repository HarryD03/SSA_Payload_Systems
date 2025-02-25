clear; clc; close all;

%% =========================
%  1) Define the list of subswath counts you want to analyze
%     For example, from 1 to 15 in steps of 1
%     Adjust as you wish.
% =========================
Nsub_list = 1:15;

% Pre-allocate an array to store the maximum feasible FOV for each Nsub
maxFOV_vs_Nsub = nan(size(Nsub_list));

%% =========================
%  2) Loop over each Nsub, run the main SAR design routine,
%     extract max feasible FOV for 10 W peak power, and store it.
% =========================
for idxN = 1:length(Nsub_list)
    Nsub = Nsub_list(idxN);  % number of subswaths for this iteration
    
    %% --- Design Variables ---
    antenna_width_vec = linspace(0.01, 10, 20);   % [m]
    range_vec         = linspace(50e3, 200e3, 30);% [m]
    res_along_vec     = 0.10;                     % [m]
    
    % Bandwidth vector (MHz converted to Hz)
    bandwidth_vec = [0.1, 0.5, 1, 10, 50, 80]*1e6; % [Hz]
    
    % Peak Power vector [W]
    p_peak_vec = [0.1, 1, 5, 10, 20, 35, 50];
    Lp = length(p_peak_vec);
    
    %% --- Fixed system parameters ---
    f = 5e9;        % Frequency [Hz]
    c = 7e8;        % Speed of light [m/s]
    v_rel = 7e3;    % Relative velocity [m/s]
    boltz_const = 1.380649e-23;  % Boltzmann constant [J/K]
    T_sys = 300;                 % System noise temp [K]
    receiver_noise_db = 2;       % [dB]
    receiver_noise    = 10^(receiver_noise_db/10);
    L_sys_db    = 5;             % System losses [dB]
    L_sys       = 10^(L_sys_db/10);
    eff_trans   = 1;             % Tx efficiency
    lambda      = c/f;           % Wavelength [m]
    quant_bits  = 2;             % bits for data quantization
    
    % Desired range resolution & pulse compression
    desired_range_resolution = 0.1;  % 10 cm
    B_eff_required = c/(2*desired_range_resolution);
    T_min = 1/B_eff_required;
    
    %% --- Create the design grid ---
    Nw   = length(antenna_width_vec); 
    Nr   = length(range_vec);        
    Nb   = length(bandwidth_vec);
    
    %% --- Initialize arrays ---
    SW_all        = zeros(Nw, Nr, Nb, Lp);
    NESZ_all      = zeros(Nw, Nr, Nb, Lp);
    data_rate_all = nan(Nw, Nr, Nb, Lp);
    mass_all      = nan(Nw, Nr, Nb, Lp);
    SNR_dB_all    = nan(Nw, Nr, Nb, Lp);
    CR_all        = nan(Nw, Nr, Nb, Lp);
    bw_used_all   = zeros(Nw, Nr, Nb, Lp);
    p_peak_store  = zeros(Nw, Nr, Nb, Lp);
    FOV_all       = nan(Nw, Nr, Nb, Lp); % store half-cone angles
    
    % Physical antenna area (depends on antenna width, but not on range/bw/power).
    A_phys = zeros(Nw, 1);
    for i = 1:Nw
        D_AT    = 2*res_along_vec; % assume along-track dimension
        A_phys(i) = antenna_width_vec(i)*D_AT;
    end
    
    %% ============================
    %  2a) Main loops over power & bandwidth
    %% ============================
    for l = 1:Lp
        p_peak_current = p_peak_vec(l);
        
        for k = 1:Nb
            bw = bandwidth_vec(k);
            t_pulse = 1/bw;
            Res_range = c*t_pulse/2;
            
            % Required pulse compression ratio
            CR_required = t_pulse / T_min;
            
            for i = 1:Nw
                for j = 1:Nr
                    R_val = range_vec(j);
                    
                    % In stripmap: swath width = lambda * R_val / antenna_width
                    SW_strip = lambda * R_val / antenna_width_vec(i);
                    
                    % In scanSAR mode: total swath = Nsub subswaths * strip width
                    if Nsub > 1
                        SW = Nsub * SW_strip;
                    else
                        SW = SW_strip;
                    end
                    SW_all(i,j,k,l) = SW;
                    
                    % Use a function for PRF bounds (some user function sarprfbounds)
                    graz_ang = lambda/2 * antenna_width_vec(i);
                    
                    if Nsub > 1
                        [PRF_min_local, PRF_max] = sarprfbounds(v_rel, res_along_vec, SW_strip, graz_ang);
                    else
                        [PRF_min_local, PRF_max] = sarprfbounds(v_rel, res_along_vec, SW, graz_ang);
                    end
                    
                    bw_used_all(i,j,k,l) = bw;
                    p_peak_store(i,j,k,l) = p_peak_current;
                    
                    if (PRF_max < PRF_min_local)
                        % Not feasible
                        NESZ_all(i,j,k,l) = NaN;
                        data_rate_all(i,j,k,l) = NaN;
                        mass_all(i,j,k,l) = NaN;
                        SNR_dB_all(i,j,k,l) = NaN;
                        CR_all(i,j,k,l) = NaN;
                        FOV_all(i,j,k,l) = NaN;
                    else
                        prf = PRF_max;
                        % In scanSAR, each subswath only sees 1 of every Nsub pulses
                        if Nsub > 1
                            duty = t_pulse * (prf / Nsub);
                        else
                            duty = t_pulse * prf;
                        end
                        P_avg = duty * p_peak_current;
                        
                        % Azimuth gain (some user function sarazgain)
                        azGain = sarazgain(R_val, lambda, v_rel, res_along_vec, prf);
                        if Nsub > 1
                            azGain = azGain / sqrt(Nsub);
                        end
                        
                        FilterGain = 1;
                        
                        % Noise Equivalent Sigma Zero
                        NESZ_all(i,j,k,l) = 10*log10( ...
                            (2*v_rel * (4*pi*R_val)^3 * boltz_const * T_sys * receiver_noise * L_sys) ...
                            / (P_avg * FilterGain * azGain * (eff_trans * 4*pi*0.7*A_phys(i)/(lambda^2))^2 * lambda^3 * Res_range) );
                        
                        % SNR
                        SNR_dB_all(i,j,k,l) = 10*log10(pi*(0.1)^2) - NESZ_all(i,j,k,l);  % for 10 cm^2 reference?
                        
                        if Nsub > 1
                            % Additional penalty for scanSAR dwell fraction
                            SNR_dB_all(i,j,k,l) = SNR_dB_all(i,j,k,l) - 10*log10(Nsub);
                        end
                        
                        % Data rate estimate
                        t_swath = 2*SW/c;  % just an example
                        data_rate_all(i,j,k,l) = 2*bw*quant_bits*t_swath*prf*duty;
                        
                        % Mass estimate
                        mass_all(i,j,k,l) = 20.5044 * antenna_width_vec(i) * (2*res_along_vec);
                        
                        % Pulse compression ratio
                        CR_all(i,j,k,l) = CR_required;
                        
                        % Half-cone angle for that swath
                        FOV_all(i,j,k,l) = asind((SW/2)/R_val);
                    end
                end
            end
        end
    end
    
    %% ============================
    %  2b) Flatten / analyze feasible designs for p_peak=10 W only
    %% ============================
    % Flatten arrays:
    SW_vec      = SW_all(:);
    FOV_vec     = FOV_all(:);
    p_peak_flat = p_peak_store(:);
    SNR_vec     = SNR_dB_all(:);
    
    feasibleIdx = ~isnan(SNR_vec) & (SNR_vec>5) & (p_peak_flat==10);
    
    % If no feasible designs, skip
    if ~any(feasibleIdx)
        fprintf('No feasible designs for Nsub=%d.\n', Nsub);
        maxFOV_vs_Nsub(idxN) = NaN;
        continue;
    end
    
    % Among feasible designs, find the maximum FOV
    FOV_feas = FOV_vec(feasibleIdx);
    maxFOV_vs_Nsub(idxN) = max(FOV_feas);
end

%% ============================
%  3) Plot the maximum feasible FOV vs the number of subswaths
% ============================
figure;
plot(Nsub_list, maxFOV_vs_Nsub, 'o-', 'LineWidth',1.5, 'MarkerSize',6);
xlabel('Number of Subswaths (N_{sub})');
ylabel('Max Feasible Half-Cone FOV (degrees)');
title('Max Feasible FOV vs. Number of Subswaths (p_{peak} = 10 W)');
grid on;
