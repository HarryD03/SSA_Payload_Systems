%% Script: Effect of Input Power on Measurement Uncertainty for Different Frequencies 
%% with PRF Computed via sarprfbounds, Including Pulse Compression, and FOV vs. Power Plot
clear; clc; close all;

%% Operational Mode Settings (scanSAR)
operational_mode = 'scanSAR';
if strcmpi(operational_mode, 'scanSAR')
    Nsub = 8;  % Number of subswaths in scanSAR mode
else
    Nsub = 1;
end

%% Fixed Design Parameters (common for all frequency cases)
antenna_width = 1.67;       % [m] Fixed antenna width
D_AT = 0.2;                 % [m] Along-track aperture (sets antenna length)
res_along = D_AT / 2;       % [m] Along-track resolution (here 0.1 m)

bw = 0.5e6;                 % [Hz] Fixed bandwidth
R_val = 50e3;               % [m] Fixed slant range (50 km)

%% Fixed System Parameters
v_rel       = 7.5e3;          % [m/s] Relative velocity
boltz_const = 1.380649e-23;   % [J/K] Boltzmann constant
T_sys       = 300;            % [K] System noise temperature - Standard
receiver_noise_db = 2;        % [dB]
receiver_noise = 10^(receiver_noise_db/10);
L_sys_db    = 5;              % [dB] System Noise - Standard 
L_sys       = 10^(L_sys_db/10);
eff_trans   = 1;              % Transmitter efficiency
A_phys      = antenna_width * D_AT;  % Physical antenna area [m^2]

%% Derived Common Parameter and Pulse Compression Calculation
t_pulse = 1 / bw;              % Transmitted pulse duration [s]

% Define desired range resolution and compute the minimum pulse width required:
desired_range_resolution = 0.1;          % [m] (10 cm resolution)
B_eff_required = 3e8 / (2 * desired_range_resolution);  % Effective bandwidth [Hz]
T_min = 1 / B_eff_required;              % Minimum pulse duration [s]
CR_required = t_pulse / T_min;           % Pulse Compression Ratio
% Use a filter gain of 1 here (or alternatively CR_required) if desired:
FilterGain = 1;

%% Define the frequencies to compare
freqVec = [5e9, 40e9, 70e9];  % Frequencies in Hz
nFreq = length(freqVec);

%% Define a vector of input peak power values (in Watts)
p_peak_vec = linspace(5, 50, 10);  % from 5W to 50W
numP = length(p_peak_vec);

%% Preallocate storage for uncertainty and FOV for each frequency scenario
uncertainty_all = zeros(nFreq, numP);
SNR_dB_all = zeros(nFreq, numP);
FOV_all = zeros(nFreq, numP);

%% Loop over the different frequencies
for f_idx = 1:nFreq
    f = freqVec(f_idx);
    lambda = 3e8 / f;  % Wavelength [m]
    
    % Compute geometry-dependent parameters
    graz_ang = (lambda/2) * antenna_width;  % approximate grazing angle
    graz_ang = 10;
    % Compute strip swath width:
    R_slant = R_val*cosd(graz_ang);
    SW_strip = lambda * R_slant / antenna_width;
    
    % Compute PRF bounds via the sarprfbounds function:
    [PRF_min_local, PRF_max] = sarprfbounds(v_rel, res_along, SW_strip, graz_ang);
    % For scanSAR mode, each subswath is illuminated less frequently:
    if Nsub > 1
        prf = PRF_min_local / Nsub;
    else
        prf = PRF_min_local;
    end
    
    % Compute swath width (for reference) and dwell time:
    if Nsub > 1
        SW = Nsub * SW_strip;
    else
        SW = SW_strip;
    end
    t_swath = 2 * SW / 3e8;  % Swath dwell time
    
    % Compute the field-of-view (FOV) [in degrees]
    % FOV = asind((SW/2) / R_val)
    R_slant = R_val*cosd(graz_ang);
    FOV = asind((SW/2) / R_slant);
    % Create a vector for FOV (constant across power values)
    FOV_all(f_idx, :) = repmat(FOV, 1, numP);
    
    % Compute the azimuth gain using a simplified function
    azGain = sarazgain(R_slant, lambda, v_rel, res_along, prf);
    if Nsub > 1
        azGain = azGain / sqrt(Nsub);
    end
    % Duty cycle based on computed PRF:
    duty = t_pulse * prf;
    
    % Loop over input peak power values for the current frequency
    for idx = 1:numP
        % Here, note that the peak power is scaled by 0.5 as in your code.
        p_peak = p_peak_vec(idx) * 0.5;
        
        % Compute average transmitted power (accounting for duty cycle)
        P_avg = duty * p_peak;
        R_slant = R_val * cosd(graz_ang);
        % Compute NESZ (Noise Equivalent Sigma Zero)
        numerator = 2 * v_rel * (4*pi*R_slant)^3 * boltz_const * T_sys * receiver_noise * L_sys;
        denominator = P_avg * FilterGain * azGain * (eff_trans * 4*pi * 0.7 * A_phys/(lambda^2))^2 * lambda^3 * (3e8 * t_pulse / 2);
        NESZ = numerator / denominator;
        NESZ_dB = 10*log10(NESZ);
        
        % Compute SNR in dB (using a reference area: pi*(0.1)^2)
        SNR_dB = 10*log10(pi * (0.1)^2) - NESZ_dB;
        SNR_dB_all(f_idx, idx) = SNR_dB;
        SNR_linear = 10^(SNR_dB/10);
        
        % Compute range uncertainty (σᵣ) and angular uncertainty (σₐₙg)
        sigma_r = 3e8 / (2 * bw * sqrt(2 * SNR_linear));
        beamwidth = lambda / D_AT;  % Synthetic aperture beamwidth
        sigma_ang = beamwidth / sqrt(2 * SNR_linear);
        
        % Average uncertainty (simple mean of range and angular uncertainties)
        uncertainty_all(f_idx, idx) = (sigma_r + sigma_ang) / 2;
    end
end

%% Plot: Measurement Uncertainty vs. Input Peak Power for Each Frequency
figure;
colors = lines(nFreq);
for f_idx = 1:nFreq
    plot(p_peak_vec, uncertainty_all(f_idx, :), 'o-', 'LineWidth', 1.5, 'Color', colors(f_idx,:));
    hold on;
end
xlabel('Input Peak Power (W)');
ylabel('Measurement Uncertainty (\sigma) [m]');
title('Measurement Uncertainty vs. Input Peak Power for Different Frequencies');
legend('5 GHz', '40 GHz', '70 GHz', 'Location', 'best');
grid on;
hold off;

%% Plot: SNR vs. Input Peak Power for Each Frequency (Optional)
figure;
for f_idx = 1:nFreq
    plot(p_peak_vec, SNR_dB_all(f_idx, :), 's-', 'LineWidth', 1.5, 'Color', colors(f_idx,:));
    hold on;
end
xlabel('Input DC Power (W)');
ylabel('SNR (dB)');
title('SNR vs. Input Peak Power for Different Frequencies');
legend('5 GHz', '40 GHz', '70 GHz', 'Location', 'best');
grid on;
hold off;

%% Plot: FOV vs. Input Peak Power for Each Frequency
% In this model, FOV is determined solely by the geometry and frequency.
% Therefore, for each frequency, FOV remains constant as a function of power.
figure;
for f_idx = 1:nFreq
    plot(p_peak_vec, FOV_all(f_idx, :), 'd-', 'LineWidth', 1.5, 'Color', colors(f_idx,:));
    hold on;
end
xlabel('Input DC Power (W)');
ylabel('Field-of-View (deg)');
title('FOV vs. Input DC Power at Allocated Frequencies');
legend('5 GHz', '40 GHz', '70 GHz', 'Location', 'best');
grid on;
hold off;

 