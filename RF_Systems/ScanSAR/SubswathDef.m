%% Script: FoV and SNR vs. Number of Sub-Apertures for Different Ranges
clear; clc;

%Chosen number of subswaths = 8

%% Fixed Design Parameters
ant_width = 3;                     % [m] Fixed antenna width
N_sub_vec = 1:15;                  % Number of sub-apertures (from 1 to 15)
R_vals = [50e3, 100e3, 150e3, 200e3];  % Fixed ranges [m]

%% System Constants
f = 5e9;                    
c = 3e8;                    
lambda = c/f;               
D_AT = 0.2;                 
quantisation_bits = 10;     
bandwidth = 1e6;            
daz_req = 10e-2;            
p_peak = 25;                
T_sys = 300;                
receiver_noise_db = 2;      
receiver_noise = 10^(receiver_noise_db/10);
L_sys_db = 5;               
L_sys = 10^(L_sys_db/10);
desired_range_resolution = 0.1;         
B_eff_required = c/(2*desired_range_resolution);  
T_min = 1/B_eff_required;               

%% Preallocate arrays for Effective FoV and SNR
% Rows correspond to N_sub and columns to different ranges.
numN = length(N_sub_vec);
numR = length(R_vals);
FOV_mat = nan(numN, numR);
SNR_mat = nan(numN, numR);

%% Parameter Sweep: For each range and each sub-aperture count
for r_idx = 1:numR
    R = R_vals(r_idx);
    for i = 1:numN
        N_sub = N_sub_vec(i);
        
        % Compute physical antenna area using SARminArea
        [A_phys, ~] = SARminArea(ant_width, D_AT, R, N_sub, lambda);
        
        % Compute SNR (in dB) and swath width (SW) using SNR_model
        [SNR_dB, SW] = SNR_model(p_peak, bandwidth, R, ant_width, lambda, daz_req, ...
                                  N_sub, T_sys, receiver_noise, L_sys, T_min, A_phys, quantisation_bits);
        SNR_mat(i, r_idx) = SNR_dB;
        
        % Effective FoV computed as atand((SW/2)/R)
        FOV_mat(i, r_idx) = atand(SW/(2*R));
    end
end

%% Plot 1: Effective FoV vs. Number of Sub-Apertures
figure; hold on; grid on;
for r_idx = 1:numR
    plot(N_sub_vec, FOV_mat(:, r_idx), '-o', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Range = %.0f km', R_vals(r_idx)/1e3));
end
xlabel('Number of Sub-Apertures');
ylabel('Effective Field of View [deg]');
title('Effective FoV vs. Number of Sub-Apertures');
legend('Location','best');
hold off;

%% Plot 2: SNR vs. Number of Sub-Apertures
figure; hold on; grid on;
for r_idx = 1:numR
    plot(N_sub_vec, SNR_mat(:, r_idx), '-o', 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Range = %.0f km', R_vals(r_idx)/1e3));
end
% Add a horizontal black dotted line at SNR = 7 dB
xLimits = xlim;
plot(xLimits, [7 7], 'k--', 'LineWidth', 1.5);
text(xLimits(2)*0.95, 7+0.5, 'SNR = 7', 'HorizontalAlignment', 'right', 'Color', 'k', 'FontSize', 10);

xlabel('Number of Sub-Apertures');
ylabel('SNR [dB]');
title('SNR vs. Number of Sub-Apertures');
legend('Location','best');
hold off;

%% -------------------------
%% Support Functions (unchanged)
%% -------------------------

function [A_phys, A_min] = SARminArea(ant_width, D_AT, disp_range, N_sub, lambda)
    c = 3e8;       % Speed of light [m/s]
    v_rel = 7.5e3;
    
    A_phys = ant_width * D_AT;  % [m^2]
    if ant_width > D_AT
        max_dimension = ant_width;
    else
        max_dimension = D_AT;
    end
    beamwidth = (lambda/ max_dimension)/2;  % half beamwidth in radians
    beamwidth = rad2deg(beamwidth);
    nu_1 = beamwidth;  % Frozen Design Condition: Look angle is half beamwidth 
    incident_angle = incident_angledeg(N_sub, beamwidth, nu_1);
    
    R_val = disp_range;
    R_slant = R_val / cosd(incident_angle); % Slant Range
    A_min = (4*v_rel*lambda*R_slant)/c * tand(incident_angle); % From 'Myth of min Area paper'
    if A_phys < A_min
        A_phys = NaN;
    end
end

function [SNR_dB, SW] = SNR_model(p_peak, bandwidth, disp_range, ant_width, lambda, daz_req, Nsub, T_sys, receiver_noise, L_sys, T_min, A_phys, quantisation_bits)
    %% Constants 
    c = 3e8;  % Speed of light [m/s]
    boltz_const = 1.380649e-23;  % Boltzmann constant [J/K]
    v_rel = 7.5e3;  % Speed of radar [m/s]
    
    p_peak_current = p_peak;
    bw = bandwidth;
    t_pulse = 1/bw;  % Transmitted pulse duration [s]
    Res_range = c*t_pulse/2;  % Uncompressed range resolution [m]   
    CR_required = t_pulse / T_min;
                        
    R_val = disp_range;
    beamwidth = (lambda/ ant_width)/2;
    beamwidth = rad2deg(beamwidth);
    if beamwidth > 90
        beamwidth = NaN;
    end
    nu_1 = beamwidth;  % Pre-freeze geometry definition 
    incident_angle = incident_angledeg(Nsub, beamwidth, nu_1);
    
    R_slant = R_val/cosd(incident_angle);
    graz_ang = 90 - incident_angle;
    SW_strip = lambda * R_slant / (ant_width*sind(graz_ang)); % SMAD p.508
    if Nsub > 1
        SW = Nsub * SW_strip;
    else
        SW = SW_strip;
    end
    
    if Nsub > 1
        [PRF_min_local, PRF_max] = sarprfbounds(v_rel, daz_req, SW_strip, graz_ang);
        prf = PRF_min_local/Nsub;
    else
        [PRF_min_local, PRF_max] = sarprfbounds(v_rel, daz_req, SW, graz_ang);
        prf = PRF_min_local;
    end
    t_swath = 2*SW/c;
    
    if PRF_max < PRF_min_local
        NESZ = NaN;
        data_rate = NaN;
        mass = NaN;
        SNR_dB = NaN;
        CR = NaN;
    else
        duty = t_pulse * prf;
        P_avg = duty * p_peak_current;
        azGain = sarazgain(R_slant, lambda, v_rel, daz_req, prf);
        FilterGain = t_pulse*bw;
        NESZ = 10*log10( (2*v_rel * (4*pi*R_slant)^3 * boltz_const * T_sys * receiver_noise * L_sys) ...
                / (P_avg * FilterGain * azGain * (4*pi * 0.7 * A_phys/(lambda^2))^2 * lambda^3 * Res_range) );
        SNR_dB = 10*log10(pi*(10e-2)^2) - NESZ;
        if Nsub > 1
            SNR_dB = SNR_dB - 10*log10(Nsub);
        end
        
        data_rate = 2 * bw * quantisation_bits * t_swath * prf * duty;
        mass = 20.5044 * ant_width * (2*daz_req);
        CR = CR_required;
    end
end

function [nu] = incident_angledeg(Nsub, beamwidth, nu_1)
    nu = nu_1;
    if Nsub ~= 1
        nu = Nsub * beamwidth;
    end 
end
