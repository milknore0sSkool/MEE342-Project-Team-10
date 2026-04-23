clc
clear
close all

%% ============================================================
% SHIGLEY / AGMA-STYLE SPUR GEAR ANALYSIS
% Units:
%   Length  -> mm
%   Force   -> N
%   Torque  -> N*m
%   Stress  -> MPa
%   Speed   -> rpm
%
% NOTE:
% 1) This script is structured around Shigley/AGMA stress equations.
% 2) You MUST supply J, St, Sc, YN, ZN from your Shigley charts/tables.
% 3) If you know your gear quality Qv and speeds, the helper at the end
%    can estimate Kv automatically. Otherwise enter Kv manually.
%% ============================================================

%% ============================================================
% INPUT DATA
%% ============================================================

gear_names = {'Gear 1','Gear 2','Gear 3','Gear 4','Gear 5','Gear 6'};

% Teeth
N = [45, 72, 45, 72, 45, 45];

% Module [mm]
m = 3.5;

% Pressure angle [deg]
phi_deg = 20;
phi = deg2rad(phi_deg);

% Pitch diameters [mm]
d_mm = [157.5, 252, 157.5, 252, 157.5, 157.5];

% Face widths [mm]
F_mm = [55, 55, 87.5, 87.5, 87.5, 87.5];

% Torque on each gear [N*m]
T_Nm = [300, 480, 480, 768, 768, 768];

% Optional speeds for each gear [rpm]
% Put your actual speeds here if you want automatic Kv from Qv
n_rpm = [0, 0, 0, 0, 0, 0];

% Gear meshes [pinion gear]
pairs = [1 2;
         3 4;
         5 6];

pair_names = {'Mesh 1-2','Mesh 3-4','Mesh 5-6'};

%% ============================================================
% MATERIAL / STRENGTH INPUTS FROM SHIGLEY TABLES / CHARTS
%% ============================================================
% These are NOT Sut or Sy.
% These are AGMA gear strength numbers from Shigley tables/figures.

% AGMA bending strength number St for each gear [MPa]
% Replace with your actual values from Shigley based on material / heat treat
St = [517, 517, 517, 517, 517, 517];   % MPa, carburized and hardened steel



% AGMA contact strength number Sc for each gear [MPa]
% Replace with your actual values from Shigley based on material / heat treat
Sc = [1170, 1170, 1170, 1170, 1170, 1170];
% Elastic coefficient Ze [sqrt(MPa)]
% steel-steel is commonly taken around this value
Ze = 191;

%% ============================================================
% GEOMETRY FACTORS FROM SHIGLEY / AGMA
%% ============================================================
% J must come from the spur-gear J chart (Shigley Fig. 14-6 / AGMA data).
% Put one J value per gear.
J = [0.33, 0.38, 0.33, 0.38, 0.33, 0.33];

% I can be computed for standard external spur gears:
% I = (cos(phi)*sin(phi)/2) * (Q/(Q+1))
% where Q = NG / NP for each mesh

%% ============================================================
% LIFE / RELIABILITY / TEMPERATURE / RIM FACTORS
%% ============================================================
% YN and ZN come from Shigley life-factor charts.
% Replace with chart values for your required life.
YN = [1, 1, 1, 1, 1, 1];   % bending stress-cycle factor
ZN = [1, 1, 1, 1, 1, 1];   % contact stress-cycle factor

% Reliability factor
KR = 1.00;   % 0.99 reliability in Shigley/AGMA baseline

% Temperature factor
KT = 1.00;   % use 1.0 if temperature is within normal range

% Rim-thickness factor for bending
KB = [1, 1, 1, 1, 1, 1];   % use >1 only for thin-rim gears

% Hardness-ratio factor for pitting
CH = [1, 1, 1];            % one per mesh

%% ============================================================
% LOAD MODIFICATION FACTORS
%% ============================================================
% Enter these based on your application and mounting.
% One value per mesh unless you want gear-by-gear size factors.

% Overload factor Ko
Ko = [1.25, 1.25, 1.25];

% Dynamic factor Kv
% Option A: manually enter
Kv = [1.10, 1.10, 1.10];

% Option B: if you know Qv and pinion speed, uncomment below
% Qv = 8;
% for k = 1:size(pairs,1)
%     p = pairs(k,1);
%     Kv(k) = agma_Kv_from_Qv(d_mm(p), n_rpm(p), Qv);
% end

% Size factor Ks
% Can be taken as 1.0 if no detrimental size effect is assumed
Ks = [1, 1, 1, 1, 1, 1];

% Load-distribution factor Km
Km = [1.30, 1.30, 1.30];

% Surface condition factor for contact
Cf = [1.00, 1.00, 1.00];

%% ============================================================
% DERIVED QUANTITIES
%% ============================================================

% Tangential transmitted load on each gear [N]
d_m = d_mm / 1000;                 % convert to m
Wt = 2 .* T_Nm ./ d_m;             % N

% Geometry factor I for each mesh
I_geom = zeros(size(pairs,1),1);

for k = 1:size(pairs,1)
    p = pairs(k,1);
    g = pairs(k,2);

    Q = N(g) / N(p);
    I_geom(k) = (cos(phi) * sin(phi) / 2) * (Q / (Q + 1));
end

%% ============================================================
% AGMA / SHIGLEY BENDING STRESS
% sigma_b = Wt*Ko*Kv*Ks*Km*KB / (F*m*J)
% Since N/mm^2 = MPa, mm-units are consistent here.
%% ============================================================

sigma_b = zeros(length(N),1);
sigma_b_all = zeros(length(N),1);
SF_bending = zeros(length(N),1);

for i = 1:length(N)

    % Find which mesh this gear belongs to
    mesh_idx = find(any(pairs == i, 2), 1, 'first');

    sigma_b(i) = (Wt(i) * Ko(mesh_idx) * Kv(mesh_idx) * Ks(i) * ...
                  Km(mesh_idx) * KB(i)) / (F_mm(i) * m * J(i));

    % AGMA allowable bending stress
    sigma_b_all(i) = (St(i) * YN(i)) / (KT * KR);

    % AGMA bending safety factor
    SF_bending(i) = sigma_b_all(i) / sigma_b(i);
end

%% ============================================================
% AGMA / SHIGLEY CONTACT STRESS
% sigma_c = Ze * sqrt( Wt*Ko*Kv*Ks*Km*Cf / (dp*F*I) )
%
% Use pinion diameter and pinion transmitted load for each mesh.
%% ============================================================

sigma_c = zeros(size(pairs,1),1);

% Two ways to report safety for pitting:
% 1) AGMA-style pitting safety factor (squared ratio)
% 2) Plain stress ratio
SH_contact_pinion = zeros(size(pairs,1),1);
SH_contact_gear   = zeros(size(pairs,1),1);
SR_contact_pinion = zeros(size(pairs,1),1);
SR_contact_gear   = zeros(size(pairs,1),1);

sigma_c_all_pinion = zeros(size(pairs,1),1);
sigma_c_all_gear   = zeros(size(pairs,1),1);

for k = 1:size(pairs,1)

    p = pairs(k,1);
    g = pairs(k,2);

    Wt_pair = Wt(p);
    dp = d_mm(p);
    Fp = F_mm(p);

    sigma_c(k) = Ze * sqrt( (Wt_pair * Ko(k) * Kv(k) * Ks(p) * Km(k) * Cf(k)) / ...
                            (dp * Fp * I_geom(k)) );

    % AGMA allowable contact stress numbers
    sigma_c_all_pinion(k) = (Sc(p) * ZN(p) * CH(k)) / (KT * KR);
    sigma_c_all_gear(k)   = (Sc(g) * ZN(g) * CH(k)) / (KT * KR);

    % Plain stress ratio
    SR_contact_pinion(k) = sigma_c_all_pinion(k) / sigma_c(k);
    SR_contact_gear(k)   = sigma_c_all_gear(k)   / sigma_c(k);

    % AGMA-style pitting safety factor
    SH_contact_pinion(k) = (sigma_c_all_pinion(k) / sigma_c(k))^2;
    SH_contact_gear(k)   = (sigma_c_all_gear(k)   / sigma_c(k))^2;
end

%% ============================================================
% PRINT RESULTS
%% ============================================================

fprintf('\n====================================================\n');
fprintf('SHIGLEY / AGMA SPUR GEAR ANALYSIS\n');
fprintf('====================================================\n\n');

fprintf('---------------- BENDING RESULTS -------------------\n');
for i = 1:length(N)
    fprintf('%s:\n', gear_names{i});
    fprintf('  AGMA bending stress, sigma_b      = %.2f MPa\n', sigma_b(i));
    fprintf('  Allowable bending stress          = %.2f MPa\n', sigma_b_all(i));
    fprintf('  AGMA bending safety factor, SF    = %.2f\n\n', SF_bending(i));
end

fprintf('---------------- CONTACT RESULTS -------------------\n');
for k = 1:size(pairs,1)
    p = pairs(k,1);
    g = pairs(k,2);

    fprintf('%s:\n', pair_names{k});
    fprintf('  Contact stress, sigma_c                 = %.2f MPa\n', sigma_c(k));
    fprintf('  Pinion allowable contact stress         = %.2f MPa\n', sigma_c_all_pinion(k));
    fprintf('  Gear allowable contact stress           = %.2f MPa\n', sigma_c_all_gear(k));
    fprintf('  Pinion stress ratio (allow/stress)      = %.2f\n', SR_contact_pinion(k));
    fprintf('  Gear stress ratio   (allow/stress)      = %.2f\n', SR_contact_gear(k));
    fprintf('  Pinion AGMA pitting safety factor, SH   = %.2f\n', SH_contact_pinion(k));
    fprintf('  Gear AGMA pitting safety factor, SH     = %.2f\n\n', SH_contact_gear(k));
end

[min_SF_bend, idx_bend] = min(SF_bending);
[min_SH_pin, idx_SH_pin] = min(SH_contact_pinion);
[min_SH_gear, idx_SH_gear] = min(SH_contact_gear);

fprintf('---------------- CRITICAL RESULTS ------------------\n');
fprintf('Minimum bending safety factor: %.2f (%s)\n', min_SF_bend, gear_names{idx_bend});
fprintf('Minimum pinion pitting safety factor: %.2f (%s)\n', min_SH_pin, pair_names{idx_SH_pin});
fprintf('Minimum gear pitting safety factor: %.2f (%s)\n', min_SH_gear, pair_names{idx_SH_gear});

%% ============================================================
% PLOTS
%% ============================================================

figure
bar(sigma_b)
title('AGMA Bending Stress')
ylabel('Stress [MPa]')
set(gca,'XTickLabel',gear_names)
grid on

figure
bar(SF_bending)
title('AGMA Bending Safety Factor')
ylabel('SF')
set(gca,'XTickLabel',gear_names)
grid on

figure
bar(sigma_c)
title('AGMA Contact Stress')
ylabel('Stress [MPa]')
set(gca,'XTickLabel',pair_names)
grid on

figure
bar([SH_contact_pinion, SH_contact_gear])
title('AGMA Pitting Safety Factor')
ylabel('SH')
legend('Pinion','Gear','Location','best')
set(gca,'XTickLabel',pair_names)
grid on

%% ============================================================
% OPTIONAL HELPER: AGMA DYNAMIC FACTOR FROM Qv
% V is converted to ft/min because the standard equation is commonly
% written that way.
%% ============================================================
function Kv = agma_Kv_from_Qv(d_mm, n_rpm, Qv)

    % pitch-line velocity in m/s
    V_ms = pi * (d_mm/1000) * n_rpm / 60;

    % convert to ft/min
    V_ftmin = V_ms * 196.8504;

    B = 0.25 * (12 - Qv)^(2/3);
    A = 50 + 56 * (1 - B);

    Kv = ((A + sqrt(V_ftmin)) / A)^B;
end