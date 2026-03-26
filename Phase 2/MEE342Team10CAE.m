clc; clear; close all;

%% ================= INPUTS =================

% ---- Material (AISI 4140 normalized) ----
Sut = 655;                  % MPa (MatWeb)
HB  = Sut / 3.45;           % ~190 HB approximation

HBP = HB;                   % pinion hardness
HBG = HB;                   % gear hardness

% ---- Gear Geometry ----
m   = 1.5;                  % module (mm)
Np  = 60;                   % pinion teeth
Ng  = 30;                   % gear teeth
phi = 20;                   % pressure angle (deg)
b   = 20;                   % face width (mm)

% ---- Derived Geometry ----
dp = m * Np;                % pinion pitch diameter (mm)
dg = m * Ng;                % gear pitch diameter (mm)
mG = Ng / Np;               % gear ratio (NOW 0.5)

% ---- Load ----
M = 105000;              % torque (N*mm)

% ---- Operating Conditions ----
rpm   = 1200;
hours = 1000;

% ---- AGMA Factors ----
SF = 2.0;                   % bending safety factor
SH = 2.0;                   % contact safety factor
KT = 1.0;                   % temperature factor
KR = 0.9;                   % reliability factor

% ---- Assumed AGMA Modifiers ----
Ko = 1.25;     % overload factor
Kv = 1.1;      % dynamic factor
Ks = 1.0;      % size factor
Km = 1.3;      % load distribution
Kb = 1.0;      % rim thickness

Cp = 191;      % elastic coefficient (MPa^0.5 for steel-steel)

%% ================= CYCLES =================

N = rpm * 60 * hours;

%% ================= ALLOWABLE STRESSES =================

% ---- BENDING ----
St = 0.533*HB + 88.3;

if N <= 1e7
    YN = 1.3558 * N^(-0.0178);
else
    YN = 1.0;
end

sigma_allow_b = (St * YN) / (KT * KR * SF);

% ---- CONTACT ----
Sc = 2.22*HB + 200;

if N <= 1e7
    ZN = 1.4488 * N^(-0.023);
else
    ZN = 1.0;
end

% Hardness ratio factor
ratio = HBP / HBG;
Aprime = 8.98e-3 * ratio - 8.29e-3;
CH = 1.0 + Aprime * (mG - 1.0);   % NOTE: <1 now

sigma_allow_c = (Sc * ZN * CH) / (KT * KR * SH);

%% ================= ACTUAL STRESSES =================

% Tangential force at pitch circle
Wt = 2*M / dp;     % N

% ---- Geometry Factors ----
J = 0.38;   % improved bending factor for ~60 teeth

I = (cosd(phi)*sind(phi)/2) * (mG/(mG+1));

% ---- BENDING STRESS ----
sigma_b = (Wt * Ko * Kv * Ks * Km * Kb) / (b * m * J);

% ---- CONTACT STRESS ----
sigma_c = Cp * sqrt( (Wt * Ko * Kv * Ks * Km) / (b * dp * I) );

%% ================= RESULTS =================

fprintf('===== AGMA GEAR ANALYSIS =====\n\n');

fprintf('--- BASIC INFO ---\n');
fprintf('HB: %.2f\n', HB);
fprintf('Cycles (N): %.3e\n', N);
fprintf('Gear Ratio (mG): %.2f\n', mG);
fprintf('Tangential Force Wt: %.2f N\n\n', Wt);

fprintf('--- BENDING ---\n');
fprintf('Actual Stress: %.2f MPa\n', sigma_b);
fprintf('Allowable Stress: %.2f MPa\n', sigma_allow_b);

if sigma_b < sigma_allow_b
    fprintf('STATUS: SAFE (BENDING)\n\n');
else
    fprintf('STATUS: FAIL (BENDING)\n\n');
end

fprintf('--- CONTACT (PITTING) ---\n');
fprintf('Actual Stress: %.2f MPa\n', sigma_c);
fprintf('Allowable Stress: %.2f MPa\n', sigma_allow_c);
fprintf('CH Factor: %.4f\n', CH);

if sigma_c < sigma_allow_c
    fprintf('STATUS: SAFE (CONTACT)\n\n');
else
    fprintf('STATUS: FAIL (CONTACT)\n\n');
end

fprintf('==============================\n');



Tmax_b = (sigma_allow_b * b * m * J * dp) / (2 * Ko * Kv * Ks * Km * Kb);   % N*mm
Tmax_c = ((sigma_allow_c / Cp)^2 * b * dp * I) / (Ko * Kv * Ks * Km);
