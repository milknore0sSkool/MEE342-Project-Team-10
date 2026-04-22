clc
clear
close all

%% ============================================================
% UPDATED GEAR ANALYSIS - MOCK CARBURIZED MATERIAL
%% ============================================================

%% INPUT DATA

gear_names = {'Gear 1','Gear 2','Gear 3','Gear 4','Gear 5','Gear 6'};

% Teeth
N = [45, 72, 45, 72, 45, 45];

% Module [mm]
m = 3.5;

% Pressure angle
phi_deg = 20;
phi = deg2rad(phi_deg);

% Pitch diameters [mm]
d_mm = [157.5, 252, 157.5, 252, 157.5, 157.5];

% UPDATED face widths [mm]
b_mm = [55, 55, 87.5, 87.5, 87.5, 87.5];

% Torque distribution [N*m]
T_Nm = [300, 480, 480, 768, 768, 768];

% Gear meshes
pairs = [1 2;
         3 4;
         5 6];

pair_names = {'Mesh 1-2','Mesh 3-4','Mesh 5-6'};

%% ============================================================
% MATERIAL (UPDATED)
%% ============================================================

Sut = 1234;   % MPa
Sy  = 986;    % MPa

% Allowables (same model you were using)
sigma_b_allow = 0.30 * Sut;   % bending allowable
sigma_c_allow = 1100;         % contact allowable

%% ============================================================
% FACTORS (kept simple)
%% ============================================================

K_total = 1.0;

%% ============================================================
% LEWIS FORM FACTOR
%% ============================================================

Y = 0.154 - 0.912 ./ N;

%% ============================================================
% TANGENTIAL LOAD
%% ============================================================

d_m = d_mm / 1000;
Wt = 2 .* T_Nm ./ d_m;

%% ============================================================
% BENDING STRESS
%% ============================================================

sigma_b = (Wt .* K_total) ./ (b_mm .* m .* Y);
FoS_bending = sigma_b_allow ./ sigma_b;

%% ============================================================
% CONTACT STRESS
%% ============================================================

Ze = 191;   % steel-steel

sigma_c = zeros(size(pairs,1),1);
FoS_contact = zeros(size(pairs,1),1);
I_geom = zeros(size(pairs,1),1);

for i = 1:size(pairs,1)

    p = pairs(i,1);
    g = pairs(i,2);

    Q = N(g) / N(p);
    I_geom(i) = (cos(phi) * sin(phi) / 2) * (Q / (Q + 1));

    Wt_pair = Wt(p);

    b_pair = b_mm(p);
    d_pair = d_mm(p);

    sigma_c(i) = Ze * sqrt( (Wt_pair * K_total) / (b_pair * d_pair * I_geom(i)) );
    FoS_contact(i) = sigma_c_allow / sigma_c(i);

end

%% ============================================================
% PRINT RESULTS
%% ============================================================

fprintf('\n====================================================\n');
fprintf('UPDATED GEAR ANALYSIS (MOCK CARBURIZED)\n');
fprintf('====================================================\n');

fprintf('Allowable bending stress = %.2f MPa\n\n', sigma_b_allow);

for i = 1:length(N)

    fprintf('%s:\n', gear_names{i});
    fprintf('  Bending stress = %.2f MPa\n', sigma_b(i));
    fprintf('  Bending FoS    = %.2f\n\n', FoS_bending(i));

end

fprintf('====================================================\n');
fprintf('CONTACT RESULTS\n');
fprintf('====================================================\n');

for i = 1:size(pairs,1)

    fprintf('%s:\n', pair_names{i});
    fprintf('  Contact stress = %.2f MPa\n', sigma_c(i));
    fprintf('  Contact FoS    = %.2f\n\n', FoS_contact(i));

end

%% ============================================================
% CRITICAL RESULTS
%% ============================================================

[min_bend_FoS, idx_bend] = min(FoS_bending);
[min_cont_FoS, idx_cont] = min(FoS_contact);

fprintf('====================================================\n');
fprintf('CRITICAL RESULTS\n');
fprintf('====================================================\n');

fprintf('Critical gear (bending): %s\n', gear_names{idx_bend});
fprintf('Min bending FoS        : %.2f\n\n', min_bend_FoS);

fprintf('Critical mesh (contact): %s\n', pair_names{idx_cont});
fprintf('Min contact FoS        : %.2f\n', min_cont_FoS);

%% ============================================================
% PLOTS
%% ============================================================

figure
bar(sigma_b)
title('Gear Bending Stress')
ylabel('MPa')
set(gca,'XTickLabel',gear_names)
grid on

figure
bar(FoS_bending)
title('Bending Factor of Safety')
ylabel('FoS')
set(gca,'XTickLabel',gear_names)
grid on

figure
bar(sigma_c)
title('Contact Stress')
ylabel('MPa')
set(gca,'XTickLabel',pair_names)
grid on

figure
bar(FoS_contact)
title('Contact Factor of Safety')
ylabel('FoS')
set(gca,'XTickLabel',pair_names)
grid on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% S-N Diagrams %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============================================================
% S-N DIAGRAM FOR GEARS (FATIGUE)
%% ============================================================

figure
hold on
grid on
set(gca,'XScale','log','YScale','log')

title('S-N Diagram for Gear Material (Mock Carburized)')
xlabel('Number of Cycles (N)')
ylabel('Stress Amplitude [MPa]')

% Material properties
Sut = 1234;   % MPa

% Define key fatigue points
N1 = 1e3;
S1 = 0.9 * Sut;

N2 = 1e6;
Se = 0.5 * Sut;

% Create S-N curve
N_curve = logspace(3,7,200);

% Basquin equation fit
b = (log10(Se) - log10(S1)) / (log10(N2) - log10(N1));
a = S1 / (N1^b);

S_curve = a * (N_curve.^b);

% Plot S-N curve
plot(N_curve, S_curve, 'LineWidth', 2)

% Plot endurance limit line
yline(Se, '--', 'Endurance Limit')

%% ============================================================
% PLOT EACH GEAR OPERATING POINT
%% ============================================================

% Assume design life (you can change this)
N_operating = 1e6;

for i = 1:length(sigma_b)
    
    scatter(N_operating, sigma_b(i), 80, 'filled')
    text(N_operating*1.1, sigma_b(i), gear_names{i})
    
end

legend('S-N Curve','Endurance Limit','Gear Operating Points')
