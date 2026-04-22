clc
clear
close all

%% ============================================================
% MATERIAL + SECTION PROPERTIES
% =============================================================
Do = 30e-3;                 % outer diameter [m]
Di = 12e-3;                 % inner diameter [m]
ro = Do/2;

I = (pi/64)*(Do^4 - Di^4);  % area moment of inertia [m^4]
J = (pi/32)*(Do^4 - Di^4);  % polar moment of inertia [m^4]

E  = 200e9;                 % elastic modulus [Pa]
Sy = 655e6;                 % yield strength [Pa] for 4140 steel
Sut = 900e6;                % approximate ultimate strength [Pa]

%% ============================================================
% GEAR TRAIN DATA
% =============================================================
% Pitch diameters [m]
d1 = 157.5e-3;
d2 = 252e-3;
d3 = 157.5e-3;
d4 = 252e-3;
d5 = 157.5e-3;
d6 = 157.5e-3;

% Torque flow [N*m]
T1 = 300;
T2 = T1*(d2/d1);    % 480
T3 = T2;            % 480
T4 = T3*(d4/d3);    % 768
T5 = T4;            % 768
T6 = T5*(d6/d5);    % 768

% Tangential gear force magnitudes [N]
F1 = 2*T1/d1;
F2 = 2*T1/d1;
F3 = 2*T3/d3;
F4 = 2*T3/d3;
F5 = 2*T5/d5;
F6 = 2*T6/d6;

%% ============================================================
% SHAFT DEFINITIONS
% =============================================================

% -------------------------
% MIDDLE SHAFT
% -------------------------
mid.name = 'Middle Shaft';
mid.L = 0.550;
mid.supports = [0.000, 0.365, 0.550];
mid.gear_positions = [0.160, 0.310, 0.420];

% Force sign convention for top-view bending
% Change these if your physical mesh directions are opposite
mid.gear_forces = [-F1, +F4, -F5];

% Torque events: [position, value]
% + torque enters shaft, - torque leaves shaft
mid.torque_events = [
    0.000, +T1;
    0.160, -T1;
    0.310, +T4;
    0.420, -T5
];

% -------------------------
% BOTTOM SHAFT
% -------------------------
bot.name = 'Bottom Shaft';
bot.L = 0.2825;
bot.supports = [0.000, 0.2825];
bot.gear_positions = [0.030, 0.165];
bot.gear_forces = [-F2, +F3];
bot.torque_events = [
    0.030, +T2;
    0.165, -T3
];

% -------------------------
% TOP SHAFT
% -------------------------
top.name = 'Top Shaft';
top.L = 0.4575;
top.supports = [0.010, 0.4525];
top.gear_positions = [0.4575];
top.gear_forces = [-F6];
top.torque_events = [
    0.000, +T6;
    0.4575, -T6
];

%% ============================================================
% ANALYZE ALL SHAFTS
% =============================================================
mid = analyze_shaft(mid, E, I, J, ro, Sy, Sut, Do);
bot = analyze_shaft(bot, E, I, J, ro, Sy, Sut, Do);
top = analyze_shaft(top, E, I, J, ro, Sy, Sut, Do);

%% ============================================================
% PRINT RESULTS
% =============================================================
print_results(mid);
print_results(bot);
print_results(top);

% Critical shaft summary
all_yield = [mid.FoS_yield_min, bot.FoS_yield_min, top.FoS_yield_min];
all_fat   = [mid.FoS_fatigue_min, bot.FoS_fatigue_min, top.FoS_fatigue_min];
names = {mid.name, bot.name, top.name};

[global_yield_min, idx_y] = min(all_yield);
[global_fat_min, idx_f]   = min(all_fat);

fprintf('\n====================================================\n');
fprintf('SYSTEM SUMMARY\n');
fprintf('====================================================\n');
fprintf('Lowest yielding FoS (von Mises): %.3f  --> %s\n', global_yield_min, names{idx_y});
fprintf('Lowest fatigue FoS (Modified Goodman): %.3f  --> %s\n', global_fat_min, names{idx_f});

%% ============================================================
% PLOTS FOR EACH SHAFT
% =============================================================
plot_shaft_results(mid);
plot_layout(mid);

plot_shaft_results(bot);
plot_layout(bot);

plot_shaft_results(top);
plot_layout(top);

%% ============================================================
% INDIVIDUAL S-N DIAGRAMS
% =============================================================
plot_SN_single(mid, Sut, Do);
plot_SN_single(bot, Sut, Do);
plot_SN_single(top, Sut, Do);

%% ============================================================
% COMBINED S-N DIAGRAM
% =============================================================
plot_SN_combined(mid, bot, top, Sut, Do);

%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function shaft = analyze_shaft(shaft, E, I, J, ro, Sy, Sut, Do)

    % -------------------------
    % Build beam FEM for vertical reactions
    % Bearings = vertical displacement fixed only, rotation free
    % -------------------------
    xn = unique(sort([shaft.supports(:); shaft.gear_positions(:)]))';
    nn = length(xn);
    ndof = 2*nn;

    K = zeros(ndof);
    F = zeros(ndof,1);

    % Assemble beam stiffness
    for e = 1:(nn-1)
        Le = xn(e+1) - xn(e);

        ke = (E*I/Le^3) * ...
            [12      6*Le    -12      6*Le;
             6*Le    4*Le^2  -6*Le    2*Le^2;
            -12     -6*Le     12     -6*Le;
             6*Le    2*Le^2  -6*Le    4*Le^2];

        dofs = [2*e-1, 2*e, 2*(e+1)-1, 2*(e+1)];
        K(dofs,dofs) = K(dofs,dofs) + ke;
    end

    % Apply point loads at gear nodes
    for i = 1:length(shaft.gear_positions)
        pos = shaft.gear_positions(i);
        idx = find(abs(xn - pos) < 1e-12, 1);
        F(2*idx - 1) = F(2*idx - 1) + shaft.gear_forces(i);
    end

    % Bearing constraints: vertical displacement fixed
    support_node_idx = zeros(size(shaft.supports));
    for i = 1:length(shaft.supports)
        support_node_idx(i) = find(abs(xn - shaft.supports(i)) < 1e-12, 1);
    end
    fixed_dofs = 2*support_node_idx - 1;
    free_dofs = setdiff(1:ndof, fixed_dofs);

    U = zeros(ndof,1);
    U(free_dofs) = K(free_dofs, free_dofs)\F(free_dofs);

    R = K*U - F;

    shaft.reaction_positions = shaft.supports;
    shaft.reactions = zeros(size(shaft.supports));
    for i = 1:length(shaft.supports)
        shaft.reactions(i) = R(2*support_node_idx(i) - 1);
    end

    % -------------------------
    % Continuous diagrams from reactions + gear forces + torque events
    % -------------------------
    shaft.x = linspace(0, shaft.L, 5000);
    x = shaft.x;

    shaft.V = zeros(size(x));
    shaft.M = zeros(size(x));
    shaft.T = zeros(size(x));
    shaft.sigma_b = zeros(size(x));
    shaft.tau_t = zeros(size(x));
    shaft.sigma_vm = zeros(size(x));
    shaft.FoS_yield = nan(size(x));
    shaft.FoS_fatigue = nan(size(x));

    for k = 1:length(x)
        xk = x(k);

        % Shear from reactions + gear forces
        Vk = 0;
        for i = 1:length(shaft.supports)
            if xk >= shaft.supports(i)
                Vk = Vk + shaft.reactions(i);
            end
        end
        for i = 1:length(shaft.gear_positions)
            if xk >= shaft.gear_positions(i)
                Vk = Vk + shaft.gear_forces(i);
            end
        end
        shaft.V(k) = Vk;

        % Moment from reactions + gear forces
        Mk = 0;
        for i = 1:length(shaft.supports)
            if xk >= shaft.supports(i)
                Mk = Mk + shaft.reactions(i)*(xk - shaft.supports(i));
            end
        end
        for i = 1:length(shaft.gear_positions)
            if xk >= shaft.gear_positions(i)
                Mk = Mk + shaft.gear_forces(i)*(xk - shaft.gear_positions(i));
            end
        end
        shaft.M(k) = Mk;

        % Torque from torque events
        Tk = 0;
        for i = 1:size(shaft.torque_events,1)
            if xk >= shaft.torque_events(i,1)
                Tk = Tk + shaft.torque_events(i,2);
            end
        end
        shaft.T(k) = Tk;

        % Stresses
        sigma_b = shaft.M(k)*ro/I;
        tau_t   = shaft.T(k)*ro/J;

        shaft.sigma_b(k) = sigma_b;
        shaft.tau_t(k)   = tau_t;
        shaft.sigma_vm(k)= sqrt(sigma_b^2 + 3*tau_t^2);

        % Yield FoS using von Mises
        if shaft.sigma_vm(k) > 1e-6
            shaft.FoS_yield(k) = Sy/shaft.sigma_vm(k);
        end
    end

    % -------------------------
    % Fatigue analysis: Modified Goodman
    % Rotating shaft assumption:
    % bending = alternating
    % torsion = mean
    % -------------------------
    Se_prime = 0.5*Sut;

    ka = 0.8;                         % machined surface
    kb = (Do*1000/7.62)^-0.107;       % size factor, Do in mm
    kc = 1.0;                         % load factor
    kd = 1.0;                         % temperature factor
    ke = 1.0;                         % reliability factor (approx)
    kf_misc = 1.0;                    % miscellaneous

    Se = Se_prime * ka * kb * kc * kd * ke * kf_misc;
    shaft.Se = Se;

    for k = 1:length(x)
        sigma_a = abs(shaft.sigma_b(k));         % fully reversed bending
        tau_m   = abs(shaft.tau_t(k));           % steady torsion
        sigma_m = sqrt(3)*tau_m;                 % equivalent mean stress

        if sigma_a > 1e-6 || sigma_m > 1e-6
            shaft.FoS_fatigue(k) = 1 / ( (sigma_a/Se) + (sigma_m/Sut) );
        end
    end

    % -------------------------
    % Summary values
    % -------------------------
    [shaft.V_max, iV] = max(abs(shaft.V));
    [shaft.M_max, iM] = max(abs(shaft.M));
    [shaft.T_max, iT] = max(abs(shaft.T));
    [shaft.sigma_vm_max, iVM] = max(abs(shaft.sigma_vm));

    shaft.x_Vmax = x(iV);
    shaft.x_Mmax = x(iM);
    shaft.x_Tmax = x(iT);
    shaft.x_sigma_vm_max = x(iVM);

    shaft.FoS_yield_min   = min(shaft.FoS_yield(isfinite(shaft.FoS_yield)));
    shaft.FoS_fatigue_min = min(shaft.FoS_fatigue(isfinite(shaft.FoS_fatigue)));

    % Max alternating bending stress for S-N point
    shaft.sigma_a_max = max(abs(shaft.sigma_b));
end

function print_results(shaft)
    fprintf('\n====================================================\n');
    fprintf('%s RESULTS\n', upper(shaft.name));
    fprintf('====================================================\n');

    for i = 1:length(shaft.supports)
        fprintf('Reaction at support x = %.4f m : %.2f N\n', shaft.supports(i), shaft.reactions(i));
    end

    fprintf('|V|max        = %.2f N at x = %.3f m\n', shaft.V_max, shaft.x_Vmax);
    fprintf('|M|max        = %.2f N*m at x = %.3f m\n', shaft.M_max, shaft.x_Mmax);
    fprintf('|T|max        = %.2f N*m at x = %.3f m\n', shaft.T_max, shaft.x_Tmax);
    fprintf('|sigma_vm|max = %.2f MPa at x = %.3f m\n', shaft.sigma_vm_max/1e6, shaft.x_sigma_vm_max);
    fprintf('Minimum FoS against yielding (von Mises) = %.3f\n', shaft.FoS_yield_min);
    fprintf('Minimum FoS against fatigue (Modified Goodman) = %.3f\n', shaft.FoS_fatigue_min);
end

function plot_shaft_results(shaft)
    x = shaft.x;

    figure('Name',[shaft.name ' - Shear'])
    plot(x, shaft.V, 'LineWidth', 2)
    grid on
    xlabel('x [m]')
    ylabel('Shear [N]')
    title([shaft.name ' - Shear Force Diagram'])

    figure('Name',[shaft.name ' - Moment'])
    plot(x, shaft.M, 'LineWidth', 2)
    grid on
    xlabel('x [m]')
    ylabel('Moment [N*m]')
    title([shaft.name ' - Bending Moment Diagram'])

    figure('Name',[shaft.name ' - Torque'])
    plot(x, shaft.T, 'LineWidth', 2)
    grid on
    xlabel('x [m]')
    ylabel('Torque [N*m]')
    title([shaft.name ' - Torque Diagram'])

    figure('Name',[shaft.name ' - Bending Stress'])
    plot(x, shaft.sigma_b/1e6, 'LineWidth', 2)
    grid on
    xlabel('x [m]')
    ylabel('\sigma_b [MPa]')
    title([shaft.name ' - Bending Stress'])

    figure('Name',[shaft.name ' - Torsional Shear'])
    plot(x, shaft.tau_t/1e6, 'LineWidth', 2)
    grid on
    xlabel('x [m]')
    ylabel('\tau [MPa]')
    title([shaft.name ' - Torsional Shear Stress'])

    figure('Name',[shaft.name ' - von Mises'])
    plot(x, shaft.sigma_vm/1e6, 'LineWidth', 2)
    grid on
    xlabel('x [m]')
    ylabel('\sigma_{vm} [MPa]')
    title([shaft.name ' - von Mises Stress'])

    FoS_y = shaft.FoS_yield;
    FoS_y(FoS_y > 10) = 10;
    figure('Name',[shaft.name ' - Yield FoS'])
    plot(x, FoS_y, 'LineWidth', 2)
    grid on
    xlabel('x [m]')
    ylabel('FoS against yielding')
    title([shaft.name ' - Yield FoS using von Mises'])
    ylim([0 10])

    FoS_f = shaft.FoS_fatigue;
    FoS_f(FoS_f > 10) = 10;
    figure('Name',[shaft.name ' - Fatigue FoS'])
    plot(x, FoS_f, 'LineWidth', 2)
    grid on
    xlabel('x [m]')
    ylabel('FoS against fatigue')
    title([shaft.name ' - Fatigue FoS using Modified Goodman'])
    ylim([0 10])
end

function plot_layout(shaft)
    figure('Name',[shaft.name ' - Layout'])
    hold on

    plot([0 shaft.L], [0 0], 'k', 'LineWidth', 6)

    % Supports in red
    for i = 1:length(shaft.supports)
        xs = shaft.supports(i);
        plot([xs xs], [-0.03 0.03], 'r', 'LineWidth', 3)
        text(xs, 0.055, sprintf('Support %d', i), 'HorizontalAlignment', 'center')
    end

    % Gears in blue
    for i = 1:length(shaft.gear_positions)
        xg = shaft.gear_positions(i);
        plot([xg xg], [-0.05 0.05], 'b', 'LineWidth', 3)
        text(xg, 0.080, sprintf('Gear %d', i), 'HorizontalAlignment', 'center')
    end

    xlabel('x [m]')
    title([shaft.name ' - Layout (Supports in Red, Gears in Blue)'])
    grid on
    axis([0 shaft.L -0.1 0.1])
    hold off
end

function plot_SN_single(shaft, Sut, Do)
    Se_prime = 0.5*Sut;

    ka = 0.8;
    kb = (Do*1000/7.62)^-0.107;
    kc = 1.0;
    kd = 1.0;
    ke = 1.0;
    kf_misc = 1.0;

    Se = Se_prime * ka * kb * kc * kd * ke * kf_misc;

    N = logspace(3, 7, 500);

    N1 = 1e3;
    S1 = 0.9*Sut;
    N2 = 1e6;
    S2 = Se;

    b = (log10(S2) - log10(S1)) / (log10(N2) - log10(N1));
    a = S1 / (N1^b);

    S = a * N.^b;

    N_ref = 1e6;
    sigma_a = shaft.sigma_a_max;

    figure('Name',[shaft.name ' - S-N'])
    loglog(N, S/1e6, 'LineWidth', 2)
    hold on
    yline(Se/1e6, '--', 'Endurance Limit')
    loglog(N_ref, sigma_a/1e6, 'ro', 'MarkerSize', 10, 'LineWidth', 2)

    grid on
    xlabel('Number of Cycles (N)')
    ylabel('Stress Amplitude [MPa]')
    title([shaft.name ' - S-N Diagram'])
    legend('S-N Curve', 'Endurance Limit', 'Operating Point', 'Location', 'southwest')
    hold off
end

function plot_SN_combined(mid, bot, top, Sut, Do)
    Se_prime = 0.5*Sut;

    ka = 0.8;
    kb = (Do*1000/7.62)^-0.107;
    kc = 1.0;
    kd = 1.0;
    ke = 1.0;
    kf_misc = 1.0;

    Se = Se_prime * ka * kb * kc * kd * ke * kf_misc;

    N = logspace(3, 7, 500);

    N1 = 1e3;
    S1 = 0.9*Sut;
    N2 = 1e6;
    S2 = Se;

    b = (log10(S2) - log10(S1)) / (log10(N2) - log10(N1));
    a = S1 / (N1^b);

    S = a * N.^b;

    N_ref = 1e6;

    figure('Name','Combined S-N Diagram')
    loglog(N, S/1e6, 'LineWidth', 2)
    hold on
    yline(Se/1e6, '--', 'Endurance Limit')

    loglog(N_ref, mid.sigma_a_max/1e6, 'ro', 'MarkerSize', 10, 'LineWidth', 2)
    loglog(N_ref, bot.sigma_a_max/1e6, 'go', 'MarkerSize', 10, 'LineWidth', 2)
    loglog(N_ref, top.sigma_a_max/1e6, 'bo', 'MarkerSize', 10, 'LineWidth', 2)

    grid on
    xlabel('Number of Cycles (N)')
    ylabel('Stress Amplitude [MPa]')
    title('S-N Diagram for All Shafts')
    legend('S-N Curve', 'Endurance Limit', ...
           mid.name, bot.name, top.name, ...
           'Location', 'southwest')
    hold off
end