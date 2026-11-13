function cathode_sizing()
% CATHODE_SIZING  Sizing calculations for the cathodes of the 150 W argon
% Hall effect thruster (MSc IRP, Cranfield University, 2026).
%
% This script repeats the calculations of the three spreadsheets used in
% the thesis:
%
%   Bare_W_Filament_Calculator.xlsx         -> filament_design
%   Langmuir_Sheet.xlsx                     -> collector_curve
%   Hollow_Cathode_Design_Sizing_Tool.xlsx  -> case_ozturk, case_ahmed,
%                                              case_montero, argon_design
%
% Run it by typing   cathode_sizing   in the MATLAB command window. It prints
% every result and saves two figures as PNG files in the current folder.
% Inputs are set at the top of each function so they can be changed without
% touching the equations.
%
% Y. Andiyappan, November 2026

bare = filament_design();

k = constants();
J_thor = 0.5 / (pi * 0.25e-3 * 0.14);
fprintf('\nW-ThO2, uniform 0.14 m wire: %.0f K\n', ...
    emission_temperature(J_thor, k.W_THO2(1), k.W_THO2(2), 1800));
fprintf('Pure W, uniform 0.14 m wire: %.0f K\n', ...
    emission_temperature(J_thor, k.PURE_W(1), k.PURE_W(2), 2500));

[G, plateaus] = collector_curve();

case_ozturk();
case_ahmed();
case_montero();
argon_design();

plot_collector(G, plateaus);
plot_filament(bare);
fprintf('\nFigures saved: collector_iv.png, filament_prediction.png\n');
end


% -------------------------------------------------------------------------
% Physical constants (CODATA 2018)
% -------------------------------------------------------------------------
function k = constants()
k.E     = 1.602176634e-19;      % elementary charge, C
k.KB    = 1.380649e-23;         % Boltzmann constant, J/K
k.ME    = 9.1093837e-31;        % electron mass, kg
k.EPS0  = 8.8541878128e-12;     % vacuum permittivity, F/m
k.SIGMA = 5.670374e-8;          % Stefan-Boltzmann constant, W/m^2/K^4
k.RGAS  = 8.314462618;          % gas constant, J/mol/K
k.NA    = 6.02214076e23;        % Avogadro number, 1/mol
k.GAMMA = 5 / 3;                % ratio of specific heats, monatomic gas
k.SCCM  = 101325e-6 / (k.RGAS * 273.15) / 60;   % mol/s per sccm

k.M_AR = 0.039948;              % argon molar mass, kg/mol
k.M_XE = 0.131293;              % xenon molar mass, kg/mol

% Emission constants [A_G in A/m^2/K^2, phi in eV]
k.PURE_W = [6.0e5, 4.54];       % Wolf (1995)
k.W_THO2 = [3.0e4, 2.63];       % thoriated tungsten
k.LAB6   = [2.9e5, 2.66];       % Lafferty (1951), as used by Goebel and Katz
k.BAO_W  = 1.2e6;               % Cronin (1981), phi = 1.67 + 2.82e-4 T
end


% -------------------------------------------------------------------------
% Basic relations used throughout
% -------------------------------------------------------------------------
function J = richardson(T, A_G, phi)
% Richardson-Dushman current density in A/m^2.
k = constants();
J = A_G .* T.^2 .* exp(-k.E * phi ./ (k.KB * T));
end

function T = emission_temperature(J, A_G, phi, T)
% Temperature that gives current density J. The equation has no closed form
% in T, so it is solved by Newton iteration on ln(J).
k = constants();
for i = 1:50
    f = log(A_G * T^2) - k.E * phi / (k.KB * T) - log(J);
    dfdT = 2 / T + k.E * phi / (k.KB * T^2);
    step = f / dfdT;
    T = T - step;
    if abs(step) < 1e-9
        break
    end
end
end

function [factor, dphi] = schottky_factor(T, V, gap)
% Increase in emission from the keeper field E = V/gap.
k = constants();
field = V / gap;
dphi = sqrt(k.E * field / (4 * pi * k.EPS0));       % in volts
factor = exp(k.E * dphi / (k.KB * T));
end

function I = child_langmuir(V, gap, area)
% Space-charge-limited current for a planar gap, in A.
k = constants();
I = (4 * k.EPS0 / 9) * sqrt(2 * k.E / k.ME) * V.^1.5 / gap^2 * area;
end

function P = choked_pressure(mdot, d_orifice, M, T_gas)
% Insert pressure for choked (continuum) flow through the orifice, Pa.
if nargin < 4
    T_gas = 1000;
end
k = constants();
area = pi * d_orifice^2 / 4;
flow = sqrt(k.GAMMA * M / (k.RGAS * T_gas)) * ...
    (2 / (k.GAMMA + 1))^((k.GAMMA + 1) / (2 * (k.GAMMA - 1)));
P = mdot / (area * flow);
end

function P = free_molecular_pressure(mdot, d_orifice, M, clausing, T_gas)
% Upper bound on insert pressure for free-molecular flow, Pa.
if nargin < 5
    T_gas = 1000;
end
k = constants();
area = pi * d_orifice^2 / 4;
P = mdot / (clausing * area) * sqrt(2 * pi * k.RGAS * T_gas / M);
end

function rho = rho_tungsten(T)
% Resistivity of tungsten, Cezairliyan and McClure (1971), ohm m.
rho = (-14.08 + 3.515e-2 * T) * 1e-8;
end

function p = pct(calc, ref)
% Percentage difference of a calculated value from a reference.
p = 100 * (calc - ref) / ref;
end


% -------------------------------------------------------------------------
% 1. Filament cathode (Chapter 4)
% -------------------------------------------------------------------------
function out = filament_design(A_G, phi, label)
k = constants();
if nargin < 3
    A_G = k.PURE_W(1);
    phi = k.PURE_W(2);
    label = 'pure W';
end

d = 0.25e-3;            % wire diameter, m
L_exposed = 0.14;       % exposed loop between the copper tube ends, m
L_emit = 0.11;          % effective emitting length, m (cooler ends emit little)
I_em = 0.50;            % required emission current, A
eps_2000 = 0.318;       % W emittance at 2000 K (Cezairliyan)
eps_3000 = 0.346;       % W emittance at 3000 K, linear between the two
rho_cold = 5.45e-8;     % resistivity at 293 K, ohm m
P_joint = 2.0;          % conduction into the two joints, W

A_emit = pi * d * L_emit;
A_rad = pi * d * L_exposed;
A_c = pi * d^2 / 4;
J = I_em / A_emit;
T = emission_temperature(J, A_G, phi, 2500);

R_hot = rho_tungsten(T) * L_exposed / A_c;
R_cold = rho_cold * L_exposed / A_c;

emis = eps_2000 + (eps_3000 - eps_2000) * (T - 2000) / 1000;
P_rad = emis * k.SIGMA * A_rad * T^4;
P_el = I_em * (phi + 2 * k.KB * T / k.E);       % electron cooling
P = P_rad + P_el + P_joint;
I_h = sqrt(P / R_hot);
V_h = I_h * R_hot;

fprintf('\n--- Filament (%s) ---\n', label);
fprintf('Required current density   %10.4g A/m2\n', J);
fprintf('Operating temperature      %10.1f K\n', T);
fprintf('R_hot / R_cold             %10.3f / %.4f ohm\n', R_hot, R_cold);
fprintf('Emittance at T             %10.3f\n', emis);
fprintf('Radiation / electron / joints  %.1f / %.2f / %.1f W\n', P_rad, P_el, P_joint);
fprintf('Heater power               %10.1f W\n', P);
fprintf('Heater current / voltage   %10.2f A / %.2f V\n', I_h, V_h);
fprintf('Current at a cold CV start %10.1f A\n', V_h / R_cold);

% Sweep in temperature for the figure
Ts = 2000:20:3000;
emis_s = eps_2000 + (eps_3000 - eps_2000) * (Ts - 2000) / 1000;
I_s = richardson(Ts, A_G, phi) * A_emit;
P_s = emis_s * k.SIGMA * A_rad .* Ts.^4 + I_s .* (phi + 2 * k.KB * Ts / k.E) + P_joint;
Ih_s = sqrt(P_s ./ (rho_tungsten(Ts) * L_exposed / A_c));

out.T = T;
out.P = P;
out.I_h = I_h;
out.V_h = V_h;
out.R_hot = R_hot;
out.R_cold = R_cold;
out.Ts = Ts;
out.I_em = I_s;
out.Ih = Ih_s;
end


% -------------------------------------------------------------------------
% 2. Collector of the cathode test rig (Chapter 6)
% -------------------------------------------------------------------------
function [G, plateaus] = collector_curve()
k = constants();
d = 0.25e-3;            % filament wire diameter, m
L_emit = 0.11;          % emitting length, m
r_a = 0.0259;           % inner radius of the copper collector, m
I_target = 0.50;

% heater current (A) and peak filament temperature (K), W-ThO2 filament
settings = [3.00 1704; 3.25 1765; 3.50 1819; 3.75 1866];

% Coaxial diode with beta^2 = 1:  I = G V^1.5
G = 8 * pi * k.EPS0 / 9 * sqrt(2 * k.E / k.ME) * L_emit / r_a;
V_target = (I_target / G)^(2 / 3);

fprintf('\n--- Collector (space-charge limit) ---\n');
fprintf('Perveance G                %10.4g A/V^1.5\n', G);
fprintf('Voltage for %.1f A          %10.1f V\n', I_target, V_target);
fprintf('Current at 300 V           %10.3f A\n', G * 300^1.5);

plateaus = zeros(size(settings, 1), 3);
for i = 1:size(settings, 1)
    I_h = settings(i, 1);
    T = settings(i, 2);
    I_sat = richardson(T, k.W_THO2(1), k.W_THO2(2)) * pi * d * L_emit;
    V_knee = (I_sat / G)^(2 / 3);
    plateaus(i, :) = [I_h, T, I_sat];
    fprintf('  %.2f A heater, %d K: plateau %.3f A from %.0f V\n', I_h, T, I_sat, V_knee);
end
end


% -------------------------------------------------------------------------
% 3. Hollow cathode: validation against published cathodes (Chapter 5)
% -------------------------------------------------------------------------
function case_ozturk()
% Ozturk et al. (2013), Bogazici LaB6 cathode. Used for the sensitivity to
% the work function, not as an exact match.
k = constants();
D_i = 2e-3;
L_i = 10e-3;
T = 1923;
I_pub = 3.0;
A = pi * D_i * L_i;
I = richardson(T, k.LAB6(1), 2.67) * A;
phi_fit = k.KB * T / k.E * log(k.LAB6(1) * T^2 * A / I_pub);
L_min = I_pub / richardson(T, k.LAB6(1), 2.67) / (pi * D_i);
% implied emissivity of the heater assembly from 107 W at 2024 K
emis = 107 / (k.SIGMA * pi * 8e-3 * 15e-3 * 2024^4);

fprintf('\n--- Case P1: Ozturk 2013 ---\n');
fprintf('Emission at 1923 K         %10.2f A   (published 3 A)\n', I);
fprintf('phi that gives 3 A         %10.3f eV\n', phi_fit);
fprintf('Minimum insert length      %10.2f mm  (built 10 mm)\n', L_min * 1e3);
fprintf('Implied heater emissivity  %10.3f\n', emis);
end

function case_ahmed()
% Ahmed et al. (2024), Surrey BaO-W cathode. Anchors the orifice model.
k = constants();
D_i = 2e-3;
L_i = 14e-3;
d_o = 1e-3;
T = 1373;
phi = 1.67 + 2.82e-4 * T;
J = k.BAO_W * T^2 * exp(-k.E * phi / (k.KB * T));
L_req = 2.5 / 5.5e4 / (pi * D_i);
P1 = choked_pressure(0.1e-6, d_o, k.M_XE);
P5 = choked_pressure(0.5e-6, d_o, k.M_XE);
I_seed = child_langmuir(250, 2e-3, pi * D_i * L_i);

fprintf('\n--- Case P2: Ahmed 2024 ---\n');
fprintf('J at 1373 K                %10.2f A/cm2  (published 5-6)\n', J / 1e4);
fprintf('Required insert length     %10.2f mm  (built 14 mm)\n', L_req * 1e3);
fprintf('Pressure at 0.1 mg/s       %10.1f Pa  (%+.1f %% vs 43 Pa)\n', P1, pct(P1, 43));
fprintf('Pressure at 0.5 mg/s       %10.1f Pa  (%+.1f %% vs 215 Pa)\n', P5, pct(P5, 215));
fprintf('Child-Langmuir at 250 V    %10.3f A\n', I_seed);
end

function case_montero()
% Montero Barriga et al. (2024), UC3M sub-ampere LaB6 cathode.
% Cases 1 and 3 of the thesis.
k = constants();
D_i = 1.3e-3;
L_i = 10e-3;
gap = 2.25e-3;
T_ign = 1575;
V_k = 300;
A = pi * D_i * L_i;

I_th = richardson(T_ign, k.LAB6(1), k.LAB6(2)) * A;
boost = schottky_factor(T_ign, V_k, gap);
I_seed = I_th * boost;
I_cl = child_langmuir(V_k, gap, A);
A_min = 0.1 * A / I_seed;
I_1773 = richardson(1773, k.LAB6(1), k.LAB6(2)) * A;

% heater: tungsten wire wound on a BN sleeve
rho = 4.8e-8 * (1 + 4.8297e-3 * 1585 + 1.663e-6 * 1585^2);
L_w = pi * (13.4e-3 + 0.5e-3) * 9;
R = rho * L_w / (pi * 0.25e-3^2);
I_h = sqrt(220 / R);
V_h = sqrt(220 * R);

fprintf('\n--- Case P3: Montero Barriga 2024 ---\n');
fprintf('Emission area              %10.2f mm2  (published 40.84)\n', A * 1e6);
fprintf('Seed current, 300 V        %10.4f A  (%+.1f %%)\n', I_seed, pct(I_seed, 0.1));
fprintf('Child-Langmuir, 300 V      %10.4f A  (%+.1f %%)\n', I_cl, pct(I_cl, 0.1));
fprintf('Minimum area for 0.1 A     %10.2f mm2\n', A_min * 1e6);
fprintf('Emission at 1773 K         %10.3f A  (rated 1 A)\n', I_1773);
fprintf('Heater R / I / V           %.3f ohm / %.2f A / %.2f V  (%+.1f %%, published 1.23 / 13.36 / 16.47)\n', ...
    R, I_h, V_h, pct(R, 1.23));
fprintf('Measured ignition (Case 3): heater 267-305 W against 220 W predicted, factor %.2f-%.2f\n', ...
    267 / 220, 305 / 220);
end


% -------------------------------------------------------------------------
% 4. Hollow cathode: argon design (Chapter 5)
% -------------------------------------------------------------------------
function argon_design()
k = constants();

% thruster current balance
P_d = 150;
V_d = 300;
mdot_anode = 0.175e-6;                      % kg/s
I_d = P_d / V_d;
I_b_max = mdot_anode * k.NA * k.E / k.M_AR;
I_b = 0.8 * I_b_max;                        % 80 % mass utilisation
I_em = I_d;                                 % emission equals discharge current
I_design = 1.5 * I_em;
I_seed = 0.1 * I_em;

% insert
D_i = 2e-3;
L_i = 10e-3;
T_op = 1800;
A = pi * D_i * L_i;
J = richardson(T_op, k.LAB6(1), k.LAB6(2));
L_min = I_design / J / (pi * D_i);
T_nom = emission_temperature(I_design / A, k.LAB6(1), k.LAB6(2), 1700);
T_dirty = emission_temperature(I_design / A, k.LAB6(1), k.LAB6(2) + 0.15, 1700);
T_seed = emission_temperature(I_seed / A, k.LAB6(1), k.LAB6(2), 1500);
T_seed_dirty = emission_temperature(I_seed / A, k.LAB6(1), k.LAB6(2) + 0.15, 1500);
I_dirty = richardson(T_op, k.LAB6(1), k.LAB6(2) + 0.15) * A;

% keeper
V_k = 400;
gap = 2.25e-3;
[~, dphi] = schottky_factor(T_seed, V_k, gap);
I_cl = child_langmuir(V_k, gap, A);

% orifice and flow
d_o = 0.5e-3;
flow = 2.0;                                 % sccm
mdot = flow * k.SCCM * k.M_AR;
P_cont = choked_pressure(mdot, d_o, k.M_AR);
P_fm = free_molecular_pressure(mdot, d_o, k.M_AR, 0.357);
mfp = k.KB * 1000 / (sqrt(2) * pi * (3.6e-10)^2 * P_cont);
Kn = mfp / d_o;

% heater: tantalum wire wound on a BN sleeve
eps_ta = 0.1991 + 4.687e-5 * T_op;          % Cezairliyan (1971)
OD = 10e-3;
L_h = 10e-3;
P_rad = eps_ta * k.SIGMA * pi * OD * L_h * T_op^4;
P_ss = 1.5 * P_rad;                         % 50 % added for conduction
P_heritage = 79 * (1800 / 1773)^4;          % Kokal (2021), scaled by T^4
P_design = 130;
d_w = 0.25e-3;
L_w = 0.24;
A_w = pi * d_w^2 / 4;
R_hot = 7.72e-7 * L_w / A_w;                % rho_Ta at 1950 K
R_cold = 1.4e-7 * L_w / A_w;
I_h = sqrt(P_design / R_hot);
V_h = sqrt(P_design * R_hot);
turns = L_w / (pi * (OD + d_w));

fprintf('\n--- Argon hollow cathode design ---\n');
fprintf('Discharge current          %10.2f A\n', I_d);
fprintf('Beam current (max / 80 %%)  %.3f / %.3f A\n', I_b_max, I_b);
fprintf('Design / seed current      %.2f / %.2f A\n', I_design, I_seed);
fprintf('J at 1800 K                %10.4g A/m2\n', J);
fprintf('Minimum insert length      %10.2f mm  (10 mm chosen)\n', L_min * 1e3);
fprintf('T for design current       %.0f K nominal, %.0f K contaminated\n', T_nom, T_dirty);
fprintf('T for seed current         %.0f K nominal, %.0f K contaminated\n', T_seed, T_seed_dirty);
fprintf('Emission at 1800 K         %.2f A nominal, %.2f A contaminated\n', J * A, I_dirty);
fprintf('Schottky lowering at 400 V %10.4f eV\n', dphi);
fprintf('Child-Langmuir at 400 V    %10.3f A\n', I_cl);
fprintf('Cathode flow               %10.3e kg/s\n', mdot);
fprintf('Insert pressure            %.0f to %.0f Pa  (%.1f to %.1f Torr), Kn = %.2f\n', ...
    P_cont, P_fm, P_cont / 133.322, P_fm / 133.322, Kn);
fprintf('Heater: radiation %.1f W, steady state %.1f W, heritage %.1f W\n', ...
    P_rad, P_ss, P_heritage);
fprintf('Design power               %10.0f W  (x%.2f analytical, x%.2f heritage)\n', ...
    P_design, P_design / P_ss, P_design / P_heritage);
fprintf('Supply with transient      %10.0f W\n', P_design * 305 / 220);
fprintf('R_hot / R_cold             %.2f / %.3f ohm\n', R_hot, R_cold);
fprintf('Heater current / voltage   %.2f A / %.2f V, %.1f turns\n', I_h, V_h, turns);
fprintf('Wire heat flux             %10.3g W/m2\n', P_design / (pi * d_w * L_w));
end


% -------------------------------------------------------------------------
% Figures (serif font, thin black lines, no titles)
% -------------------------------------------------------------------------
function style_axes(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 9, 'LineWidth', 0.8, ...
    'TickDir', 'in', 'Box', 'on', 'XColor', 'k', 'YColor', 'k');
end

function plot_collector(G, plateaus)
V = linspace(0, 700, 1401);
I_sc = G * V.^1.5;

fig = figure('Units', 'inches', 'Position', [1 1 4.8 3.4], 'Color', 'w', ...
    'PaperPositionMode', 'auto');
ax = axes(fig);
hold(ax, 'on');
h = plot(ax, V, I_sc, 'k-', 'LineWidth', 1.0);
names = {'Space-charge limit'};
styles = {'--', '-.', ':', '-'};
shades = {'k', 'k', 'k', [0.4 0.4 0.4]};
for i = 1:size(plateaus, 1)
    I_h = plateaus(i, 1);
    T = plateaus(i, 2);
    I_sat = plateaus(i, 3);
    h(end + 1) = plot(ax, V, min(I_sc, I_sat), 'LineStyle', styles{i}, ...
        'Color', shades{i}, 'LineWidth', 1.0);
    names{end + 1} = sprintf('%.2f A heater (%d K)', I_h, T);
end
plot(ax, [0 700], [0.5 0.5], 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
xlim(ax, [0 700]);
ylim(ax, [0 0.8]);
xlabel(ax, 'Collector voltage, V_c (V)');
ylabel(ax, 'Collector current, I_c (A)');
legend(ax, h, names, 'Box', 'off', 'Location', 'northwest', 'FontSize', 8);
style_axes(ax);
print(fig, 'collector_iv.png', '-dpng', '-r300');
close(fig);
end

function plot_filament(bare)
fig = figure('Units', 'inches', 'Position', [1 1 4.8 3.4], 'Color', 'w', ...
    'PaperPositionMode', 'auto');
ax = axes(fig);
yyaxis(ax, 'left');
h1 = semilogy(ax, bare.Ts, bare.I_em, 'k-', 'LineWidth', 1.0);
ylim(ax, [1e-3 10]);
ylabel(ax, 'Emission current, I_{em} (A)');
yyaxis(ax, 'right');
h2 = plot(ax, bare.Ts, bare.Ih, 'k--', 'LineWidth', 1.0);
ylim(ax, [4 9.5]);
ylabel(ax, 'Heating current, I_h (A)');
xlim(ax, [2000 3000]);
xlabel(ax, 'Filament temperature, T (K)');
legend(ax, [h1 h2], {'Emission current', 'Heating current'}, 'Box', 'off', ...
    'Location', 'northwest', 'FontSize', 8);
style_axes(ax);
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
print(fig, 'filament_prediction.png', '-dpng', '-r300');
close(fig);
end
