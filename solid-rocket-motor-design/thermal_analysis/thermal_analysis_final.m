%% thermal_analysis_final.m — Stark-100k Cooling Jacket
%  POLIMI Space Propulsion 2025-26 | Group: I Maggici
%
% =========================================================================
% README
% =========================================================================
% PURPOSE:
%   Steady-state 1-D thermal analysis of the water-cooled jacket for the
%   Stark-100k solid rocket motor (100 kN thrust, non-aluminised AP/HTPB,
%   BATES grain, initial chamber pressure 63.77 bar).
%   The code verifies that the cooling jacket design keeps all wall layers
%   within operational limits and checks four safety criteria:
%     (1) cavitation avoidance in the coolant channels
%     (2) onset of nucleate boiling (ONB) — Bergles-Rohsenow criterion
%     (3) critical heat flux (CHF) — Zuber + subcooling correction
%     (4) mechanical stress on the CuCrZr liner — Von Mises criterion
%
% CODE STRUCTURE:
%   S1   — Pure scalar inputs (no formulas, no function calls)
%   S2   — Derived scalar quantities
%   S3   — Axial grid and conical nozzle profile
%   S4   — Gas thermodynamic and transport properties (NASA CEA interpolation)
%   S5   — Cooling channel geometry (variable width along axis)
%   S6   — External feedline pressure losses
%   S7   — 1-D axial marching thermal solver
%   S8   — Post-processing: safety checks and feed system sizing
%            S8a: jacket hydraulics and cavitation check
%            S8b: ONB check (Bergles-Rohsenow 1964)
%            S8c: CHF check (Zuber 1958 + Griffith subcooling)
%            S8d: mechanical stress on CuCrZr liner (Von Mises)
%            S8e: N2 blow-down sizing (adiabatic model, conservative)
%   S8.5 — Parametric sweep: TBC thickness vs T_wg and T_TBC
%   S8.6 — Sensitivity analysis: T_wg vs coolant mass flow rate
%   S8.7 — Sensitivity analysis: T_wg vs number of channels Nch
%   S8.8 — TBC vs no-TBC comparison (proves TBC necessity)
%   S9   — Figures (Figs 1–12)
%
% KEY MODELS:
%   - Bartz (1957) + iterative sigma factor          [hot-gas side h_g]
%   - Gnielinski (1976) + Petukhov friction factor   [coolant side h_c]
%   - Shah & London (1978) AR correction             [deep-channel Nu penalty]
%   - IAPWS-IF97 Region 1                            [liquid water properties]
%   - Wagner-Pruss (1993)                            [saturation curve]
%   - Bergles-Rohsenow (1964)                        [ONB threshold]
%   - Zuber (1958) + Griffith                        [pool CHF + subcooling]
%   - Von Mises criterion                            [liner stress]
%
% REQUIREMENTS:
%   MATLAB 2025b / Octave 11.1.x — no additional toolboxes required.
%
% INPUT SOURCE TAGS:
%   [ASSIGN] = from the assignment specification
%   [CEA]    = from NASA CEA2 (AP/HTPB, O/F=5.32, shifting-EQ to throat,
%              frozen after throat)
%   [LIT]    = from peer-reviewed literature (reference cited in comment)
%   [DESIGN] = engineering design choice, justified in the report
% =========================================================================

clear; clc; close all;

%% S1 — INPUTS  (scalars only; no formulas, no function calls)

% --- Combustion gas properties — NASA CEA, AP/HTPB non-aluminised, O/F=5.32
% Thermodynamic properties (gamma, T, p, Cp, M) from shifting-equilibrium
% composition in the convergent section, frozen from the throat onward.
% Transport properties (mu, Pr) from the frozen-reactions branch in both
% sections, as required for boundary-layer heat-transfer calculations
% (Sutton & Biblarz 9th ed. Sec. 8.3): in the thermal BL the chemical
% relaxation time >> flow time, so the mixture behaves as frozen.
T0      = 2852.29;      % [K]        stagnation temperature              [CEA]
Mmol    = 23.683;       % [kg/kmol]  mean molar mass                     [CEA]
p0      = 63.77e5;      % [Pa]       initial chamber stagnation pressure [DESIGN]
cst     = 1542.0;       % [m/s]      characteristic velocity c*          [CEA]
cp0_fr  = 1899.1;       % [J/(kgK)]  frozen Cp at chamber                [CEA]
Pr0_fr  = 0.6268;       % [-]        frozen Pr at chamber                [CEA]
mu0_fr  = 9.1736e-5;    % [Pa*s]     frozen viscosity at chamber         [CEA]

% --- Engine (assignment) -------------------------------------------------
F       = 100e3;        % [N]        sea-level thrust                 [ASSIGN]
Isp     = 240.4;        % [s]        sea-level specific impulse       [ASSIGN]
Itot    = 2.5e6;        % [N*s]      total impulse target             [ASSIGN]

% --- Nozzle geometry -----------------------------------------------------
AR_in       = 2.0;      % [-]        area ratio at jacket inlet       [ASSIGN]
AR_out      = 2.0;      % [-]        area ratio at jacket outlet      [ASSIGN]
a_conv_deg  = 45.0;     % [deg]      convergent half-angle            [DESIGN]
a_div_deg   = 15.0;     % [deg]      divergent half-angle             [DESIGN]
RcRt        = 1.5;      % [-]        throat curvature ratio Rc/Rt     [DESIGN]

% --- Thermal protection layers -------------------------------------------
% TBC top coat: YSZ (Yttria-Stabilised Zirconia)
use_TBC      = true;
k_TBC        = 1.5;            % [W/(mK)]  YSZ conductivity              [LIT: Marple JTST 2007]
T_lim_TBC    = 1400+273.15;    % [K]       YSZ phase-stability limit      [LIT: Chen 2025]
t_TBC        = 72e-6;          % [m]       YSZ thickness                  [DESIGN]

% Bond coat: NiCoCrAlY (ensures adhesion and CTE compatibility)
use_bond    = true;
t_bond      = 0.035e-3;        % [m]       bond coat thickness            [LIT: Marple 2007]
T_lim_bond  = 1050+273.15;     % [K]       NiCoCrAlY operating limit
k_bond      = 20.0;            % [W/(mK)]  NiCoCrAlY conductivity         [LIT: Clyne & Leevers]

% Structural liner: CuCrZr (precipitation-hardened copper alloy)
t_wall      = 1e-3;            % [m]       liner thickness                [DESIGN]
k_wall      = 320.0;           % [W/(mK)]  CuCrZr conductivity ~400 C    [LIT: ITER MPH]
T_lim_wall  = 773;             % [K]       CuCrZr over-ageing limit (500 C) [LIT: ITER MPH]

% --- Cooling channel geometry (channel width w_ch computed axially in S5) ---
% Variable-width strategy: fin thickness is held constant; channel width
% varies with local nozzle circumference (Strategy A, Sec. 3.3.2 of report).
h_ch        = 1.5e-3;          % [m]       channel height (fin depth)     [DESIGN]
Nch         = 90;              % [-]       number of channels             [DESIGN]
t_fin       = 1.5e-3;          % [m]       fin thickness (constant)       [DESIGN]

% --- Coolant: water, single pass; properties updated node-by-node (IAPWS-IF97)
T_in        = 18+273.15;       % [K]       inlet temperature (18 C)       [ASSIGN]
mdot_w      = 7.6;             % [kg/s]    total coolant mass flow rate   [DESIGN]
rho_w_ref   = 998.2;           % [kg/m^3]  reference density at T_in      [IAPWS]

% --- External hydraulic circuit ------------------------------------------
K_inlet     = 0.5;    % [-]   manifold-to-channel entry (sharp contraction)
K_outlet    = 1.0;    % [-]   channel-to-manifold exit (sudden expansion)
K_valve     = 0.1;    % [-]   ball valve (fully open)
K_orifice   = 5.0;    % [-]   calibration orifice
K_elbow     = 0.9;    % [-]   four 90-degree elbows combined
K_filter    = 2.0;    % [-]   cartridge filter
K_ext_sum   = K_valve + K_orifice + K_elbow + K_filter;
L_ext       = 2.0;    % [m]   feedline length
D_ext       = 25e-3;  % [m]   feedline inner diameter
eps_ext     = 1.5e-6; % [m]   pipe roughness (smooth stainless steel)
SF_cav      = 1.4;    % [-]   cavitation safety factor

% Dump tank: vented vessel at the jacket outlet — sets the lower bound
% on coolant static pressure inside the marching solver.
p_dump      = 1.5e5;  % [Pa]  effective dump-tank back-pressure (1.5 bar) [DESIGN]

% --- CuCrZr mechanical properties [ITER MPH] -----------------------------
lambda_E    = 17e-6;   % [1/K]  thermal expansion coefficient
E_wall      = 128e9;   % [Pa]   Young's modulus
nu_wall     = 0.34;    % [-]    Poisson's ratio
sigma_y0    = 450e6;   % [Pa]   yield strength at 20 C
dSY_dT      = -0.44e6; % [Pa/K] yield-strength temperature slope (20–400 C, linear)
SF_mec      = 1.2;     % [-]    mechanical safety factor

% --- Gas-band thermal radiation ------------------------------------------
% Non-aluminised AP/HTPB: no condensed-phase particles; radiation governed
% by CO2/H2O spectral bands only. Conservative upper bound (S&B Tab. 8-3).
eps_g       = 0.10;         % [-]          effective gas emissivity      [LIT: S&B Tab.8-3]
sigma_SB    = 5.6704e-8;    % [W/(m^2K^4)] Stefan-Boltzmann constant

% --- CHF safety factor ---------------------------------------------------
% SF = 1.5 is appropriate for forced-convection single-phase water cooling
% (subcooled throughout: ONB L1 and L2 both PASS).
% Reference: Huzel & Hwang NASA SP-125 Sec. 4.3; Celata et al. (1994).
SF_CHF_req  = 1.5;     % [-]   minimum required CHF safety factor        [DESIGN]

% --- Solver control parameters -------------------------------------------
N           = 200;     % [-]   number of axial nodes
N_sig       = 15;      % [-]   maximum fixed-point iterations per node
tol_sig     = 0.5;     % [K]   convergence tolerance on T_hg and T_wc

% --- NASA CEA gas-property tables (jacket region: AR = 1.0 to 2.0) ------
% Thermodynamics: shifting-equilibrium composition (convergent, x<=0) and
% frozen composition (divergent, x>0).  Transport properties (mu, Pr) from
% the frozen-reactions branch throughout, as appropriate for BL calculations.
% Two separate AR grids match the CEA output stations for each branch.
%
% Convergent branch — subsonic, AR from throat (1.0) to jacket inlet (2.0)
AR_tab_eq = [1.0,    1.05,   1.1,    1.2,    1.5,    2.0];    % [-] area ratio

T_eq    = [2598.28, 2694.19, 2724.21, 2757.98, 2801.34, 2826.29]; % [K]   static T
p_eq    = [37.704,  46.920,  50.222,  54.203,  59.761,  63.207]*1e5; % [Pa] static p
gam_eq  = [1.2076,  1.2029,  1.2014,  1.1998,  1.1977,  1.1965];  % [-]  gamma (EQ)
Cp_eq   = [1875.1,  1884.7,  1887.5,  1890.7,  1894.6,  1896.8];  % [J/(kgK)] Cp frozen
mu_eq   = [8.5925,  8.8139,  8.8827,  8.9598,  9.0584,  9.1149]*1e-5; % [Pa·s] mu frozen
Pr_eq   = [0.6306,  0.6292,  0.6288,  0.6283,  0.6276,  0.6272];  % [-]  Pr frozen
M_eq_tab= [1.000,   0.781,   0.701,   0.600,   0.439,   0.313];   % [-]  Mach (subsonic)

% Divergent branch — supersonic, AR from throat (1.0) to jacket outlet (2.0)
AR_tab_fr = [1.0,    1.05,   1.15,   1.5,    2.0];             % [-] area ratio

T_fr    = [2598.28, 2470.88, 2363.69, 2145.11, 1958.71];         % [K]   static T
p_fr    = [37.704,  28.080,  21.804,  12.758,  7.8401]*1e5;      % [Pa]  static p
gam_fr  = [1.2076,  1.2138,  1.2187,  1.2274,  1.2333];          % [-]   gamma (frozen)
Cp_fr   = [1875.1,  1861.5,  1849.0,  1820.4,  1792.1];          % [J/(kgK)] Cp frozen
mu_fr   = [8.5925,  8.2944,  8.0398,  7.5091,  7.0429]*1e-5;     % [Pa·s] mu frozen
Pr_fr   = [0.6306,  0.6323,  0.6335,  0.6356,  0.6365];          % [-]   Pr frozen
M_fr_tab= [1.000,   1.243,   1.427,   1.774,   2.063];           % [-]   Mach (supersonic)

% --- N2-pressurised regulated feed system --------------------------------
% A dome-loaded pressure regulator maintains p_tank constant at the water
% vessel until the N2 bottle pressure drops to the lock-up value
% p_N2_final = p_tank + dp_regulator. The N2 vessel is sized for the full
% test window using an ADIABATIC blow-down model (see S8e for details).
dt_extra    = 4.0;     % [s]      extra cooling time (2 s pre + 2 s post-shutdown)
p_tank      = 52e5;    % [Pa]     regulated water-tank pressure          [DESIGN]
p_N2_init   = 300e5;   % [Pa]     N2 vessel initial pressure             [DESIGN]
dp_regulator= 5e5;     % [Pa]     dome-loaded regulator lock-up drop     [DESIGN]
T_gas       = 293.15;  % [K]      N2 initial temperature (ambient)
M_N2        = 28.0134; % [kg/kmol] N2 molar mass
gamma_N2    = 1.40;    % [-]      N2 isentropic exponent (diatomic ideal gas)
feed_model  = 'adiabatic'; % ['adiabatic'|'isothermal'] — adiabatic is design case

%% S2 — DERIVED SCALAR QUANTITIES

g0      = 9.80665;                    % [m/s^2]      standard gravity
Runiv   = 8314.5;                     % [J/(kmol·K)] universal gas constant
Rgas    = Runiv / Mmol;               % [J/(kgK)]    specific gas constant

mdot_prop = F / (Isp * g0);           % [kg/s]  propellant mass flow rate
t_burn    = Itot / F;                 % [s]     nominal burn time

% Throat area from the initial chamber pressure p0 and characteristic
% velocity c* (both from the ballistic design).  Gives At = 102.6 cm²
% consistent with Table 4.1 of the report.
At  = mdot_prop * cst / p0;          % [m^2]  throat area  →  102.6 cm²
Rt  = sqrt(At / pi);                  % [m]    throat radius
Dt  = 2 * Rt;                         % [m]    throat diameter
Rc  = RcRt * Rt;                      % [m]    throat curvature radius (for Bartz)

a_conv = a_conv_deg * pi/180;                % [rad]
a_div  = a_div_deg  * pi/180;               % [rad]
R_in   = Rt * sqrt(AR_in);                   % [m]     jacket inlet radius
R_out  = Rt * sqrt(AR_out);                  % [m]     jacket outlet radius

% Conical-only lengths (no throat arc correction at this stage)
L_conv = (R_in  - Rt) / tan(a_conv);
L_div  = (R_out - Rt) / tan(a_div);
L_tot  = L_conv + L_div;

% Thermal resistances per unit area [m^2 K/W]
R_bond = (t_bond / k_bond) * use_bond;
R_wall = t_wall  / k_wall;

mdot_ch = mdot_w / Nch;                      % [kg/s]  per-channel mass flow

[rho_in, ~, mu_in, ~] = iapws_if97_r1(T_in, p_tank);  % only rho and mu needed (S6)

fprintf('=== Stark-100k Thermal Analysis ===\n');
fprintf('Dt=%.1f mm  At=%.2f cm^2  mdot=%.2f kg/s  t_b=%.1f s\n', ...
        Dt*1e3, At*1e4, mdot_prop, t_burn);

%% S3 — NOZZLE GEOMETRY  (conical profile, throat at x=0)
% The axial grid runs from x=-L_conv (jacket inlet, subsonic AR=2) to
% x=+L_div (jacket outlet, supersonic AR=2).  Throat at x=0.
% n_conv nodes are allocated proportionally to the convergent length.

n_conv = max(3, round(N * L_conv / L_tot));
n_div  = N - n_conv + 1;
x_conv = linspace(-L_conv, 0,     n_conv);
x_div  = linspace(0,       L_div, n_div);
x      = [x_conv, x_div(2:end)];
N      = length(x);

% Conical wall radius r(x)
r = zeros(1, N);
for i = 1:N
    if x(i) <= 0
        r(i) = Rt + abs(x(i)) * tan(a_conv);
    else
        r(i) = Rt + x(i) * tan(a_div);
    end
end

AR      = (r / Rt).^2;                       % [-]  local area ratio A(x)/At
dx_vec  = diff(x);
ds      = sqrt(dx_vec.^2 + diff(r).^2);      % [m]  arc-length element (wall slant)
r_avg   = (r(1:end-1) + r(2:end)) / 2;       % [m]  node-average radius
C_wall  = 2*pi*r;                             % [m]  wall circumference
[~, i_th] = min(abs(x));                      % [-]  throat node index
x_mm      = x * 1e3;                          % [mm] for plots

%% S4 — GAS PROPERTIES  (interpolation on local area ratio AR(x))
% Convergent (x<=0): shifting-equilibrium tables (AR_tab_eq).
% Divergent (x>0):  frozen-composition tables (AR_tab_fr).
% Adiabatic wall temperature T_aw from S&B Eq. 8-16:
%   T_aw = T0*(1 + r_rec*(gam-1)/2*M^2) / (1 + (gam-1)/2*M^2)
% Turbulent recovery factor: r_rec = Pr^(1/3).

gam_loc = zeros(1,N); Cp_loc = zeros(1,N);
Pr_loc  = zeros(1,N); mu_loc = zeros(1,N);
T_s     = zeros(1,N); p_s    = zeros(1,N);
M       = zeros(1,N);

for i = 1:N
    if x(i) <= 0
        % Convergent: subsonic stations (AR_tab_eq)
        gam_loc(i) = interp1(AR_tab_eq,gam_eq, AR(i),'linear','extrap');
        Cp_loc(i)  = interp1(AR_tab_eq,Cp_eq,  AR(i),'linear','extrap');
        Pr_loc(i)  = interp1(AR_tab_eq,Pr_eq,  AR(i),'linear','extrap');
        mu_loc(i)  = interp1(AR_tab_eq,mu_eq,  AR(i),'linear','extrap');
        T_s(i)     = interp1(AR_tab_eq,T_eq,   AR(i),'linear','extrap');
        p_s(i)     = interp1(AR_tab_eq,p_eq,   AR(i),'linear','extrap');
        M(i)       = interp1(AR_tab_eq,M_eq_tab,AR(i),'linear','extrap');
    else
        % Divergent: supersonic stations (AR_tab_fr)
        gam_loc(i) = interp1(AR_tab_fr,gam_fr, AR(i),'linear','extrap');
        Cp_loc(i)  = interp1(AR_tab_fr,Cp_fr,  AR(i),'linear','extrap');
        Pr_loc(i)  = interp1(AR_tab_fr,Pr_fr,  AR(i),'linear','extrap');
        mu_loc(i)  = interp1(AR_tab_fr,mu_fr,  AR(i),'linear','extrap');
        T_s(i)     = interp1(AR_tab_fr,T_fr,   AR(i),'linear','extrap');
        p_s(i)     = interp1(AR_tab_fr,p_fr,   AR(i),'linear','extrap');
        M(i)       = interp1(AR_tab_fr,M_fr_tab,AR(i),'linear','extrap');
    end
end

rho_s  = p_s * Mmol ./ (Runiv * T_s);        % [kg/m^3]  ideal-gas density
a_s    = sqrt(gam_loc * Rgas .* T_s);         % [m/s]     speed of sound
r_rec  = Pr_loc.^(1/3);                       % [-]       turbulent recovery factor
xi     = 1 + (gam_loc-1)/2 .* M.^2;           % [-]       stagnation ratio
T_aw   = T0 .* (1 + r_rec.*(gam_loc-1)/2.*M.^2) ./ xi;  % [K]  adiabatic wall T

%% S5 — CHANNEL GEOMETRY  (Strategy: variable channel width)
% The channel pitch tracks the local nozzle circumference:
%   pitch(x) = C_wall(x) / Nch
%   w_ch(x)  = pitch(x) - t_fin   (fin thickness constant for structural uniformity)
% The wetted perimeter is 3-sided: P_wet = w_ch + 2*h_ch (open-top channel).

pitch_arr  = C_wall / Nch;
w_ch_arr   = pitch_arr - t_fin;

if w_ch_arr(i_th) <= 0
    error('Channels do not fit at throat: pitch=%.2f mm < t_fin=%.2f mm.', ...
          pitch_arr(i_th)*1e3, t_fin*1e3);
end
if w_ch_arr(i_th) < 0.5e-3
    warning('w_ch at throat=%.2f mm < 0.5 mm CNC manufacturing limit.', ...
            w_ch_arr(i_th)*1e3);
end

A_ch_arr   = w_ch_arr * h_ch;                 % [m^2]  channel cross-section
P_wet_arr  = w_ch_arr + 2*h_ch;               % [m]    3-sided wetted perimeter
Dh_arr     = 4 * A_ch_arr ./ P_wet_arr;       % [m]    hydraulic diameter
AR_ch_arr  = h_ch ./ w_ch_arr;                % [-]    aspect ratio h/w

fprintf('\nChannel geometry | Inlet  | Throat | Outlet\n');
fprintf('w_ch [mm]        | %6.2f | %6.2f | %6.2f\n', ...
        w_ch_arr(1)*1e3, w_ch_arr(i_th)*1e3, w_ch_arr(end)*1e3);
fprintf('Dh   [mm]        | %6.2f | %6.2f | %6.2f\n', ...
        Dh_arr(1)*1e3, Dh_arr(i_th)*1e3, Dh_arr(end)*1e3);

%% S6 — FEEDLINE PRESSURE BUDGET
% Computes the coolant static pressure available at the jacket inlet after
% all external line losses (distributed friction + concentrated losses).
% This value is then passed to the marching solver as the inlet boundary
% condition.  Cavitation is NOT checked here: it is checked inside the
% jacket (S8a), where both the local static pressure and the local coolant
% temperature are known at every node.

A_ext = pi * D_ext^2 / 4;                           % [m^2]  feedline cross-section
v_ext = mdot_w / (rho_in * A_ext);                   % [m/s]  feedline velocity
Re_ext = rho_in * v_ext * D_ext / mu_in;             % [-]    feedline Reynolds number

% Haaland friction factor (explicit, valid for all Re > 4000)
f_ext = (-1.8*log10((eps_ext/(3.7*D_ext))^1.11 + 6.9/Re_ext))^(-2);

dp_ext_dist = f_ext * (L_ext/D_ext) * 0.5 * rho_in * v_ext^2;  % distributed losses
dp_ext_conc = 0.5 * rho_in * v_ext^2 * K_ext_sum;              % concentrated losses
dp_feedline = dp_ext_dist + dp_ext_conc;                        % total feedline drop
p_inlet     = p_tank - dp_feedline;                             % jacket inlet pressure

fprintf('\nFeedline pressure budget:\n');
fprintf('  v_ext = %.2f m/s  Re = %.0f  f = %.4f\n', v_ext, Re_ext, f_ext);
fprintf('  dp_distributed = %.3f bar  dp_concentrated = %.3f bar\n', ...
        dp_ext_dist/1e5, dp_ext_conc/1e5);
fprintf('  dp_feedline = %.3f bar  =>  p_inlet = %.2f bar\n', ...
        dp_feedline/1e5, p_inlet/1e5);

%% S7 — 1-D AXIAL MARCHING THERMAL SOLVER
% At each node the solver:
%  (0) Updates local static pressure by subtracting friction losses from
%      the previous element (forward Euler, Petukhov friction factor).
%  (1) Evaluates IAPWS-IF97 water properties at (T_w(i), p_node).
%  (2) Computes flow quantities (v, Re).
%  (3) Runs a fixed-point loop on the hot-gas-side surface temperature
%      T_hg and the cold-wall temperature T_wc until convergence:
%       - Bartz sigma (S&B Eq.8-23) uses T_hg
%       - Gnielinski Nu + Petukhov f_D uses T_wc
%       - Fin efficiency (adiabatic tip, Incropera Eq.3.86/3.91)
%       - Gas-band radiation absorbed into T_aw_eff = T_aw + q_rad/h_g
%       - Series thermal resistance network: gas|TBC|bond|wall|coolant
%  (4) Marches coolant enthalpy: T_w(i+1) from finite-difference balance.

% Bartz (1957) prefactor: function of throat geometry and chamber conditions
% only — constant along the nozzle axis.  Local h_g(x) is then obtained as
% h_g = Bpf * AR(x)^(-0.9) * sigma(T_hg), with sigma evaluated iteratively
% inside the fixed-point loop at each node.
Bpf = (0.026/Dt^0.2) * (Dt/Rc)^0.1 ...
    * (mu0_fr^0.2 * cp0_fr / Pr0_fr^0.6) ...
    * (p0/cst)^0.8;        % [W/(m^2 K)]  Bartz equation prefactor (S&B Eq. 8-23)

R_TBC_arr = (t_TBC / k_TBC) * use_TBC * ones(1, N);

% --- Main solver call ----------------------------------------------------
[T_w, T_wc, T_wg, T_bond, T_TBC, q_flux, h_g, sigma, ...
 rho_w_arr, cp_w_arr, mu_w_arr, k_w_arr, Pr_w_arr, ...
 v_w_arr, Re_arr, Nu_arr, h_c_arr, hc_eff_arr, eta_fin_arr, eta_o_arr, ...
 ~, ~, p_static, q_rad_arr, ~] = ...
    marching_solver(N, N_sig, tol_sig, T_aw, T_s, T0, xi, AR, Bpf, ...
                    R_TBC_arr, R_bond, R_wall, ...
                    A_ch_arr, P_wet_arr, Dh_arr, r, AR_ch_arr, ...
                    mdot_w, mdot_ch, Nch, h_ch, t_fin, k_wall, ...
                    r_avg, ds, T_in, p_inlet, p_dump, K_inlet, K_outlet, ...
                    eps_g, sigma_SB);

[q_max_val, i_qmax] = max(q_flux);

fprintf('\nThermal solver results:\n');
fprintf('  Bpf = %.2f kW/(m^2K)\n', Bpf/1e3);
fprintf('  T_wg  [C] max = %.1f  (limit %.0f C  %s)\n', ...
        max(T_wg)-273.15, T_lim_wall-273, margin_str(max(T_wg),T_lim_wall));
fprintf('  T_TBC [C] max = %.1f  (limit %.0f C  %s)\n', ...
        max(T_TBC)-273.15, T_lim_TBC-273, margin_str(max(T_TBC),T_lim_TBC));
fprintf('  T_bond[C] max = %.1f  (limit %.0f C  %s)\n', ...
        max(T_bond)-273.15, T_lim_bond-273, margin_str(max(T_bond),T_lim_bond));
fprintf('  q_max = %.2f MW/m^2 @ x=%.1f mm\n', q_max_val/1e6, x(i_qmax)*1e3);
fprintf('  q_rad @ throat = %.3f MW/m^2 (%.1f%% of q_conv)\n', ...
        q_rad_arr(i_th)/1e6, 100*q_rad_arr(i_th)/max(q_flux(i_th),1e-10));
fprintf('  Coolant DeltaT = %.1f K\n', T_w(end)-T_in);

%% S8 — POST-PROCESSING

% -------------------------------------------------------------------------
% 8a — Jacket cavitation check (node-by-node)
% -------------------------------------------------------------------------
% Cavitation occurs when the local static pressure drops below the local
% vapour pressure.  It is meaningful to check this only where both the
% static pressure AND the coolant temperature are known simultaneously,
% i.e., inside the jacket where the marching solver couples hydraulics and
% energy.
%
% The feedline merely establishes p_inlet = p_tank - dp_feedline (S6).
% Inside the jacket, p_static(i) decreases due to friction while T_w(i)
% increases due to heat absorption — both effects reduce the cavitation
% margin.  The critical point is wherever min(p_static / (SF_cav*p_vap))
% occurs along the jacket axis.
%
% Criterion:  p_static(x) > SF_cav * p_vap(T_w(x))  at every node.

p_vap_arr = arrayfun(@p_sat_wagner, T_w);               % p_vap at local T_w(x)
cav_margin = p_static ./ (SF_cav * p_vap_arr);          % cavitation margin, all nodes
[cav_min_val, i_cav_min] = min(cav_margin);

fprintf('\nJacket hydraulics and cavitation:\n');
fprintf('  p_static [bar] inlet/throat/outlet = %.2f / %.2f / %.2f\n', ...
        p_static(1)/1e5, p_static(i_th)/1e5, p_static(end)/1e5);
fprintf('  Cavitation margin min = %.2f @ x=%.1f mm  %s\n', ...
        cav_min_val, x(i_cav_min)*1e3, ternary(cav_min_val>1,'PASS','FAIL'));

% -------------------------------------------------------------------------
% 8b — ONB check (Bergles-Rohsenow 1964)
% -------------------------------------------------------------------------
% Criterion L1: T_sat(p_static) > T_wc  -> wall is subcooled -> safe
% Criterion L2: actual superheat dT_sup < dT_ONB_req  -> nucleation
%               cannot be activated (Bergles-Rohsenow threshold).

T_sat_arr    = arrayfun(@T_sat_wagner, p_static);
ONB_margin_1 = T_sat_arr - T_wc;               % [K]  >0 = safe (L1)
% Surface tension and latent heat evaluated at T_sat (saturation state),
% as required by the Bergles-Rohsenow and Zuber correlations.
% Using T_w (bulk, 18-44 C) here would give sigma ~2x higher than the
% correct value at T_sat (~263 C at 50 bar), causing a ~40% error on
% dT_ONB_req and ~19% on q_CHF_Zuber.
sigma_LG     = 0.0756 - 1.42e-4*(T_sat_arr - 273.15); % [N/m] at T_sat
h_fg         = 2.501e6 - 2430*(T_sat_arr - 273.15); % [J/kg] latent heat at T_sat
rho_V        = p_static ./ (461.5 * T_sat_arr); % [kg/m^3] vapour density @ p_static
dT_sup       = T_wc - T_sat_arr;               % [K]   wall superheat

q_ONB_BR = (k_w_arr .* h_fg .* rho_V) ./ (8*sigma_LG .* T_sat_arr) .* dT_sup.^2;
q_ONB_BR(dT_sup <= 0) = Inf;                   % subcooled wall: ONB impossible

dT_ONB_req = sqrt(max(0, 8*sigma_LG.*T_sat_arr.*q_flux./(k_w_arr.*h_fg.*rho_V)));

% Flow regime is fully characterised by L1 and L2 criteria above.
% (regime vector removed — not used in current figures)

[ONB1_min, i_ONB1] = min(ONB_margin_1);
ratio_BR = q_flux ./ q_ONB_BR;
ratio_BR_max = max(ratio_BR(isfinite(ratio_BR)));
if isempty(ratio_BR_max), ratio_BR_max = 0; end

fprintf('\nONB check:\n');
fprintf('  L1 (T_sat-T_wc) min = %+.1f K @ x=%.1f mm  %s\n', ...
        ONB1_min, x(i_ONB1)*1e3, ternary(ONB1_min>0,'PASS','NB possible'));
fprintf('  L2 (B-R ratio)  max = %.3g  %s\n', ...
        ratio_BR_max, ternary(ratio_BR_max<1,'PASS (single-phase)','SNB regime'));

% -------------------------------------------------------------------------
% 8c — CHF (Critical Heat Flux): Zuber + Griffith subcooling correction
% -------------------------------------------------------------------------
% (a) q_CHF_pool = Zuber + h_c*DeltaT_sub  [conservative bound, informative only]
% (b) q_CHF_fc   = q_CHF_pool * E_fc        [forced-conv. realistic, design check]
% Rationale: the Zuber pool-boiling limit is physically inapplicable to forced-
% convection single-phase flow (ONB L1 and L2 both PASS).  The FC bound (Celata
% enhancement) is the appropriate criterion for the actual operating regime.
% Heat flux reference: channel-wall area (q_flux converted via area ratio).

q_CHF_Zuber = (pi/24) .* h_fg .* rho_V .* ...
              (sigma_LG.*g0.*(rho_w_arr-rho_V)./rho_V.^2).^0.25;
DeltaT_sub  = max(0, T_sat_arr - T_w);         % [K]  bulk subcooling
q_CHF_pool  = q_CHF_Zuber + h_c_arr .* DeltaT_sub;  % conservative bound

% Celata-type forced-convection enhancement factor (capped at 3)
Pe_w  = Re_arr .* Pr_w_arr;
Ja_sub= cp_w_arr .* DeltaT_sub ./ max(h_fg,1);
E_fc  = min(3.0, 1.0 + 0.25*Ja_sub.*(Pe_w/1e5).^0.5);
q_CHF_fc = q_CHF_pool .* E_fc;

% Convert gas-side heat flux to channel-wall basis
q_flux_ch  = q_flux .* (2*pi*r) ./ (Nch .* P_wet_arr);

SF_CHF_arr = q_CHF_fc ./ max(q_flux_ch, 1e-6);    % SF based on FC (design criterion)
[SF_CHF_min_val, i_SF_CHF] = min(SF_CHF_arr);

fprintf('\nCHF check:\n');
fprintf('  q_CHF_fc @ throat = %.2f MW/m^2  q_flux_ch = %.4f MW/m^2\n', ...
        q_CHF_fc(i_th)/1e6, q_flux_ch(i_th)/1e6);
fprintf('  SF_CHF (FC) min = %.1f @ x=%.1f mm  %s (req > %.1f)\n', ...
        SF_CHF_min_val, x(i_SF_CHF)*1e3, ...
        ternary(SF_CHF_min_val>SF_CHF_req,'PASS','CHECK'), SF_CHF_req);
fprintf('  SF_CHF (pool, conservative) min = %.1f  (informative only)\n', ...
        min(q_CHF_pool./max(q_flux_ch,1e-6)));

% -------------------------------------------------------------------------
% 8d — Mechanical stress on the CuCrZr liner
% -------------------------------------------------------------------------
% The liner carries ONLY the differential pressure |p_coolant - p_gas|;
% the global chamber pressure is reacted by the outer backing shell
% (Huzel & Hwang, NASA SP-125, Sec. 8.4).
%
% Two stress contributions on the liner:
%   sigma_mem  = |Dp| * r / t_wall          (membrane / hoop stress)
%   sigma_bend = 0.5 * |Dp| * (w_ch/t_wall)^2  (plate bending at fin roots)
% Thermal biaxial stress: sigma_th = lambda_E * E * dT_wall / (1 - nu)
% Von Mises criterion with temperature-degraded yield strength.

Dp_liner   = p_static - p_s;
abs_Dp     = abs(Dp_liner);
sigma_mem  = abs_Dp .* r / t_wall;
sigma_bend = 0.5 .* abs_Dp .* (w_ch_arr ./ t_wall).^2;
sigma_p    = max(sigma_mem, sigma_bend);
dT_wall    = T_wg - T_wc;
sigma_th   = lambda_E * E_wall * dT_wall / (1 - nu_wall);
sigma_VM   = sqrt(sigma_p.^2 + sigma_th.^2 - sigma_p.*sigma_th);

T_wall_avg = 0.5*(T_wg + T_wc);
sigma_y    = max(50e6, sigma_y0 + dSY_dT*(T_wall_avg - 293.15));
% Floor on sigma_VM prevents division by zero at nodes with negligible
% pressure differential and temperature gradient (e.g. very first node).
mec_margin = sigma_y ./ (SF_mec * max(sigma_VM, 1e3));
[mec_min_val, i_mec_min] = min(mec_margin);

fprintf('\nLiner mechanical stress:\n');
fprintf('  sigma_mem max = %.1f MPa  sigma_th max = %.1f MPa  sigma_VM max = %.1f MPa\n', ...
        max(sigma_mem)/1e6, max(sigma_th)/1e6, max(sigma_VM)/1e6);
fprintf('  Mech. margin min = %.2f @ x=%.1f mm  %s\n', ...
        mec_min_val, x(i_mec_min)*1e3, ternary(mec_min_val>1,'PASS','FAIL'));

% -------------------------------------------------------------------------
% 8e — Feed system sizing: water tank and N2 pressurant bottle
% -------------------------------------------------------------------------
% The N2 blow-down bottle must keep p_N2 >= p_N2_final = p_tank + dp_reg
% throughout the full test window (burn + pre/post margins).
%
% Two limiting thermodynamic models bracket the real expansion:
%
%   ISOTHERMAL (pV = const): the gas stays at constant temperature T_gas,
%   exchanging heat freely with the vessel walls.  This is the OPTIMISTIC
%   model: pressure decays slowly, so a smaller bottle suffices.
%   Equation: V_N2_iso = V_water * p_final / (p_init - p_final)
%
%   ADIABATIC (pV^gamma = const): no heat exchange occurs during expansion.
%   The gas cools as it expands, so pressure decays faster than isothermal.
%   A LARGER bottle is therefore needed to meet the same p_final.
%   Equation: V_N2_adi = V_water / [(p_init/p_final)^(1/gamma) - 1]
%   with V_N2_adi > V_N2_iso since (p_init/p_final)^(1/gamma) < p_init/p_final.
%
%   DESIGN CHOICE — adiabatic: for a test of only ~29 s the N2 has little
%   time to exchange heat with the steel walls, making the adiabatic model
%   both more physically representative and more conservative (larger bottle).
%
% NOTE on Fig 12: each curve is plotted using its own correctly-sized bottle.
%   The adiabatic curve therefore starts with MORE gas and appears to decay
%   MORE SLOWLY in the plot, while the isothermal curve (smaller bottle)
%   decays FASTER.  Both reach p_N2_final at the end of the test by design.
%   The gap between the two curves illustrates the conservatism of the
%   adiabatic model over the entire test window.

t_burn_total = t_burn + dt_extra;             % [s]   total cooling system operation time
m_water      = mdot_w * t_burn_total;         % [kg]  total water mass required
V_water      = m_water / rho_w_ref;           % [m^3] total water volume

p_N2_final   = p_tank + dp_regulator;         % [Pa]  N2 lock-up pressure
R_N2         = Runiv / M_N2;                  % [J/(kgK)] N2 specific gas constant

% Bottle volumes — both models computed for comparison
V_N2_iso = V_water * p_N2_final / (p_N2_init - p_N2_final);       % isothermal (lower bound)
V_N2_adi = V_water / ((p_N2_init/p_N2_final)^(1/gamma_N2) - 1);   % adiabatic  (design)
m_N2_iso = p_N2_init * V_N2_iso / (R_N2 * T_gas);
m_N2_adi = p_N2_init * V_N2_adi / (R_N2 * T_gas);

% Select design values
switch lower(feed_model)
    case 'isothermal', V_N2 = V_N2_iso; m_N2 = m_N2_iso;
    case 'adiabatic',  V_N2 = V_N2_adi; m_N2 = m_N2_adi;
end

fprintf('\nFeed system sizing:\n');
fprintf('  t_op = %.1f s  |  m_water = %.1f kg  |  V_water = %.1f L\n', ...
        t_burn_total, m_water, V_water*1e3);
fprintf('  Adiabatic (design): V_N2 = %.1f L  m_N2 = %.2f kg\n', V_N2_adi*1e3, m_N2_adi);
fprintf('  Isothermal  (ref) : V_N2 = %.1f L  m_N2 = %.2f kg\n', V_N2_iso*1e3, m_N2_iso);
fprintf('  Lock-up pressure p_N2_final = %.0f bar\n', p_N2_final/1e5);

% Pre-compute N2 pressure evolution for Fig 12
% Each model uses its own sized volume, so both curves reach p_N2_final at t_burn_total.
t_vec       = linspace(0, t_burn_total, 300);
V_expelled  = mdot_w * t_vec / rho_w_ref;                              % [m^3] water expelled
p_N2_iso_t  = p_N2_init * V_N2_iso ./ (V_N2_iso + V_expelled);        % isothermal
p_N2_adi_t  = p_N2_init * (V_N2_adi ./ (V_N2_adi + V_expelled)).^gamma_N2; % adiabatic

%% S8.5 — TBC THICKNESS PARAMETRIC SWEEP
% Sweeps t_TBC from 0 to 250 µm and records T_wg_max and T_TBC_max.
% Justifies the design thickness as the minimum that keeps T_TBC < T_lim_YSZ.

fprintf('\nTBC sweep running...\n');
t_TBC_vec   = linspace(0, 250e-6, 30);
T_wg_max_v  = zeros(size(t_TBC_vec));
T_TBC_max_v = zeros(size(t_TBC_vec));

for kt = 1:length(t_TBC_vec)
    R_TBC_sw = (t_TBC_vec(kt)/k_TBC)*use_TBC*ones(1,N);
    [~,~,T_wg_sw,~,T_TBC_sw] = marching_solver( ...
        N, N_sig, tol_sig, T_aw, T_s, T0, xi, AR, Bpf, ...
        R_TBC_sw, R_bond, R_wall, A_ch_arr, P_wet_arr, Dh_arr, r, AR_ch_arr, ...
        mdot_w, mdot_ch, Nch, h_ch, t_fin, k_wall, ...
        r_avg, ds, T_in, p_inlet, p_dump, K_inlet, K_outlet, eps_g, sigma_SB);
    T_wg_max_v(kt)  = max(T_wg_sw);
    T_TBC_max_v(kt) = max(T_TBC_sw);
end

feas_mask = T_TBC_max_v < T_lim_TBC;
if any(feas_mask)
    [~,ik] = min(T_wg_max_v(feas_mask));
    idx_f  = find(feas_mask); i_opt = idx_f(ik);
    fprintf('  Optimum: t_TBC=%.0f um -> T_wg=%.0f C  T_TBC=%.0f C\n', ...
            t_TBC_vec(i_opt)*1e6, T_wg_max_v(i_opt)-273.15, T_TBC_max_v(i_opt)-273.15);
end
fprintf('  Design:  t_TBC=%.0f um -> T_wg=%.0f C  T_TBC=%.0f C\n', ...
        t_TBC*1e6, max(T_wg)-273.15, max(T_TBC)-273.15);

%% S8.6 — COOLANT MASS FLOW SENSITIVITY
% Identifies the minimum mdot_w that keeps T_wg below T_lim_wall.

fprintf('\nmdot_w sweep running...\n');
mdot_vec = linspace(2.0, 12.0, 20);
T_wg_max_mdot = zeros(size(mdot_vec));

for km = 1:length(mdot_vec)
    mdot_ch_k = mdot_vec(km)/Nch;
    [~,~,T_wg_md] = marching_solver( ...
        N, N_sig, tol_sig, T_aw, T_s, T0, xi, AR, Bpf, ...
        R_TBC_arr, R_bond, R_wall, A_ch_arr, P_wet_arr, Dh_arr, r, AR_ch_arr, ...
        mdot_vec(km), mdot_ch_k, Nch, h_ch, t_fin, k_wall, ...
        r_avg, ds, T_in, p_inlet, p_dump, K_inlet, K_outlet, eps_g, sigma_SB);
    T_wg_max_mdot(km) = max(T_wg_md);
end
mdot_min_safe = mdot_vec(find(T_wg_max_mdot < T_lim_wall,1,'first'));
if isempty(mdot_min_safe), mdot_min_safe = NaN; end
fprintf('  Min safe mdot = %.1f kg/s  |  design %.1f kg/s (margin %.0f K)\n', ...
        mdot_min_safe, mdot_w, T_lim_wall-max(T_wg));

%% S8.7 — NUMBER OF CHANNELS Nch SENSITIVITY
% Sweeps Nch with t_fin and h_ch fixed (Geometry rebuilt each time).
% Reports T_wg_max, maximum coolant velocity, and total pressure drop.

fprintf('\nNch sweep running...\n');
Nch_vec = [60, 90, 120, 150, 180];
T_wg_max_Nch = zeros(size(Nch_vec));
v_max_Nch    = zeros(size(Nch_vec));
dp_tot_Nch   = zeros(size(Nch_vec));

for kN = 1:length(Nch_vec)
    Nch_k     = Nch_vec(kN);
    pitch_k   = 2*pi*r / Nch_k;
    w_ch_k    = pitch_k - t_fin;
    A_ch_k    = w_ch_k * h_ch;
    P_wet_k   = w_ch_k + 2*h_ch;
    Dh_k      = 4*A_ch_k ./ P_wet_k;
    AR_ch_k   = h_ch ./ w_ch_k;
    mdot_ch_k = mdot_w / Nch_k;
    [~,~,T_wg_Nk,~,~,~,~,~,~,~,~,~,~,v_Nk,~,~,~,~,~,~,~,~,p_Nk] = ...
        marching_solver(N, N_sig, tol_sig, T_aw, T_s, T0, xi, AR, Bpf, ...
            R_TBC_arr, R_bond, R_wall, A_ch_k, P_wet_k, Dh_k, r, AR_ch_k, ...
            mdot_w, mdot_ch_k, Nch_k, h_ch, t_fin, k_wall, ...
            r_avg, ds, T_in, p_inlet, p_dump, K_inlet, K_outlet, eps_g, sigma_SB);
    T_wg_max_Nch(kN) = max(T_wg_Nk);
    v_max_Nch(kN)    = max(v_Nk);
    dp_tot_Nch(kN)   = p_inlet - min(p_Nk);
end
fprintf('  Nch | T_wg_max[C] | v_max[m/s] | dp[bar]\n');
for kN = 1:length(Nch_vec)
    fprintf('  %3d |   %6.1f    |   %5.1f   |  %5.2f\n', ...
        Nch_vec(kN), T_wg_max_Nch(kN)-273.15, v_max_Nch(kN), dp_tot_Nch(kN)/1e5);
end

%% S8.8 — TBC vs NO-TBC COMPARISON
% Re-runs the solver without TBC and bond coat to quantify the thermal
% benefit and prove that TBC is mandatory to stay below T_lim_wall.

fprintf('\nNo-TBC run...\n');
R_TBC_noTBC = zeros(1,N);
[~,~,T_wg_noTBC,~,~,q_noTBC] = marching_solver( ...
    N, N_sig, tol_sig, T_aw, T_s, T0, xi, AR, Bpf, ...
    R_TBC_noTBC, 0, R_wall, A_ch_arr, P_wet_arr, Dh_arr, r, AR_ch_arr, ...
    mdot_w, mdot_ch, Nch, h_ch, t_fin, k_wall, ...
    r_avg, ds, T_in, p_inlet, p_dump, K_inlet, K_outlet, eps_g, sigma_SB);

dT_TBC = max(T_wg_noTBC) - max(T_wg);
dq_TBC = max(q_noTBC)    - max(q_flux);
fprintf('  With TBC:    T_wg_max=%.1f C  q_max=%.2f MW/m^2\n', ...
        max(T_wg)-273.15, max(q_flux)/1e6);
fprintf('  Without TBC: T_wg_max=%.1f C  q_max=%.2f MW/m^2\n', ...
        max(T_wg_noTBC)-273.15, max(q_noTBC)/1e6);
fprintf('  TBC benefit: dT=%.1f K  dq=%.2f MW/m^2 (%.0f%%)\n', ...
        dT_TBC, dq_TBC/1e6, 100*dq_TBC/max(q_noTBC));
if max(T_wg_noTBC) > T_lim_wall
    fprintf('  --> Without TBC: T_wg=%.0f C EXCEEDS CuCrZr limit=%.0f C. TBC IS REQUIRED.\n', ...
            max(T_wg_noTBC)-273.15, T_lim_wall-273.15);
end

%% S9 — FIGURES

% Colour palette (dark-background friendly)
c_red    = [0.95 0.20 0.10];
c_orange = [0.95 0.55 0.05];
c_yellow = [0.90 0.85 0.10];
c_green  = [0.15 0.90 0.30];
c_cyan   = [0.05 0.85 0.95];
c_blue   = [0.20 0.55 0.95];
c_violet = [0.75 0.20 0.95];
c_white  = [1.00 1.00 1.00];
c_black  = [0.00 0.00 0.00];

% Fine grid for smooth interpolated plots
x_fine    = linspace(x(1), x(end), 5*N);
x_fine_mm = x_fine * 1e3;
T_w_f    = interp1(x, T_w,    x_fine,'linear');
T_wc_f   = interp1(x, T_wc,   x_fine,'linear');
T_wg_f   = interp1(x, T_wg,   x_fine,'linear');
T_bond_f = interp1(x, T_bond, x_fine,'linear');
T_TBC_f  = interp1(x, T_TBC,  x_fine,'linear');
q_f      = interp1(x, q_flux, x_fine,'linear');
q_rad_f  = interp1(x, q_rad_arr, x_fine,'linear');

% -------------------------------------------------------------------
% Fig 1 — Nozzle profile, area ratio, variable channel width
% -------------------------------------------------------------------
figure('Name','Fig 1 — Nozzle & Channel Geometry','Position',[40 40 1400 420]);
subplot(1,3,1);
    plot(x_mm, r*1e3, 'Color',c_blue,'LineWidth',2.5); hold on;
    plot(x_mm,-r*1e3, 'Color',c_blue,'LineWidth',2.5);
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; axis equal;
    xlabel('x [mm]'); ylabel('r [mm]'); title('Nozzle profile r(x)');
subplot(1,3,2);
    plot(x_mm, AR,'Color',c_violet,'LineWidth',2.5); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    yline(1,'Color',c_white,'LineStyle',':','LineWidth',1.0);
    hold off; grid on; box on;
    xlabel('x [mm]'); ylabel('A/A_t [-]'); title('Area ratio AR(x)');
subplot(1,3,3);
    plot(x_mm, w_ch_arr*1e3,'Color',c_orange,'LineWidth',2.5,'DisplayName','w_{ch}(x)'); hold on;
    yline(t_fin*1e3,'Color',c_white,'LineStyle','--','LineWidth',1.2, ...
          'DisplayName',sprintf('t_{fin}=%.1f mm',t_fin*1e3));
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on;
    xlabel('x [mm]'); ylabel('w_{ch} [mm]'); title('Channel width w_{ch}(x)');
    legend('Location','best');

% -------------------------------------------------------------------
% Fig 2 — CEA gas thermodynamic state (8 subplots)
% -------------------------------------------------------------------
figure('Name','Fig 2 — Gas Thermodynamic State','Position',[40 480 1600 700]);
subplot(2,4,1);
    plot(x_mm,M,'Color',c_red,'LineWidth',2); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    yline(1,'Color',c_white,'LineStyle',':','LineWidth',1.0);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('M [-]'); title('Mach Number');
subplot(2,4,2);
    plot(x_mm,T_s-273.15,'Color',c_orange,'LineWidth',2); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('T_s [C]'); title('Static Temperature T_s');
subplot(2,4,3);
    plot(x_mm,p_s/1e5,'Color',c_yellow,'LineWidth',2); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('p_s [bar]'); title('Static Pressure p_s');
subplot(2,4,4);
    plot(x_mm,rho_s,'Color',c_green,'LineWidth',2); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('\rho_s [kg/m^3]'); title('Gas Density \rho_s');
subplot(2,4,5);
    plot(x_mm,a_s,'Color',c_cyan,'LineWidth',2); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('a_s [m/s]'); title('Speed of Sound');
subplot(2,4,6);
    plot(x_mm,Cp_loc,'Color',c_blue,'LineWidth',2); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('C_p [J/(kgK)]'); title('Specific Heat C_p');
subplot(2,4,7);
    plot(x_mm,Pr_loc,'Color',c_violet,'LineWidth',2); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('Pr [-]'); title('Prandtl Number');
subplot(2,4,8);
    plot(x_mm,mu_loc*1e6,'Color',c_red,'LineWidth',2); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('\mu [\muPa\cdots]'); title('Dynamic Viscosity');

% -------------------------------------------------------------------
% Fig 3 — HTC profiles and fin efficiency
% -------------------------------------------------------------------
figure('Name','Fig 3 — HTC and Fin Efficiency','Position',[80 80 1300 500]);
subplot(1,3,1);
    semilogy(x_mm,h_g/1e3,        'Color',c_red, 'LineWidth',2.5,'DisplayName','h_g'); hold on;
    semilogy(x_mm,hc_eff_arr/1e3, 'Color',c_blue,'LineWidth',2.5,'DisplayName','h_{c,eff}');
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('HTC [kW/(m^2K)]');
    title('HTC Comparison h_g vs h_{c,eff}'); legend('Location','best');
subplot(1,3,2);
    plot(x_mm,eta_fin_arr,'Color',c_green, 'LineWidth',2.5,'DisplayName','\eta_{fin}'); hold on;
    plot(x_mm,eta_o_arr,  'Color',c_orange,'LineWidth',2.5,'DisplayName','\eta_o');
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('\eta [-]');
    title('Fin Efficiencies \eta_{fin}, \eta_o'); legend('Location','best');
subplot(1,3,3);
    plot(x_mm,sigma,'Color',c_violet,'LineWidth',2.5); hold on;
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2);
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('\sigma [-]');
    title('Bartz correction factor \sigma');

% -------------------------------------------------------------------
% Fig 4 — Wall heat flux (convective + radiative)
% -------------------------------------------------------------------
figure('Name','Fig 4 — Wall Heat Flux','Position',[100 100 900 450]);
plot(x_fine_mm,(q_f+q_rad_f)/1e6,'Color',c_red,   'LineWidth',2.5,'DisplayName','q_{tot}'); hold on;
plot(x_fine_mm,q_f/1e6,           'Color',c_orange,'LineWidth',2.0,'DisplayName','q_{conv}');
plot(x_fine_mm,q_rad_f/1e6,       'Color',c_yellow,'LineWidth',2.0,'DisplayName','q_{rad}');
scatter(x(i_qmax)*1e3,q_max_val/1e6,100,'v','filled','MarkerFaceColor',c_yellow, ...
        'HandleVisibility','off');
text(x(i_qmax)*1e3+2,q_max_val/1e6,sprintf(' %.2f MW/m^2',q_max_val/1e6), ...
     'Color',c_yellow,'FontWeight','bold');
xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
hold off; grid on; box on;
xlabel('x [mm]'); ylabel('q'''' [MW/m^2]');
title('Wall heat flux — convective + radiative'); legend('Location','best');

% -------------------------------------------------------------------
% Fig 5 — Temperature stack (all layers)
% -------------------------------------------------------------------
figure('Name','Fig 5 — Temperature Stack','Position',[120 120 1100 520]);
hold on;
    plot(x_mm,T_aw-273.15,'Color',c_white,'LineStyle','--','LineWidth',1.5,'DisplayName','T_{aw}');
    if use_TBC
        plot(x_fine_mm,T_TBC_f-273.15,'Color',c_red,   'LineWidth',2.5,'DisplayName','T_{TBC} (YSZ)');
        yline(T_lim_TBC-273.15,'Color',c_red,'LineStyle',':','LineWidth',1.5,'HandleVisibility','off');
    end
    if use_bond
        plot(x_fine_mm,T_bond_f-273.15,'Color',c_orange,'LineWidth',2.0,'DisplayName','T_{bond}');
        yline(T_lim_bond-273.15,'Color',c_orange,'LineStyle',':','LineWidth',1.5,'HandleVisibility','off');
    end
    plot(x_fine_mm,T_wg_f-273.15,'Color',c_violet,'LineWidth',2.5,'DisplayName','T_{wg} (CuCrZr hot)');
    plot(x_fine_mm,T_wc_f-273.15,'Color',c_blue,  'LineWidth',2.5,'DisplayName','T_{wc} (CuCrZr cold)');
    plot(x_fine_mm,T_w_f-273.15, 'Color',c_cyan,  'LineWidth',2.5,'DisplayName','T_w (bulk)');
    yline(T_lim_wall-273.15,'Color',c_violet,'LineStyle',':','LineWidth',1.5,'HandleVisibility','off');
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
hold off;
xlabel('x [mm]'); ylabel('T [C]'); title('Temperature Stack — All Layers');
legend('Location','best'); grid on; box on;

% -------------------------------------------------------------------
% Fig 6 — Safety margins (cavitation, ONB L1, ONB L2, mechanical)
% -------------------------------------------------------------------
figure('Name','Fig 6 — Safety Margins','Position',[140 140 1400 680]);

subplot(2,2,1);  % Cavitation
    semilogy(x_mm,p_static/1e5,        'Color',c_blue,  'LineWidth',2.5,'DisplayName','p_{static}'); hold on;
    semilogy(x_mm,SF_cav*p_vap_arr/1e5,'Color',c_red,   'LineStyle','--','LineWidth',2.0,'DisplayName','SF\cdot p_{vap}');
    semilogy(x_mm,p_vap_arr/1e5,       'Color',c_orange,'LineStyle',':','LineWidth',1.5,'DisplayName','p_{vap}');
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('p [bar] (log)');
    title(sprintf('Cavitation — margin min=%.2f',min(cav_margin))); legend('Location','east');

subplot(2,2,2);  % ONB L1: bulk subcooling
    plot(x_mm,ONB_margin_1,'Color',c_green,'LineWidth',2.5,'DisplayName','T_{sat}-T_{wc}'); hold on;
    yline(0,'Color',c_red,'LineStyle','--','LineWidth',1.5,'DisplayName','Sat. limit');
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('T_{sat}-T_{wc} [K]');
    title(sprintf('ONB L1 — subcooling min=%.0f K',ONB1_min)); legend('Location','best');

subplot(2,2,3);  % ONB L2: Bergles-Rohsenow superheat threshold
    plot(x_mm,dT_ONB_req,'Color',c_yellow,'LineWidth',2.5,'DisplayName','\DeltaT_{ONB,req}'); hold on;
    plot(x_mm,dT_sup,    'Color',c_cyan,  'LineWidth',2.5,'DisplayName','\DeltaT_{sup}');
    yline(0,'Color',c_red,'LineStyle','--','LineWidth',1.5,'HandleVisibility','off');
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('\DeltaT [K]');
    title(sprintf('ONB L2 (Bergles-Rohsenow)'));
    legend('Location','best');

subplot(2,2,4);  % Mechanical margin (Von Mises)
    plot(x_mm,mec_margin,'Color',c_violet,'LineWidth',2.5,'DisplayName','\sigma_y/(SF\cdot\sigma_{VM})'); hold on;
    yline(1,'Color',c_red,'LineStyle','--','LineWidth',1.5,'DisplayName','Failure');
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('margin [-]');
    title(sprintf('Mechanical margin (min=%.2f)',min(mec_margin))); legend('Location','best');

% -------------------------------------------------------------------
% Fig 7 — TBC thickness parametric sweep
% -------------------------------------------------------------------
t_TBC_um = t_TBC_vec * 1e6;
figure('Name','Fig 7 — TBC Sweep','Position',[160 160 1000 450]);
subplot(1,2,1);
    plot(t_TBC_um,T_wg_max_v-273.15,'Color',c_violet,'LineWidth',2.5,'DisplayName','T_{wg,max}'); hold on;
    yline(T_lim_wall-273.15,'Color',c_red,'LineStyle','--','LineWidth',2.0, ...
          'DisplayName',sprintf('T_{lim,CuCrZr}=%.0f C',T_lim_wall-273.15));
    xline(t_TBC*1e6,'Color',c_green,'LineStyle','--','LineWidth',1.5, ...
          'DisplayName',sprintf('Design %.0f \\mum',t_TBC*1e6));
    hold off; grid on; box on;
    xlabel('t_{TBC} [\mum]'); ylabel('T_{wg,max} [C]');
    title('T_{wg,max} vs t_{TBC}'); legend('Location','best');
subplot(1,2,2);
    plot(t_TBC_um,T_TBC_max_v-273.15,'Color',c_red,'LineWidth',2.5,'DisplayName','T_{TBC,max}'); hold on;
    yline(T_lim_TBC-273.15,'Color',c_orange,'LineStyle','--','LineWidth',2.0, ...
          'DisplayName',sprintf('T_{lim,YSZ}=%.0f C',T_lim_TBC-273.15));
    xline(t_TBC*1e6,'Color',c_green,'LineStyle','--','LineWidth',1.5, ...
          'DisplayName',sprintf('Design %.0f \\mum',t_TBC*1e6));
    hold off; grid on; box on;
    xlabel('t_{TBC} [\mum]'); ylabel('T_{TBC,max} [C]');
    title('T_{TBC,max} vs t_{TBC}'); legend('Location','best');

% -------------------------------------------------------------------
% Fig 8 — Coolant mass flow sensitivity
% -------------------------------------------------------------------
figure('Name','Fig 8 — mdot Sensitivity','Position',[180 180 800 450]);
plot(mdot_vec,T_wg_max_mdot-273.15,'Color',c_orange,'LineWidth',2.5,'DisplayName','T_{wg,max}'); hold on;
yline(T_lim_wall-273.15,'Color',c_red,'LineStyle','--','LineWidth',2.0, ...
      'DisplayName',sprintf('T_{lim,CuCrZr}=%.0f C',T_lim_wall-273.15));
xline(mdot_w,'Color',c_green,'LineStyle','--','LineWidth',1.5, ...
      'DisplayName',sprintf('Design %.1f kg/s',mdot_w));
if ~isnan(mdot_min_safe)
    xline(mdot_min_safe,'Color',c_yellow,'LineStyle',':','LineWidth',1.5, ...
          'DisplayName',sprintf('Min safe %.1f kg/s',mdot_min_safe));
end
hold off; grid on; box on;
xlabel('mdot_w [kg/s]'); ylabel('T_{wg,max} [C]');
title('T_{wg,max} vs coolant mass flow rate'); legend('Location','best');

% -------------------------------------------------------------------
% Fig 9 — TBC vs no-TBC comparison
% -------------------------------------------------------------------
figure('Name','Fig 9 — TBC vs no-TBC','Position',[200 200 1200 450]);
subplot(1,2,1);
    plot(x_mm,T_wg-273.15,       'Color',c_violet,'LineWidth',2.5, ...
         'DisplayName',sprintf('With TBC (%.0f\\mum)',t_TBC*1e6)); hold on;
    plot(x_mm,T_wg_noTBC-273.15, 'Color',c_red,   'LineWidth',2.5,'DisplayName','No TBC');
    yline(T_lim_wall-273.15,'Color',c_red,'LineStyle','--','LineWidth',1.5, ...
          'DisplayName',sprintf('T_{lim}=%.0f C',T_lim_wall-273.15));
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('T_{wg} [C]');
    title(sprintf('T_{wg}(x): TBC benefit = %.0f K',dT_TBC)); legend('Location','best');
subplot(1,2,2);
    plot(x_mm,q_flux/1e6,  'Color',c_violet,'LineWidth',2.5,'DisplayName','With TBC'); hold on;
    plot(x_mm,q_noTBC/1e6, 'Color',c_red,   'LineWidth',2.5,'DisplayName','No TBC');
    xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
    hold off; grid on; box on; xlabel('x [mm]'); ylabel('q'''' [MW/m^2]');
    title(sprintf('Wall heat flux: dq=%.2f MW/m^2 (%.0f%%)', ...
                  dq_TBC/1e6, 100*dq_TBC/max(q_noTBC)));
    legend('Location','best');

% -------------------------------------------------------------------
% Fig 10 — Coolant bulk temperature evolution
% -------------------------------------------------------------------
figure('Name','Fig 10 — Coolant Temperature','Position',[220 220 800 450]);
plot(x_fine_mm,T_w_f-273.15,'Color',c_cyan,'LineWidth',2.5,'DisplayName','T_w (bulk)'); hold on;
scatter(x_mm(1),  T_w(1)-273.15,  80,'o','filled','MarkerFaceColor',c_green, ...
        'DisplayName',sprintf('Inlet  %.1f C',T_w(1)-273.15));
scatter(x_mm(end),T_w(end)-273.15,80,'s','filled','MarkerFaceColor',c_red, ...
        'DisplayName',sprintf('Outlet %.1f C',T_w(end)-273.15));
xline(0,'Color',c_white,'LineStyle','--','LineWidth',1.2,'HandleVisibility','off');
hold off; grid on; box on;
xlabel('x [mm]'); ylabel('T_w [C]');
title(sprintf('Coolant bulk temperature  \\DeltaT=%.1f K',T_w(end)-T_in));
legend('Location','best');

% -------------------------------------------------------------------
% Fig 12 — N2 pressurant pressure during the test
%
% Each curve is computed using its own correctly-sized bottle volume, so
% both reach p_N2_final (lock-up pressure) exactly at t = t_burn_total.
% The adiabatic curve starts higher (larger bottle) and decays more slowly
% than the isothermal curve (smaller bottle, faster pressure drop).
% This is the physical consequence of the sizing: the adiabatic model
% requires more gas to compensate for the cooling-induced pressure loss,
% resulting in a safer (higher) pressure margin throughout the test.
% -------------------------------------------------------------------
figure('Name','Fig 12 — N2 Blow-down Pressure','Position',[240 80 900 450]);
plot(t_vec, p_N2_adi_t/1e5, 'Color',c_orange,'LineWidth',2.5, ...
     'DisplayName','Adiabatic (design, conservative)'); hold on;
plot(t_vec, p_N2_iso_t/1e5, 'Color',c_cyan,  'LineWidth',2.5, ...
     'DisplayName','Isothermal (reference, lower bound)');
yline(p_N2_final/1e5, 'Color',c_red,'LineStyle','--','LineWidth',2.0, ...
      'DisplayName',sprintf('Lock-up pressure = %.0f bar', p_N2_final/1e5));
xline(t_burn, 'Color',c_white,'LineStyle','--','LineWidth',1.2, ...
      'DisplayName',sprintf('t_{burn} = %.0f s', t_burn));
% Annotate final bottle pressure at end of test
text(t_burn_total*0.97, p_N2_adi_t(end)/1e5 + 5, ...
     sprintf('%.0f bar', p_N2_adi_t(end)/1e5), 'Color',c_orange, ...
     'FontSize',9,'HorizontalAlignment','right');
text(t_burn_total*0.97, p_N2_iso_t(end)/1e5 - 8, ...
     sprintf('%.0f bar', p_N2_iso_t(end)/1e5), 'Color',c_cyan, ...
     'FontSize',9,'HorizontalAlignment','right');
hold off; grid on; box on;
xlabel('t [s]'); ylabel('p_{N_2} [bar]');
title(sprintf(['N_2 blow-down pressure — ' ...
               'adiabatic: V=%.1f L, m=%.1f kg  |  ' ...
               'isothermal: V=%.1f L, m=%.1f kg'], ...
              V_N2_adi*1e3, m_N2_adi, V_N2_iso*1e3, m_N2_iso));
legend('Location','northeast');

% -------------------------------------------------------------------
% Fig 11 — Radial temperature profile at throat (schematic, not to scale)
%
% Shows the temperature distribution from the hot-gas free stream across
% all layers to the coolant bulk.  Widths are schematic (not to scale) to
% make thin coatings (TBC, bond) visible alongside the CuCrZr liner.
% All temperature values are taken directly from the marching solver at
% the throat node.
%
% Profile shapes:
%   Gas BL:     T(xi) = T_aw - (T_aw-T_TBC)*xi^3   [turbulent, Pr~0.63]
%   Solid layers: linear (Fourier, constant k per layer)
%   Coolant BL: T(xi) = T_wc + (T_w-T_wc)*xi^(1/3) [turbulent, high-Pr water]
% -------------------------------------------------------------------
figure('Name','Fig 11 — Radial T Profile at Throat','Position',[250 250 1080 600]);

% Interface temperatures at throat node [degC]
T_aw_r    = T_aw(i_th)   - 273.15;   % adiabatic wall T (free-stream recovery)
T_TBC_r   = T_TBC(i_th)  - 273.15;   % TBC hot face (gas-facing)
T_bond_r  = T_bond(i_th) - 273.15;   % bond coat hot face (= TBC cold face)
T_wg_r    = T_wg(i_th)   - 273.15;   % CuCrZr hot face (= bond coat cold face)
T_wc_r    = T_wc(i_th)   - 273.15;   % CuCrZr cold face (coolant-side)
T_w_r     = T_w(i_th)    - 273.15;   % coolant bulk temperature

% Schematic display widths [arbitrary units] — chosen for visual clarity
% Order: gas BL | TBC | bond coat | CuCrZr wall | coolant BL | coolant bulk
W_sch  = [2.2,  0.90,  0.45,  2.8,  1.5,  1.6];
x_e    = [0, cumsum(W_sch)];   % x-edges of each region

% y-axis range with margins
dT_rng = T_aw_r - T_w_r;
T_lo   = T_w_r   - 0.06 * dT_rng;
T_hi   = T_aw_r  + 0.10 * dT_rng;

% Layer background colours (semi-transparent)
c_bg_gas  = [0.95 0.55 0.25];   % warm orange — hot gas
c_bg_tbc  = [0.90 0.85 0.15];   % yellow      — YSZ ceramic
c_bg_bond = [0.65 0.65 0.90];   % light blue  — NiCoCrAlY bond
c_bg_wall = [0.25 0.80 0.85];   % cyan        — CuCrZr
c_bg_cool = [0.20 0.50 0.95];   % blue        — coolant
c_prof    = [0.95 0.20 0.10];   % red         — temperature profile line

hold on;

% --- Background shading for each region ----------------------------------
fill([x_e(1),x_e(2),x_e(2),x_e(1)],[T_lo,T_lo,T_hi,T_hi], ...
     c_bg_gas, 'FaceAlpha',0.18,'EdgeColor','none');
fill([x_e(2),x_e(3),x_e(3),x_e(2)],[T_lo,T_lo,T_hi,T_hi], ...
     c_bg_tbc, 'FaceAlpha',0.30,'EdgeColor','none');
fill([x_e(3),x_e(4),x_e(4),x_e(3)],[T_lo,T_lo,T_hi,T_hi], ...
     c_bg_bond,'FaceAlpha',0.30,'EdgeColor','none');
fill([x_e(4),x_e(5),x_e(5),x_e(4)],[T_lo,T_lo,T_hi,T_hi], ...
     c_bg_wall,'FaceAlpha',0.22,'EdgeColor','none');
fill([x_e(5),x_e(7),x_e(7),x_e(5)],[T_lo,T_lo,T_hi,T_hi], ...
     c_bg_cool,'FaceAlpha',0.18,'EdgeColor','none');

% --- Vertical interface lines --------------------------------------------
for ki = 2:6
    xline(x_e(ki),'Color',[0.65 0.65 0.65],'LineStyle','-','LineWidth',1.0, ...
          'HandleVisibility','off');
end

% --- Temperature profile (piecewise) ------------------------------------
N_bl = 300;

% 1. Gas boundary layer: turbulent, steep gradient at wall
%    xi=0 (free stream) -> T_aw;  xi=1 (wall) -> T_TBC
%    Shape: T = T_aw - (T_aw-T_TBC)*xi^3  [flat near free stream, steep at wall]
xi_g  = linspace(0, 1, N_bl);
x_gp  = x_e(1) + xi_g * W_sch(1);
T_gp  = T_aw_r - (T_aw_r - T_TBC_r) .* xi_g.^3;
plot(x_gp, T_gp, 'Color',c_prof,'LineWidth',2.8,'HandleVisibility','off');

% 2. TBC (YSZ): linear conduction, T_TBC -> T_bond
plot([x_e(2),x_e(3)],[T_TBC_r,T_bond_r],'Color',c_prof,'LineWidth',2.8,'HandleVisibility','off');

% 3. Bond coat (NiCoCrAlY): linear, T_bond -> T_wg
plot([x_e(3),x_e(4)],[T_bond_r,T_wg_r],'Color',c_prof,'LineWidth',2.8,'HandleVisibility','off');

% 4. CuCrZr liner: linear, T_wg -> T_wc
plot([x_e(4),x_e(5)],[T_wg_r,T_wc_r],'Color',c_prof,'LineWidth',2.8,'HandleVisibility','off');

% 5. Coolant boundary layer: high-Pr turbulent, steep gradient at wall
%    xi=0 (wall) -> T_wc;  xi=1 (bulk) -> T_w
%    Shape: T = T_wc + (T_w-T_wc)*xi^(1/3)  [steep at wall, flat toward bulk]
xi_c  = linspace(0, 1, N_bl);
x_cp  = x_e(5) + xi_c * W_sch(5);
T_cp  = T_wc_r + (T_w_r - T_wc_r) .* xi_c.^(1/3);
plot(x_cp, T_cp, 'Color',c_prof,'LineWidth',2.8,'HandleVisibility','off');

% 6. Coolant bulk: isothermal
plot([x_e(6),x_e(7)],[T_w_r,T_w_r],'Color',c_prof,'LineWidth',2.8,'HandleVisibility','off');

% --- Interface markers (hollow circles) ----------------------------------
T_pts = [T_aw_r, T_TBC_r, T_bond_r, T_wg_r, T_wc_r, T_w_r];
x_pts = [x_e(1),  x_e(2),  x_e(3),  x_e(4),  x_e(5), x_e(6)];
scatter(x_pts, T_pts, 75, 'o', 'filled', ...
        'MarkerFaceColor',c_white,'MarkerEdgeColor',c_prof, ...
        'LineWidth',1.8,'HandleVisibility','off');

% --- Light horizontal reference lines at key temperatures ----------------
lw_r = 0.7;
for kT = T_pts
    yline(kT,'Color',[0.55 0.55 0.55],'LineStyle',':','LineWidth',lw_r,'HandleVisibility','off');
end

% --- Temperature annotations (right side of plot) -----------------------
x_ann = x_e(7) + 0.08;
fsz   = 8.5;
text(x_ann, T_aw_r,   sprintf('T_{aw} = %.0f °C',  T_aw_r),   'FontSize',fsz,'Color',c_black,'FontWeight','bold');
text(x_ann, T_TBC_r,  sprintf('T_{TBC} = %.0f °C', T_TBC_r),  'FontSize',fsz,'Color',c_black);
text(x_ann, T_bond_r, sprintf('T_{bond} = %.0f °C',T_bond_r), 'FontSize',fsz,'Color',c_black);
text(x_ann, T_wg_r,   sprintf('T_{wg} = %.0f °C',  T_wg_r),   'FontSize',fsz,'Color',c_black);
text(x_ann, T_wc_r,   sprintf('T_{wc} = %.0f °C',  T_wc_r),   'FontSize',fsz,'Color',c_black);
text(x_ann, T_w_r,    sprintf('T_w = %.0f °C',     T_w_r),    'FontSize',fsz,'Color',c_black);

% --- Layer name labels at top of each region ----------------------------
T_lbl = T_hi - 0.025 * dT_rng;
fsz_l = 9;
text((x_e(1)+x_e(2))/2, T_lbl, {'Hot gas','BL'}, ...
     'HorizontalAlignment','center','FontSize',fsz_l,'Color',c_black,'FontWeight','bold');
text((x_e(2)+x_e(3))/2, T_lbl, {'YSZ','TBC'}, ...
     'HorizontalAlignment','center','FontSize',fsz_l,'Color',c_black,'FontWeight','bold');
text((x_e(3)+x_e(4))/2, T_lbl, {'Bond'}, ...
     'HorizontalAlignment','center','FontSize',fsz_l,'Color',c_black,'FontWeight','bold');
text((x_e(4)+x_e(5))/2, T_lbl, {'CuCrZr'}, ...
     'HorizontalAlignment','center','FontSize',fsz_l,'Color',c_black,'FontWeight','bold');
text((x_e(5)+x_e(6))/2, T_lbl, {'Coolant','BL'}, ...
     'HorizontalAlignment','center','FontSize',fsz_l,'Color',c_black,'FontWeight','bold');
text((x_e(6)+x_e(7))/2, T_lbl, {'Bulk'}, ...
     'HorizontalAlignment','center','FontSize',fsz_l,'Color',c_black,'FontWeight','bold');

% --- Physical thickness annotations at bottom of solid layers ------------
T_bot = T_lo + 0.018 * dT_rng;
text((x_e(2)+x_e(3))/2, T_bot, sprintf('%.0f µm',t_TBC*1e6), ...
     'HorizontalAlignment','center','FontSize',8,'Color',c_black);
text((x_e(3)+x_e(4))/2, T_bot, sprintf('%.0f µm',t_bond*1e6), ...
     'HorizontalAlignment','center','FontSize',8,'Color',c_black);
text((x_e(4)+x_e(5))/2, T_bot, sprintf('%.1f mm',t_wall*1e3), ...
     'HorizontalAlignment','center','FontSize',8,'Color',c_black);

% --- Arrows indicating heat flow direction (left to right) ---------------
annotation_x1 = (x_e(1) + 0.12) / x_e(7);
annotation_x2 = (x_e(1) + 0.50) / x_e(7);
annotation_y  = 0.08;
annotation('arrow',[annotation_x1, annotation_x2],[annotation_y, annotation_y], ...
           'Color',c_black,'LineWidth',1.5,'HeadWidth',8,'HeadLength',8);
text(x_e(1)+0.05, T_lo+0.07*dT_rng, 'q \rightarrow', ...
     'FontSize',9,'Color',c_black,'FontWeight','bold');

hold off;
box on; grid off;
set(gca,'XTick',[],'XTickLabel',{});
ylabel('Temperature [°C]','FontSize',11);
title(sprintf(['Radial temperature profile at throat  (x = %.1f mm)  —  ' ...
               'schematic widths, not to scale;  solver temperatures'], x(i_th)*1e3), ...
      'FontSize',10);
xlim([x_e(1), x_e(7) + 1.4]);
ylim([T_lo, T_hi]);


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function [T_w, T_wc, T_wg, T_bond, T_TBC, q_flux, h_g, sigma, ...
          rho_w_arr, cp_w_arr, mu_w_arr, k_w_arr, Pr_w_arr, ...
          v_w_arr, Re_arr, Nu_arr, h_c_arr, hc_eff_arr, eta_fin_arr, eta_o_arr, ...
          mu_wall_arr, visc_corr_arr, p_static_arr, q_rad_arr, n_iter_arr] = ...
    marching_solver(N, N_sig, tol_sig, T_aw, T_s_gas, T0, xi, AR, Bpf, ...
                    R_TBC_arr, R_bond, R_wall, ...
                    A_ch_arr, P_wet_arr, Dh_arr, r, AR_ch_arr, ...
                    mdot_w, mdot_ch, Nch, h_ch, t_fin, k_wall, ...
                    r_avg, ds, T_in, p_inlet, p_dump, K_inlet, K_outlet, ...
                    eps_g, sigma_SB)
% MARCHING_SOLVER  1-D axial thermal marching with coupled static pressure.
%
% At each node i the solver performs 5 steps:
%  (0) Pressure update: subtract Petukhov friction loss of element (i-1)->i.
%      Floor at p_dump (dump-tank back-pressure).
%  (1) IAPWS-IF97 Region 1: water properties at (T_w(i), p_node).
%  (2) Flow quantities: velocity v_i, Reynolds Re_i.
%  (3) Fixed-point iteration on T_hg (hot-gas-facing surface) and T_wc:
%       - Bartz sigma (S&B Eq.8-23) evaluated at T_hg.
%       - Gnielinski Nu with Petukhov f_D; Shah-London -5% for AR_ch>3.
%       - Fin efficiency (adiabatic tip, Incropera Eq.3.86/3.91).
%       - Radiation: q_rad = eps_g*sigma_SB*(T_s^4-T_hg^4); absorbed into
%         T_aw_eff = T_aw + q_rad/h_g (effective driving temperature).
%       - Series resistance network: 1/h_g | R_TBC | R_bond | R_wall | 1/hc_eff.
%       - Convergence tolerance: max(|dT_hg|, |dT_wc|) < tol_sig.
%  (4) Enthalpy march: T_w(i+1) = T_w(i) + q*2*pi*r_avg*ds/(mdot_w*cp).
%  K_outlet loss applied post-hoc to the last node.

[T_w,T_wc,T_wg,T_bond,T_TBC]          = deal(zeros(1,N));
[q_flux,h_g,sigma,q_rad_arr]           = deal(zeros(1,N));
[rho_w_arr,cp_w_arr,mu_w_arr,k_w_arr]  = deal(zeros(1,N));
[Pr_w_arr,v_w_arr,Re_arr,Nu_arr]       = deal(zeros(1,N));
[h_c_arr,hc_eff_arr]                   = deal(zeros(1,N));
[eta_fin_arr,eta_o_arr]                = deal(zeros(1,N));
[mu_wall_arr,visc_corr_arr]            = deal(zeros(1,N));
p_static_arr = zeros(1,N);
n_iter_arr   = zeros(1,N);

T_w(1) = T_in;

% Concentrated inlet loss (K_inlet) applied before first node
[rho_in_ref,~,~,~] = iapws_if97_r1(T_in, p_inlet);
v_in_ref   = mdot_ch / (rho_in_ref * A_ch_arr(1));
dp_K_inlet = 0.5 * rho_in_ref * v_in_ref^2 * K_inlet;
p_run      = max(p_dump, p_inlet - dp_K_inlet);

for i = 1:N

    % (0) Static pressure at node i
    if i > 1
        Re_prev = max(1, Re_arr(i-1));
        f_prev  = max(1e-3, (0.790*log(Re_prev)-1.64)^(-2));  % Petukhov
        dp_el   = f_prev * (rho_w_arr(i-1)*v_w_arr(i-1)^2/2) ...
                  / Dh_arr(i-1) * ds(i-1);
        p_run   = max(p_dump, p_run - dp_el);
    end
    p_static_arr(i) = p_run;

    % (1) Water properties: IAPWS-IF97 Region 1 at (T_w(i), p_run)
    [rho_i,cp_i,mu_i,k_i] = iapws_if97_r1(T_w(i), p_run);
    rho_w_arr(i)=rho_i; cp_w_arr(i)=cp_i;
    mu_w_arr(i) =mu_i;  k_w_arr(i) =k_i;
    Pr_i = mu_i*cp_i/k_i; Pr_w_arr(i)=Pr_i;

    % (2) Flow quantities
    v_i = mdot_ch / (rho_i*A_ch_arr(i));
    Re_i = rho_i*v_i*Dh_arr(i)/mu_i;
    v_w_arr(i)=v_i; Re_arr(i)=Re_i;

    % (3) Petukhov friction factor for Gnielinski
    f_D = max(1e-3, (0.790*log(max(1,Re_i))-1.64)^(-2));

    % (4) Fixed-point iteration on T_hg (hot-gas surface) and T_wc
    T_hg_guess = 0.85*T_aw(i);   % initial guess: hot-gas surface
    T_wc_guess = T_w(i) + 50;    % initial guess: cold wall

    % Safe defaults (cover the degenerate N_sig=0 case)
    sig_k=1; h_g_k=Bpf; q_k=0; q_rad_k=0;
    T_wc_k=T_w(i); T_wg_k=T_w(i); T_bond_k=T_w(i); T_TBC_k=T_w(i);
    Nu_k=0; h_c_k=0; hc_eff_k=1; mu_wall_k=mu_i; visc_k=1;
    eta_fin_k=1; eta_o_k=1; kk_end=N_sig;

    for kk = 1:N_sig

        % Bartz sigma correction (S&B Eq.8-23)
        sig_k = (T_hg_guess/(2*T0)*xi(i)+0.5)^(-0.68) * xi(i)^(-0.12);
        h_g_k = Bpf * AR(i)^(-0.9) * sig_k;

        % Gnielinski Nu (valid for Re > 3000, Pr > 0.5)
        T_wc_eval = min(620, max(274, T_wc_guess));
        [~,~,mu_wall_k,~] = iapws_if97_r1(T_wc_eval, p_run);
        visc_k = 1;  % (mu/mu_wall)^0.11 disabled — conservative choice
        num_gn = (f_D/8)*(Re_i-1000)*Pr_i;
        den_gn = 1 + 12.7*sqrt(f_D/8)*(Pr_i^(2/3)-1);
        Nu_k   = max(4, num_gn/den_gn * visc_k);  % floor at laminar limit (Re<1000 guard)
        if AR_ch_arr(i) > 3          % Shah & London correction for deep channels
            Nu_k = Nu_k * 0.95;
        end
        h_c_k = Nu_k * k_i / Dh_arr(i);

        % Fin efficiency (adiabatic tip, Incropera Eq.3.86/3.91)
        m_fin_k   = sqrt(2*h_c_k/(k_wall*t_fin));
        eta_fin_k = tanh(m_fin_k*h_ch)/(m_fin_k*h_ch);
        f_fin_k   = 2*h_ch/P_wet_arr(i);
        eta_o_k   = 1 - f_fin_k*(1-eta_fin_k);
        % Effective HTC normalised to hot-gas area (2*pi*r per unit length)
        hc_eff_k  = h_c_k*eta_o_k*(Nch*P_wet_arr(i))/(2*pi*r(i));

        % Gas-band radiation uses local static gas T (S&B Sec. 8.1)
        q_rad_k    = max(0, eps_g*sigma_SB*(T_s_gas(i)^4 - T_hg_guess^4));
        T_aw_eff_k = T_aw(i) + q_rad_k/h_g_k;  % effective driving temperature

        % Series thermal resistance network
        Rtot_k = 1/h_g_k + R_TBC_arr(i) + R_bond + R_wall + 1/hc_eff_k;
        q_k    = (T_aw_eff_k - T_w(i)) / Rtot_k;

        % Interface temperatures (cold to hot)
        T_wc_k   = T_w(i)   + q_k/hc_eff_k;
        T_wg_k   = T_wc_k   + q_k*R_wall;
        T_bond_k = T_wg_k   + q_k*R_bond;
        T_TBC_k  = T_bond_k + q_k*R_TBC_arr(i);

        % Hot-gas-facing surface = outermost layer
        if R_TBC_arr(i) > 0
            T_hg_k = T_TBC_k;   % with TBC
        else
            T_hg_k = T_wg_k;    % no TBC: copper hot face
        end

        % Convergence check
        if max(abs(T_hg_k-T_hg_guess), abs(T_wc_k-T_wc_guess)) < tol_sig
            kk_end = kk; break;
        end
        T_hg_guess = T_hg_k;
        T_wc_guess = T_wc_k;
        kk_end = kk;
    end % fixed-point loop

    % Store converged values
    h_g(i)=h_g_k; sigma(i)=sig_k; q_flux(i)=q_k; q_rad_arr(i)=q_rad_k;
    T_wc(i)=T_wc_k; T_wg(i)=T_wg_k; T_bond(i)=T_bond_k; T_TBC(i)=T_TBC_k;
    Nu_arr(i)=Nu_k; h_c_arr(i)=h_c_k; hc_eff_arr(i)=hc_eff_k;
    eta_fin_arr(i)=eta_fin_k; eta_o_arr(i)=eta_o_k;
    mu_wall_arr(i)=mu_wall_k; visc_corr_arr(i)=visc_k;
    n_iter_arr(i)=kk_end;

    % (5) Enthalpy march for T_w(i+1)
    if i < N
        T_w(i+1) = T_w(i) + q_flux(i)*2*pi*r_avg(i)*ds(i)/(mdot_w*cp_i);
    end

end % main axial loop

% Apply K_outlet loss to last node (post-hoc, floor at p_dump)
p_static_arr(N) = max(p_dump, p_static_arr(N) ...
                  - 0.5*rho_w_arr(N)*v_w_arr(N)^2*K_outlet);

end % marching_solver


function [rho, cp, mu, k] = iapws_if97_r1(T, p)
% IAPWS_IF97_R1  Liquid water: rho and cp from IAPWS-IF97 Region 1
%                (34-term dimensionless Gibbs free energy equation);
%                mu and k from NIST tables at 50 bar, linear interp on T
%                (error < 1% over the operating range).
%  Inputs:  T [K] (273.15 to 623.15),  p [Pa] (<= 100 MPa)
%  Outputs: rho [kg/m^3], cp [J/(kgK)], mu [Pa*s], k [W/(mK)]

    T = min(623.15, max(273.15, T));
    p = min(100e6,  max(611.0,  p));

    R_spec=461.526; p_ref=16.53e6; T_ref=1386.0;
    pi_=p/p_ref; tau=T_ref/T;

    I = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,5,8,8,21,23,29,30,31,32];
    J = [-2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,3,17,-4,0,6,-5,-2,10,-8,-11,-6,...
         -29,-31,-38,-39,-40,-41];
    n = [0.14632971213167,-0.84548187169114,-3.7563603672040, ...
         3.3855169168385,-0.95791963387872,0.15772038513228, ...
        -0.016616417199501,0.00081214629983568,0.00028319080123804, ...
        -0.00060706301565874,-0.018990068218419,-0.032529748770505, ...
        -0.021841717175414,-5.283835796993e-5,-4.7184321073267e-4, ...
        -3.0001780793026e-4,4.7661393906987e-5,-4.4141845330846e-6, ...
        -7.2694996297594e-16,-3.1679644845054e-5,-2.8270797985312e-6, ...
        -8.5205128120103e-10,-2.2425281908e-6,-6.5171222895601e-7, ...
        -1.4341729937924e-13,-4.0516996860117e-7,-1.2734301741641e-9, ...
        -1.7424871230634e-10,-6.8762131295531e-19,1.4478307828521e-20, ...
         2.6335781662795e-23,-1.1947622640071e-23,1.8228094581404e-24, ...
        -9.3537087292458e-26];

    dpi=7.1-pi_; dtau=tau-1.222;
    gamma_pi=0; gamma_tautau=0;
    for m=1:34
        Im=I(m); Jm=J(m); nm=n(m);
        if Im==0, base_pi=0; else, base_pi=-Im*dpi^(Im-1); end
        if Jm==0||Jm==1, base_tautau=0;
        else, base_tautau=Jm*(Jm-1)*dtau^(Jm-2); end
        gamma_pi     = gamma_pi     + nm*base_pi*dtau^Jm;
        gamma_tautau = gamma_tautau + nm*dpi^Im*base_tautau;
    end
    rho = p / ((pi_*gamma_pi)*R_spec*T);
    cp  = -tau^2 * gamma_tautau * R_spec;

    % mu and k: NIST tables at 50 bar, linear interpolation on T [C]
    T_tab  = [1,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,220,240,260,280,300];
    mu_tab = [1.731,1.307,1.002,0.798,0.653,0.547,0.467,0.404,0.354,0.315,0.282, ...
              0.232,0.196,0.170,0.150,0.134,0.122,0.112,0.104,0.098,0.092];
    k_tab  = [0.574,0.587,0.600,0.616,0.631,0.644,0.654,0.663,0.670,0.675,0.680, ...
              0.685,0.685,0.680,0.670,0.660,0.645,0.625,0.600,0.568,0.530];
    T_C = min(300, max(1, T-273.15));
    mu  = interp1(T_tab, mu_tab, T_C, 'linear') * 1e-3;
    k   = interp1(T_tab, k_tab,  T_C, 'linear');
end


function psat = p_sat_wagner(T)
% P_SAT_WAGNER  Saturation pressure via Wagner & Pruss (1993) Eq.2.5.
    Tc=647.096; pc=22.064e6;
    a1=-7.85951783; a2=1.84408259; a3=-11.7866497;
    a4=22.6807411;  a5=-15.9618719; a6=1.80122502;
    T=min(Tc-1e-6, max(273.16,T));
    theta=1-T/Tc;
    psat=pc*exp((Tc/T)*(a1*theta+a2*theta^1.5+a3*theta^3 ...
                        +a4*theta^3.5+a5*theta^4+a6*theta^7.5));
end


function Tsat = T_sat_wagner(p)
% T_SAT_WAGNER  Numerical inversion of p_sat_wagner via Newton-Raphson.
    Tc=647.096; pc=22.064e6;
    p=min(pc-1, max(611.657,p));
    T=min(Tc-1, max(274, 372.76+27*log10(p/101325)));  % initial guess
    for k=1:40
        ps=p_sat_wagner(T);
        if abs(ps-p)/p<1e-9, break; end
        dps=(p_sat_wagner(T+0.01)-p_sat_wagner(T-0.01))/0.02;
        T=min(Tc-0.5, max(274, T+sign((p-ps)/dps)*min(abs((p-ps)/dps),20)));
    end
    Tsat=T;
end


function s = margin_str(T_val, T_lim)
% MARGIN_STR  Returns a formatted margin string for console output.
    if T_val < T_lim
        s = sprintf('OK (margin %.0f K)', T_lim-T_val);
    else
        s = sprintf('EXCEEDS by %.0f K!', T_val-T_lim);
    end
end


function out = ternary(cond, a, b)
% TERNARY  Inline ternary operator.
    if cond, out=a; else, out=b; end
end
