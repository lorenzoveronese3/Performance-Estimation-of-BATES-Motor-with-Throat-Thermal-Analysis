%% ========================================================================
%  STARK-100k  -  UNIFIED DESIGN SCRIPT (quasi-steady chamber model)
%  ========================================================================
%  Combines:
%    - Outer Pc_in bisection to enforce max(Pc) = Pc_target_max (70 bar)
%    - Detailed nozzle sizing with CEA equilibrium data
%    - BATES grain design with nested bisection on web and ri/R
%    - Quasi-steady chamber model (mdot_gen = mdot_nozz instantaneously)
%    - Thermal protections sizing
%
%  Chamber model: Pc(x) = (rho_p * a * Ab(x) * c* / At)^(1/(1-n))
%  Time advancement: x(i+1) = x(i) + rb(i)*dt, with rb = a*Pc^n
%  ========================================================================
clear; close all; clc;

%% ========================================================================
%  CONSTANTS AND REQUIREMENTS
%  ========================================================================
g0      = 9.80665;       % gravitational acceleration    [m/s^2]
Ru      = 8314.3;        % universal gas constant        [J/(kmol.K)]
Pa_sl   = 1.01325 * 1e5; % atmospheric pressure (SL)     [Pa]

F_req   = 100e3;         % required thrust               [N]
It_req  = 2.5e6;         % required total impulse        [N.s]

%% ========================================================================
%  PROPELLANT DATA (AP/HTPB 84.18/15.82, for grain density)
%  ========================================================================
rho_AP   = 1950; % AP density [kg/m^3]
rho_HTPB = 915; % HTPB density [kg/m^3]
Wt_AP    = 0.8418; % Weight percentage of AP
Wt_HTPB  = 0.1582; % Weight percentage of HTPB
rho_p    = 1/((Wt_AP/rho_AP)+(Wt_HTPB/rho_HTPB)); % Density of the propellant [kg/m^3]

%% ========================================================================
%  THERMOCHEMISTRY FROM CEA (SHIFTING EQUILIBRIUM, O/F = 5.32)
%  ========================================================================
gamma_chamber = 1.1953;  % Specific heat ratio in the combustion chamber from CEA [-]
Mmol_chamber  = 23.683;  % Combustion gas molar mass in the chamber from CEA [kg/kmol]
Tc            = 2852.3;  % Chamber temperature from CEA [K]
rho_chamber   = 6.6911;  % Chamber gas density from CEA [kg/m^3]
R_gas         = Ru / Mmol_chamber; % Specific gas constant of the combustion products [J/(kg K)]

c_star_ideal  = 1542.0; % Ideal characteristic velocity from CEA [m/s]
CF_ideal      = 1.5754; % Ideal thrust coefficient from CEA [-]
A_ratio       = 8.1223; % Optimal nozzle expansion ratio Ae/At from CEA [-]

mu_gas_exit   = 0.53147e-4; % Exhaust gas dynamic viscosity at nozzle exit from CEA [Pa s]

eps_cea_sup = [1.0,    1.05,   1.15,   1.5,    2.0,    4.0,    6.0,    8.0,    8.1223]; % Supersonic nozzle area-ratio points A/A* from CEA [-]
rho_cea_sup = [4.1463, 3.2505, 2.6401, 1.7034, 1.1468, 0.49129,0.30884,0.22394,0.22021]; % Gas density along the supersonic branch [kg/m^3]
v_cea_sup   = [1047.9, 1272.8, 1431.6, 1700.9, 1894.2, 2211.1, 2344.8, 2425.4, 2429.7]; % Flow velocity along the supersonic branch [m/s]
a_cea_sup   = [1047.9, 1024.0, 1003.2, 958.8,  918.2,  832.8,  787.9,  757.7,  756.2]; % Local speed of sound along the supersonic branch [m/s]

eps_cea_sub = [10.0,   6.0,    4.0,    2.0,    1.5,    1.1,    1.0]; % Subsonic nozzle area-ratio points A/A* from CEA [-]
rho_cea_sub = [6.6793, 6.6581, 6.6161, 6.3727, 6.0812, 5.2605, 4.1463]; % Gas density along the subsonic convergent branch [kg/m^3]
v_cea_sub   = [64.54,  108.27, 163.95, 340.98, 476.27, 750.77, 1047.9]; % Flow velocity along the subsonic convergent branch [m/s]
a_cea_sub   = [1093.9, 1093.6, 1093.0, 1089.4, 1084.9, 1071.0, 1047.9]; % Local speed of sound along the subsonic branch [m/s]

alpha_deg = 15; % Divergent nozzle half-angle [deg]
beta_deg  = 45; % Convergent nozzle half-angle [deg]
lambda    = (1 + cosd(alpha_deg)) / 2; % Conical nozzle divergence correction factor [-]

N_disc = 150; % Number of points used to discretize each nozzle branch [-]

eps_sub_local = linspace(10.0, 1.0, N_disc); % Local area-ratio grid from convergent inlet to throat [-]
rho_sub_local = interp1(eps_cea_sub, rho_cea_sub, eps_sub_local, 'pchip'); % Interpolated density distribution on the subsonic branch [kg/m^3]
v_sub_local   = interp1(eps_cea_sub, v_cea_sub,   eps_sub_local, 'pchip'); % Interpolated velocity distribution on the subsonic branch [m/s]
q_sub_local   = 0.5 .* rho_sub_local .* (v_sub_local.^2); % Dynamic pressure distribution on the subsonic branch [Pa]

eps_sup_local = linspace(1.0, A_ratio, N_disc); % Local area-ratio grid from throat to nozzle exit [-]
rho_sup_local = interp1(eps_cea_sup, rho_cea_sup, eps_sup_local, 'pchip');  % Interpolated density distribution on the supersonic branch [kg/m^3]
v_sup_local   = interp1(eps_cea_sup, v_cea_sup,   eps_sup_local, 'pchip');  % Interpolated velocity distribution on the supersonic branch [m/s]
q_sup_local   = 0.5 .* rho_sup_local .* (v_sup_local.^2); % Dynamic pressure distribution on the supersonic branch [Pa]

%% ========================================================================
%  VIEILLE'S LAW FIT
%  ========================================================================
P_data  = [10.1  9.7 10.3  30.2 31.0 29.8 ...
           50.2 51.0 50.3  70.2 69.0 69.8 ...
           91.2 89.1 89.0];                          % Effective test pressure dataset [bar]
rb_data = [ 4.0  3.8  4.1   5.6  6.0  5.7 ...
            7.0  7.2  7.1   8.4  8.3  8.6 ...
            8.8  9.0  9.2];                          % Burning rate dataset [mm/s]

coeffs  = polyfit(log(P_data), log(rb_data), 1); % Linear regression in log-log space to identify Vieille's law coefficients
n_v     = coeffs(1);   % Pressure exponent of Vieille's law [-]
a_v     = exp(coeffs(2));   % Pre-exponential coefficient of Vieille's law [mm/s / bar^n]
rb_fit_pts = a_v * P_data.^n_v;    % Burning rates predicted by the fitted Vieille's law [mm/s]
R2 = 1 - sum((rb_data - rb_fit_pts).^2) / sum((rb_data - mean(rb_data)).^2);  % Coefficient of determination of the log-log fit [-]

% SI conversion of the Vieille coefficient:
%   rb[m/s]  = a_SI * Pc[Pa]^n
%   a_v is in [mm/s / bar^n], so:  a_SI = a_v * 1e-3 * (1e5)^(-n)
a_SI = a_v * 1e-3 * (1e5)^(-n_v);          % [m/s / Pa^n]

%% ========================================================================
%  BATES GEOMETRY FUNCTION
%  ========================================================================
grain_Ab = @(x, ri_g, R_g, L0_g) ...
    2*pi*(ri_g + x).*(L0_g - 2*x) + 2*pi*(R_g.^2 - (ri_g + x).^2); % Grain burning surface expressed as a handle function

tb_req = It_req / F_req; % Required burning time [s]
dt     = 0.01;     % time step for quasi-steady integration [s]

%% ========================================================================
%  OUTER BISECTION ON Pc_in
%  ========================================================================
Pc_target_max = 70;     % Maximum operating pressure target (MEOP) [bar]
Pc_in_low    = 52;      % Lower bound of the bisection bracket on Pc_in [bar]
Pc_in_high   = 77;      

for iter_Pc = 1:100      % Outer bisection loop, max 100 iterations

    Pc_in = (Pc_in_low + Pc_in_high) / 2;         % Current bisection midpoint [bar]       
    rb_nom = a_v * Pc_in^n_v;                     % Nominal rb at Pc_in [mm/s]
    p_c    = Pc_in * 1e5;                         % Pc_in in [Pa]

    %% --- NOZZLE SIZING --------------------------------------------------
    A_t = 0.0095;                                 % First guess no losses 
    F_error = 1000;                               % Break condition

    for iter_nozzle = 1:50                        % Inner loop: iterates At until thrust matches F_req = 100kN
        R_t = sqrt(A_t / pi);

        R_sub_start = sqrt((A_t * 10.0) / pi);                % Convergent inlet radius at eps=10, starting area considered for friction calculation (effect neglectible in convergent part, confirming hypothesis of eps 10) [m]          
        L_conv = (R_sub_start - R_t) / tand(beta_deg);        % Length convergent part from inlet to throat [m]
        A_local_sub = eps_sub_local .* A_t;                     % Local cross-sectional area on subsonic branch    [m^2]
        dA_w_sub    = abs(diff(A_local_sub)) ./ sind(beta_deg); % Elemental wall area on convergent for friction calculation [m^2]

        R_e   = sqrt((A_t * A_ratio) / pi);           % Exit radius at optimal expansion ratio [m]
        L_div = (R_e - R_t) / tand(alpha_deg);         % Divergent length from throat to exit (15° chosen as industry common use) [m]
        A_local_sup = eps_sup_local .* A_t;             % Local cross-sectional area on supersonic branch [m^2]
        dA_w_sup    = abs(diff(A_local_sup)) ./ sind(alpha_deg);   % Elemental wall area on divergent for friction calculation [m^2]

        L_total = L_conv + L_div;       % Total nozzle length used as characteristic length for Re [m]
        Re_L    = (rho_cea_sup(end) * v_cea_sup(end) * L_total) / mu_gas_exit;    % Reynolds number at nozzle exit (to overestimate Re and losses) [-]
        Cf      = 0.074 / (Re_L^0.2);           % Skin friction coefficient from turbulent flat-plate correlation [-]

        F_fric_conv  = sum(Cf .* q_sub_local(1:end-1) .* dA_w_sub);    % Friction force on convergent wall [N]
        F_fric_div   = sum(Cf .* q_sup_local(1:end-1) .* dA_w_sup);     % Friction force on divergent wall [N]
        F_fric_total = F_fric_conv + F_fric_div;          % Total viscous friction loss [N]

        delta_CF_visc = F_fric_total / (p_c * A_t);         % Viscous correction to thrust coefficient [-]
        CF_real       = (CF_ideal * lambda) - delta_CF_visc;      % Real thrust coefficient after losses (friction and div) [-]
        T_actual      = p_c * A_t * CF_real;                    % Actual thrust produced at current At [N]

        F_error = F_req - T_actual;                  % Thrust error to be driven to zero [N]
        if abs(F_error) < 0.1, break; end           % Convergence check: error below 0.1 N
        A_t = A_t * (F_req / T_actual);               % Rescale At proportionally to correct thrust
    end

    CF = CF_real; % Real thrust coefficient after applying divergence losses and viscous/friction losses [-]
    c_eff = CF * c_star_ideal; % Effective exhaust velocity, obtained as CF*c*, used to evaluate the specific impulse [m/s]
    Isp = c_eff / g0; % Effective specific impulse of the motor, computed from the effective exhaust velocity [s]
    At = A_t; % Final nozzle throat area obtained after the iterative nozzle sizing loop [m^2]
    Ae = A_local_sup(end); % Nozzle exit area, corresponding to the last point of the supersonic nozzle branch [m^2]
    Rt = sqrt(At / pi); % Nozzle throat radius, obtained from the throat area assuming a circular section [m]
    Re_exit = sqrt(Ae / pi); % Nozzle exit radius, obtained from the exit area assuming a circular section [m]
    L_nozzle = L_total; % Total nozzle length used in the preliminary nozzle-loss calculation [m]
    mdot_in = p_c * At / c_star_ideal; % Required motor mass flow rate at the selected chamber pressure, from c* = Pc*At/mdot [kg/s]

    %% --- BATES GEOMETRY -------------------------------------------------
    %  For quasi-steady model, imposing Pc(0) = Pc_in gives Ab(0) directly:
    %    Ab_req = At * Pc_in^(1-n) * 1e5^(1-n) / (rho_p * a * c*)
    %  which is equivalent to: mdot_in / (rho_p * rb_nom)
    Ab_req    = mdot_in / (rho_p * rb_nom * 1e-3);         % Required initial burning surface for mdot equilibrium [m^2]
    web_naive = rb_nom * 1e-3 * tb_req;                    % First estimate of web thickness: rb * tb [m]
    web_low   = web_naive * 0.70;                      % Lower bound of web bisection bracket [m]
    web_high  = web_naive * 1.30;                   % Upper bound of web bisection bracket [m]
    tol_web   = 1e-6;                        % Tolerance on burning-time error [s]

    for iter_web = 1:300             % Outer bisection on web thickness, max 300 iterations
        web_try = (web_low + web_high) / 2;         % Current web candidate (first is naive) [m]
          
        % Inner bisection on ri/R, finds the port-to-OD ratio that gives Ab(0) = Ab_req
        riR_lo = 0.05;     % Lower bound for ri/R [-]
        riR_hi = 0.95;     % Upper bound for ri/R [-]
        for iter_riR = 1:200              
            riR_mid = (riR_lo + riR_hi) / 2;              % Current ri/R candidate [-]
            R_k     = web_try / (1 - riR_mid);           % Outer radius from web and ri/R [m]
            ri_k    = riR_mid * R_k;                        % Inner (port) radius [m]
            L0_k    = 4*ri_k + 3*web_try;              % BATES grain length from neutrality condition [m]
            Ab0_k   = grain_Ab(0, ri_k, R_k, L0_k);           % Initial burning surface for this geometry [m^2]

            if Ab0_k > Ab_req, riR_hi = riR_mid;            % Surface too large: increase ri/R to reduce Ab
            else,              riR_lo = riR_mid; end           % Surface too small: increase ri/R
            if abs(Ab0_k - Ab_req)/Ab_req < 1e-12, break; end            % Convergence: relative error below 1e-12
        end
        riR_try = riR_mid;                           % Converged ri/R ratio [-]
        R_try   = web_try / (1 - riR_try);              % Corresponding outer radius [m]
        ri_try  = riR_try * R_try;                     % Corresponding inner radius [m]
        L0_try  = 4*ri_try + 3*web_try;                  % Corresponding grain length [m]

        % Quasi-steady integration to check if this geometry burns for tb_req seconds (all SI units)
        x_tmp  = 0;              % Burnt web thickness [m]
        t_tmp  = 0;             % Time burning current [s]
        while x_tmp < web_try && t_tmp < tb_req*3                  % Integrate until burnout or time limit
            Ab_tmp = grain_Ab(x_tmp, ri_try, R_try, L0_try);           % Current burning surface [m^2]
            Pc_tmp = (rho_p * a_SI * Ab_tmp * c_star_ideal / At)^(1/(1-n_v));   % Quasi-steady Pc [Pa]
            rb_tmp = a_SI * Pc_tmp^n_v;                                          % Current burning rate [m/s]
            x_tmp  = x_tmp + rb_tmp * dt;       % Advance burnt web thickness [m]
            t_tmp  = t_tmp + dt;                    % Advance time [s]
        end
        tb_calc = t_tmp;          % Burning time produced by this web candidate [s]

        if tb_calc < tb_req, web_low  = web_try;          % Grain burns too fast: increase web
        else,                web_high = web_try; end         % Grain burns too slow: decrease web
        if abs(tb_calc - tb_req) < tol_web, break; end         % Convergence: tb error below tolerance
    end

    web       = web_try;           % Converged web thickness from bisection [m]
    ri_over_R = riR_try;                 % Converged port-to-OD ratio from bisection [-]
    R_grain   = web / (1 - ri_over_R);       % Outer radius of the grain [m]
    ri        = ri_over_R * R_grain;            % Inner (port) radius of the grain [m]
    L0        = 4*ri + 3*web;            % Initial grain length from BATES neutrality condition [m]
    
    % Physical nozzle lengths (convergent starts at chamber diameter)
    eps_contraction = pi * R_grain^2 / At;        % real contraction ratio
    R_chamber = R_grain;                           % convergent inlet = grain OD
    L_conv_phys = (R_chamber - Rt) / tand(beta_deg);  % physical convergent [m]
    L_div_phys  = L_div;                              % divergent unchanged from nozzle sizing [m]
    L_nozzle_phys = L_conv_phys + L_div_phys;         % total physical nozzle length  [m]
    
    %% --- FINAL QUASI-STEADY TIME INTEGRATION ----------------------------
    N = floor(tb_req/dt) + 1;           % Maximum number of time steps for pre-allocation
    x_vec  = zeros(1, N);             % Burnt web thickness at each time step [m]
    Pc_vec = zeros(1, N);               % Chamber pressure at each time step [Pa]
    rb_vec = zeros(1, N);            % Burning rate at each time step [m/s]
    Ab_vec = zeros(1, N);      % Burning surface at each time step [m^2]
    t_vec  = zeros(1, N);               % Time at each time step [s]

    x_vec(1)  = 0;             % Initial burnt thickness: zero at ignition [m]
    Ab_vec(1) = grain_Ab(0, ri, R_grain, L0);        % Initial burning surface from grain geometry [m^2]
    Pc_vec(1) = (rho_p * a_SI * Ab_vec(1) * c_star_ideal / At)^(1/(1-n_v));   % Initial equilibrium Pc [Pa]
    rb_vec(1) = a_SI * Pc_vec(1)^n_v;                                          % Initial burning rate at Pc(0) [m/s]
    t_vec(1)  = 0;                                                  % Initial time [s]

    i = 2;
    while i <= N
        x_vec(i) = x_vec(i-1) + rb_vec(i-1) * dt;                 % Advance burnt thickness by one Euler step [m]
        if x_vec(i) >= web, break; end                         % Stop integration when full web is consumed

        Ab_vec(i) = grain_Ab(x_vec(i), ri, R_grain, L0);            % Updated burning surface [m^2]
        Pc_vec(i) = (rho_p * a_SI * Ab_vec(i) * c_star_ideal / At)^(1/(1-n_v));       % Updated equilibrium Pc [Pa]
        rb_vec(i) = a_SI * Pc_vec(i)^n_v;              % Updated burning rate [m/s]
        t_vec(i)  = t_vec(i-1) + dt;                  % Advance time [s]
        i = i + 1;
    end

    % Trim: remove unused pre-allocated entries beyond burnout
    t  = t_vec(1:i-1)';
    x  = x_vec(1:i-1)';
    Pc = Pc_vec(1:i-1)' / 1e5;                      % Chamber pressure converted to bar [bar]
    rb = rb_vec(1:i-1)';                            % [m/s]
    Ab = Ab_vec(1:i-1)';

    %% --- UPDATE OUTER BISECTION ------------------------------------------
    if max(Pc) > Pc_target_max
        Pc_in_high = Pc_in;            % Pressure overshot: lower the upper bound
    else
        Pc_in_low = Pc_in;            % Pressure undershot: raise the lower bound
    end
    if abs(max(Pc) - Pc_target_max) < 1e-4, break; end            % Convergence: max(Pc) within 1e-4 bar of target
end

fprintf('Pc_in = %.4f bar  =>  max(Pc) = %.4f bar  (%d outer iterations)\n\n', ...
         Pc_in, max(Pc), iter_Pc);

%% ========================================================================
%  FINAL PROFILES
%  ========================================================================
F         = CF * Pc * 1e5 .* At;          % Thrust profile over time [N]
mdot_gen  = rho_p * Ab .* rb;              % Mass flow rate generated by grain combustion [kg/s]
mdot_nozz = Pc * 1e5 .* At / c_star_ideal;      % Mass flow rate through the nozzle [kg/s]

Mp      = trapz(t, mdot_nozz);             % Total propellant mass expelled, integrated over time [kg]
It      = trapz(t, F);                   % Total impulse, integrated thrust over time [N.s]
Isp_eff = It / (Mp * g0);                % Effective specific impulse [s]
tb      = t(end);                         % Actual burning time from integration [s]
Ab0     = Ab(1);                             % Initial burning surface [m^2]
Ab_web  = Ab(end);                           % Final burning surface at web burnout [m^2]
Pc_eq_check = (rho_p * Ab0 * a_SI * c_star_ideal / At)^(1/(1-n_v)) / 1e5;  % Equilibrium Pc at t=0, cross-check [bar]

%% ========================================================================
%  THERMAL PROTECTIONS
%  ========================================================================
%  Ablative insulation for the combustion chamber walls.
%  The chamber is divided into three zones (forward dome, cylindrical wall,
%  aft dome) each subject to a different erosion regression rate.
%  Thickness is sized as:  t_insul = rdot * tb * SF  (with a minimum floor).
%  ========================================================================
mat_name    = 'EPDM + silica filler'; % Insulation material name
rho_insul   = 1200; % Insulation material density [kg/m^3]
k_insul     = 0.3;  % Thermal conductivity of the insulation [W/(m·K)]
T_decomp    = 700;  % Decomposition temperature of the insulation [K]
SF_insul    = 1.5;  % Safety factor applied to insulation thickness [-]
t_min_insul = 1.5;  % Minimum allowable insulation thickness [mm]

rdot_low    = 0.08; % Regression rate for the forward dome (low heat flux zone) [mm/s]
rdot_med    = 0.10; % Regression rate for the cylindrical wall (moderate heat flux zone) [mm/s]
rdot_high   = 0.15; % Regression rate for the aft dome / pre-nozzle region (high heat flux zone) [mm/s]

% Insulation thickness per zone: max between erosion-based sizing and structural minimum
t_insul_fwd = max(rdot_low  * tb * SF_insul, t_min_insul); % Forward dome insulation thickness [mm]
t_insul_cyl = max(rdot_med  * tb * SF_insul, t_min_insul); % Cylindrical wall insulation thickness [mm]
t_insul_aft = max(rdot_high * tb * SF_insul, t_min_insul); % Aft dome insulation thickness [mm]
t_insul_rad = t_min_insul; % Radial wall insulation thickness, set to structural minimum [mm]

% Wetted surface areas of the combustion chamber zones (simplified cylindrical + flat ends)
S_cyl = 2 * pi * R_grain * L0; % Lateral surface area of the cylindrical chamber wall [m^2]
S_fwd = pi * R_grain^2; % Surface area of the forward dome (flat disk approximation) [m^2]
S_aft = pi * R_grain^2; % Surface area of the aft dome (flat disk approximation) [m^2]

% Insulation mass per zone (thickness converted from mm to m via 1e-3)
M_insul_cyl = rho_insul * S_cyl * t_insul_cyl * 1e-3; % Cylindrical wall insulation mass [kg]
M_insul_fwd = rho_insul * S_fwd * t_insul_fwd * 1e-3; % Forward dome insulation mass [kg]
M_insul_aft = rho_insul * S_aft * t_insul_aft * 1e-3; % Aft dome insulation mass [kg]
M_insul_tot = M_insul_cyl + M_insul_fwd + M_insul_aft; % Total combustion chamber insulation mass [kg]

% Liner (HTPB-based, case bonding + plasticizer migration barrier)
rho_liner = 950;                                        % Liner density [kg/m^3]
t_liner   = 1.0;                                        % Liner thickness [mm]
M_liner   = rho_liner * (S_cyl + S_fwd + S_aft) * t_liner * 1e-3;  % Total liner mass [kg]

%% ========================================================================
%  PRINT RESULTS
%  ========================================================================
fprintf('=============================================================\n');
fprintf('  STARK-100k  -  UNIFIED DESIGN RESULTS\n');
fprintf('  Thermochemistry: NASA CEA equilibrium (O/F = 5.32)\n');
fprintf('  Nozzle: conical 15/45 deg, discretized CEA-based friction\n');
fprintf('  Chamber model: quasi-steady (instantaneous equilibrium)\n');
fprintf('=============================================================\n\n');

fprintf('===== PART A: DESIGN PARAMETERS ============================\n\n');

fprintf('--- Propellant data ---\n');
fprintf('  rho_prop = %.4f kg/m^3 (Grain density)\n\n', rho_p);

fprintf('--- Thermochemistry (CEA equilibrium, chamber) ---\n');
fprintf('  gamma_c     = %.4f (Specific heat ratio in the combustion chamber)\n', gamma_chamber);
fprintf('  Mmol_c      = %.3f kg/kmol (Combustion gas molar mass in the chamber)\n', Mmol_chamber);
fprintf('  Tc          = %.1f K (Chamber temperature)\n', Tc);
fprintf('  rho_c       = %.4f kg/m^3 (Chamber gas density)\n', rho_chamber);
fprintf('  R_gas       = %.2f J/(kg*K) (Specific gas constant of the combustion products)\n', R_gas);
fprintf('  c*_ideal    = %.1f m/s (Ideal characteristic velocity)\n', c_star_ideal);
fprintf('  CF_ideal    = %.4f (Ideal thrust coefficient)\n', CF_ideal);
fprintf('  epsilon     = %.4f   (Optimal nozzle expansion Pe = Pa)\n', A_ratio);
fprintf('  mu_gas_exit = %.3e Pa*s (Exhaust gas dynamic viscosity at nozzle exit)\n\n', mu_gas_exit);

fprintf('--- Vieille''s law ---\n');
fprintf('  n = %.5f [-] (Pressure exponent)\n', n_v);
fprintf('  a = %.4f [mm/(s*bar^n)] (Burning-rate coefficient)\n', a_v);
fprintf('  a_SI = %.4e [m/(s*Pa^n)] (Burning-rate coefficient expressed in SI units)\n', a_SI);
fprintf('  rb = %.4f * Pc^%.4f   [mm/s, bar] (Burning rate as a function of chamber pressure)\n', a_v, n_v);
fprintf('  rb(Pc_in)  = %.2f mm/s  (Burning rate at Pc_in = %.2f bar)\n', rb_nom, Pc_in);
fprintf('  R^2         = %.6f (Coefficient of determination of the log-log fit)\n\n', R2);

fprintf('--- Nozzle losses and coefficients ---\n');
fprintf('  alpha       = %d deg  (Divergent half-angle)\n', alpha_deg);
fprintf('  beta        = %d deg  (Convergent half-angle)\n', beta_deg);
fprintf('  lambda      = %.4f (Divergence loss factor)\n', lambda);
fprintf('  Cf_skin     = %.5f (Average viscous friction coefficient)\n', Cf);
fprintf('  F_fric_conv = %.2f N (Wall-friction loss contribution from the subsonic convergent branch)\n', F_fric_conv);
fprintf('  F_fric_div  = %.2f N (Wall-friction loss contribution from the supersonic divergent branch)\n', F_fric_div);
fprintf('  dCF_visc    = %.4f (Viscous drag contribution)\n', delta_CF_visc);
fprintf('  CF_real     = %.4f   (Real thrust coefficient)\n', CF);
fprintf('  c_eff       = %.1f m/s (Effective exhaust velocity)\n', c_eff);
fprintf('  Isp         = %.2f s (Effective specific impulse)\n\n', Isp);

fprintf('--- Nozzle geometry ---\n');
fprintf('  At       = %.2f cm^2 (Throat area) (Dt = %.1f mm)\n', At*1e4, 2*Rt*1e3);
fprintf('  Ae       = %.2f cm^2 (Exit area) (De = %.1f mm)\n', Ae*1e4, 2*Re_exit*1e3);
fprintf('  eps_c    = %.2f (Contraction ratio)\n', eps_contraction);
fprintf('  L_conv   = %.1f mm  (Physical length of the convergent, chamber excluded)\n', L_conv_phys*1e3);
fprintf('  L_div    = %.1f mm (Physical length of the divergent)\n', L_div_phys*1e3);
fprintf('  L_total  = %.1f mm  (Physical convergent+divergent length)\n', L_nozzle_phys*1e3);
fprintf('  Re_L     = %.2e (Reynolds number evaluated using exit-flow properties and the total nozzle length as characteristic length)\n', Re_L);
fprintf('  mdot_in = %.2f kg/s (Required motor mass flow rate)\n', mdot_in);
fprintf('  tb_req   = %.2f s (Required burning time)\n\n', tb_req);

fprintf('===== PART B: GRAIN GEOMETRY (BATES, N=1) ==================\n\n');
fprintf('  ri/R      = %.6f (Ratio between initial internal radius and external radius of the grain)\n', ri_over_R);
fprintf('  web       = %.2f mm (Thickness of the grain)\n', web*1e3);
fprintf('  R_grain   = %.1f mm (Initial external radius of the grain) (D_ext = %.0f mm)\n', R_grain*1e3, 2*R_grain*1e3);
fprintf('  ri        = %.1f mm (Initial internal radius of the grain) (D_bore = %.0f mm)\n', ri*1e3, 2*ri*1e3);
fprintf('  L0        = %.1f mm (Initial length of the grain) (BATES: 4ri+3w)\n', L0*1e3);
fprintf('  Ab(0)     = %.2f cm^2 (Initial burning surface)\n', Ab0*1e4);
fprintf('  Ab(web)   = %.2f cm^2 (Final burning surface) (Neutrality check: Ab(0) ≈ Ab(web))\n', Ab_web*1e4);
fprintf('  Pc_eq(0)  = %.4f bar (Initial equilibrium chamber pressure)\n\n', Pc_eq_check);

fprintf('===== PART C: TIME-DOMAIN (QUASI-STEADY) ===================\n\n');
fprintf('  Pc_in   = %.4f bar (Initial chamber pressure selected by the outer bisection)\n', Pc_in);
fprintf('  Pc(0)    = %.4f bar   Pc(tf) = %.4f bar (Initial and final chamber pressure)\n', Pc(1), Pc(end));
fprintf('  Pc_max   = %.4f bar (Maximum chamber pressure reached during firing) (target: %.1f bar)\n\n', max(Pc), Pc_target_max);

fprintf('  F(0)     = %.2f kN    F(tf) = %.2f kN (Initial and final thrust)\n', F(1)/1e3, F(end)/1e3);
fprintf('  F_max    = %.2f kN (Maximum thrust)\n\n', max(F)/1e3);

fprintf('  rb(0)    = %.3f mm/s  rb(tf) = %.3f mm/s (Initial and final burning rate)\n', rb(1)*1e3, rb(end)*1e3);
fprintf('  rb_max   = %.3f mm/s (Maximum burning rate)\n\n', max(rb)*1e3);

fprintf('  tb       = %.4f s (Burning time) (req: %.2f s)\n', tb, tb_req);
fprintf('  Mp       = %.2f kg (Total mass expelled from the nozzle)\n', Mp);
fprintf('  It       = %.4f MN.s (Total produced impulse) (req: %.3f MN.s)\n', It/1e6, It_req/1e6);
fprintf('  Isp_eff  = %.2f s (Effective specific impulse)\n\n', Isp_eff);

fprintf('===== PART D: THERMAL PROTECTIONS ===========================\n\n');

fprintf('--- Material: %s ---\n', mat_name);
fprintf('  rho_insul   = %d kg/m^3 (Insulation material density)\n', rho_insul);
fprintf('  k_insul     = %.2f W/(m*K) (Thermal conductivity of the insulation)\n', k_insul);
fprintf('  T_decomp    = %d K (Decomposition temperature of the insulation)\n', T_decomp);
fprintf('  SF_insul    = %.1f (Safety factor applied to insulation thickness)\n', SF_insul);
fprintf('  t_min_insul = %.1f mm (Minimum allowable insulation thickness)\n\n', t_min_insul);
fprintf('--- Liner ---\n');
fprintf('  rho_liner = %d kg/m^3 (Liner material density)\n', rho_liner);
fprintf('  t_liner   = %.1f mm (Liner thickness)\n', t_liner);
fprintf('  M_liner   = %.1f kg (Total liner mass)\n\n', M_liner);

fprintf('--- Insulation thicknesses (uniform per zone) ---\n');
fprintf('  Zone                     rdot [mm/s]   t_insul [mm]\n');
fprintf('  -------------------------------------------------------\n');
fprintf('  Forward dome             rdot_fwd   = %.2f      t_insul_fwd = %.1f\n', rdot_low,  t_insul_fwd);
fprintf('  Cylindrical wall         rdot_cyl   = %.2f      t_insul_cyl = %.1f\n', rdot_med,  t_insul_cyl);
fprintf('  Aft dome (pre-nozzle)    rdot_aft   = %.2f      t_insul_aft = %.1f\n', rdot_high, t_insul_aft);
fprintf('  Radial wall              -                      t_insul_rad = %.1f (structural minimum)\n\n', t_insul_rad);

fprintf('--- Wetted surface areas ---\n');
fprintf('  S_cyl = %.4f m^2 (Lateral surface area of the cylindrical chamber wall)\n', S_cyl);
fprintf('  S_fwd = %.4f m^2 (Surface area of the forward dome)\n', S_fwd);
fprintf('  S_aft = %.4f m^2 (Surface area of the aft dome)\n\n', S_aft);

fprintf('--- Insulation mass ---\n');
fprintf('  M_insul_cyl = %.1f kg (Cylindrical wall insulation mass)\n', M_insul_cyl);
fprintf('  M_insul_fwd = %.1f kg (Forward dome insulation mass)\n', M_insul_fwd);
fprintf('  M_insul_aft = %.1f kg (Aft dome insulation mass)\n', M_insul_aft);
fprintf('  M_insul_tot = %.1f kg (Total chamber insulation mass) (%.1f%% of Mp)\n\n', ...
         M_insul_tot, M_insul_tot/Mp*100);
fprintf('=============================================================\n');

%% ========================================================================
%  FIGURES
%  ========================================================================

% --- FIGURE 1: Design parameters & nozzle ---
f1 =figure('Position', [50 50 1400 500], 'Color', 'w');
set(f1, 'WindowStyle', 'docked');

subplot(1,3,1);
hold on; grid on; box on;
P_fit = linspace(5, 100, 200);
plot(P_data, rb_data, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', [0.8 0.2 0.2]);
plot(P_fit, a_v*P_fit.^n_v, 'b-', 'LineWidth', 2);
plot(Pc_in, rb_nom, 'gs', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 10);
xlabel('Pressure [bar]'); ylabel('r_b [mm/s]');
title(sprintf('Vieille: r_b = %.3f p^{%.3f}  (R^2=%.4f)', a_v, n_v, R2));
legend('Data', 'Fit', sprintf('%.1f bar', Pc_in), 'Location', 'southeast');

subplot(1,3,2);
hold on; grid on; box on;
eps_sub_plot = linspace(eps_contraction, 1.0, N_disc);
A_local_sub_plot = eps_sub_plot .* At;
x_coords_sub_plot = -(sqrt(A_local_sub_plot/pi) - Rt) ./ tand(beta_deg);
x_coords_sup_plot = (sqrt(A_local_sup/pi) - Rt) ./ tand(alpha_deg);
x_full    = [x_coords_sub_plot, x_coords_sup_plot];
R_profile = sqrt([A_local_sub_plot, A_local_sup]./pi) * 1e3;
plot(x_full*1e3, R_profile, 'k-', 'LineWidth', 2);
plot(x_full*1e3, -R_profile, 'k-', 'LineWidth', 2);
xlabel('Axial position [mm]'); ylabel('Radius [mm]');
title('Nozzle profile'); axis equal; set(gca, 'FontSize', 10);

subplot(1,3,3);
hold on; grid on; box on;
a_sub_loc = interp1(eps_cea_sub, a_cea_sub, eps_sub_local, 'pchip');
a_sup_loc = interp1(eps_cea_sup, a_cea_sup, eps_sup_local, 'pchip');
Mach_full = [v_sub_local ./ a_sub_loc, v_sup_local ./ a_sup_loc];
plot(x_full*1e3, Mach_full, 'b-', 'LineWidth', 2);
yline(1.0, 'r--', 'LineWidth', 1.5);
xlabel('Axial position [mm]'); ylabel('Mach');
title('Mach evolution'); set(gca, 'FontSize', 10);

sgtitle('PART A  |  Design Parameters & Nozzle', ...
        'FontSize', 14, 'FontWeight', 'bold');

% --- FIGURE 2: Time-domain profiles ---
f2 =figure('Position', [50 50 1400 900], 'Color', 'w');
set(f2, 'WindowStyle', 'docked');
subplot(2,3,1);
hold on; grid on; box on;
plot(t, F/1e3, 'b-', 'LineWidth', 2.5);
plot(t(1), F(1)/1e3, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(t(end), F(end)/1e3, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
yline(F(1)/1e3, '--', 'Color', [0.6 0.6 0.6]);
xlabel('Time [s]'); ylabel('F [kN]'); title('Thrust profile');
legend('F(t)', sprintf('F(0)=%.1f kN', F(1)/1e3), ...
       sprintf('F(t_f)=%.1f kN', F(end)/1e3), 'Location', 'south');
set(gca, 'FontSize', 10);

subplot(2,3,2);
hold on; grid on; box on;
plot(t, Pc, 'b-', 'LineWidth', 2.5);
plot(t(1), Pc(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(t(end), Pc(end), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
yline(Pc_target_max, 'r--', 'LineWidth', 1.5);
yline(Pc_in, '--', 'Color', [0.6 0.6 0.6]);
xlabel('Time [s]'); ylabel('P_c [bar]');
title(sprintf('Chamber pressure (Pc_{max}=%.2f bar)', max(Pc)));
legend('P_c(t)', sprintf('P_c(0)=%.2f', Pc(1)), ...
       sprintf('P_c(t_f)=%.2f', Pc(end)), ...
       sprintf('Pc_{target}=%.1f', Pc_target_max), 'Location', 'south');
set(gca, 'FontSize', 10);

subplot(2,3,3);
hold on; grid on; box on;
plot(t, rb*1e3, 'Color', [0.8 0.2 0.2], 'LineWidth', 2.5);
yline(rb_nom, '--', 'Color', [0.6 0.6 0.6]);
xlabel('Time [s]'); ylabel('r_b [mm/s]');
title('Local burning rate');
legend('r_b(t)', sprintf('r_b(%.1f bar)=%.2f', Pc_in, rb_nom), 'Location', 'south');
set(gca, 'FontSize', 10);

subplot(2,3,4);
hold on; grid on; box on;
plot(t, Ab*1e4, 'Color', [0.2 0.6 0.2], 'LineWidth', 2.5);
plot(t(1), Ab(1)*1e4, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(t(end), Ab(end)*1e4, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
yline(Ab(1)*1e4, '--', 'Color', [0.6 0.6 0.6]);
xlabel('Time [s]'); ylabel('A_b [cm^2]');
title('Burning area');
legend('A_b(t)', sprintf('A_b(0)=%.0f', Ab(1)*1e4), ...
       sprintf('A_b(t_f)=%.0f', Ab(end)*1e4), 'Location', 'south');
set(gca, 'FontSize', 10);

subplot(2,3,5);
hold on; grid on; box on;
plot(t, mdot_gen,  'Color', [0.8 0.5 0],  'LineWidth', 2);
plot(t, mdot_nozz, 'Color', [0.6 0.2 0.8], 'LineWidth', 2);
yline(mdot_in, '--', 'Color', [0.6 0.6 0.6]);
xlabel('Time [s]'); ylabel('Mass flow [kg/s]');
title('Mass flow rates');
legend('mdot_{gen}', 'mdot_{nozz}', 'Location', 'south');
set(gca, 'FontSize', 10);

subplot(2,3,6);
hold on; grid on; box on;
fill([t; flipud(t)], [F/1e3; F(1)/1e3*ones(size(F))], ...
     [0.3 0.5 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(t, F/1e3, 'b-', 'LineWidth', 2);
yline(F(1)/1e3, '--', 'Color', [0.6 0.6 0.6]);
xlabel('Time [s]'); ylabel('F [kN]');
title(sprintf('Total impulse = %.4f MN{\\cdot}s', It/1e6));
legend('Excess vs nominal', 'F(t)', 'F_{nom}', 'Location', 'south');
set(gca, 'FontSize', 10);

sgtitle(sprintf(['BATES (N=1, quasi-steady)  |  r_i/R=%.4f  |  ' ...
    'P_c(0)\\approx%.1f bar  |  t_b=%.2f s  |  I_t=%.4f MN{\\cdot}s'], ...
    ri_over_R, Pc(1), tb, It/1e6), ...
    'FontSize', 13, 'FontWeight', 'bold');


%% ========================================================================
%  NOZZLE ABLATIVE THERMAL PROTECTION
%  ========================================================================
%  Ablative insulation on the nozzle sections NOT covered by the ceramic
%  throat insert (designed by colleagues):
%    - Convergent ablative: from chamber entrance (eps ≈ 62.5) down to eps = 2
%    - Divergent ablative:  from eps = 2 out to the nozzle exit (eps = 8.12)
%
%  The ceramic insert covers the throat region from eps = 2 (convergent
%  side) through the throat to eps = 2 (divergent side).
%
%  Material: EPDM + silica filler (same as combustion chamber insulation)
%  Sizing method:  t_insul = rdot * tb * SF  (with minimum floor)
%
%  Regression rates are higher than in the chamber because gas velocities
%  in the nozzle are substantially greater, which increases convective
%  heat transfer to the walls.
%
%  Ref: Sutton & Biblarz, "Rocket Propulsion Elements", 9th ed.
%  ========================================================================

%% NOZZLE ABLATIVE GEOMETRY

% Slant (wall) lengths of the two ablative frustum sections, from the nozzle geometry
l_conv_abl_wall = 0.032737;    % Convergent ablative section wall (slant) length [m]
l_div_abl_wall  = 0.089440;    % Divergent ablative section wall (slant) length [m]

% Radii at the boundaries between ablative and ceramic sections
r_t     = sqrt(At / pi);                % Throat radius [m]
r_eps2  = r_t * sqrt(2);                % Radius at area ratio eps = 2 (ceramic/ablative interface) [m]
r_inlet = r_t * sqrt(62.4948);          % Radius at the convergent inlet (chamber-side entrance) [m]
r_exit  = sqrt(Ae / pi);                % Nozzle exit radius [m]

% Wetted wall areas of the ablative conical frustums (lateral surface of a truncated cone)
A_wall_conv_abl = pi * (r_inlet + r_eps2) * l_conv_abl_wall;   % Convergent ablative wetted wall area [m^2]
A_wall_div_abl  = pi * (r_eps2  + r_exit) * l_div_abl_wall;    % Divergent ablative wetted wall area [m^2]

%% NOZZLE ABLATIVE REGRESSION RATES
% Assigned higher than chamber values due to elevated gas velocity in the nozzle.
% Convergent (subsonic, accelerating flow): moderate-to-high convective load.
% Divergent (supersonic, but density drops away from throat): moderate convective load.
rdot_conv_noz = 0.12;          % Ablative regression rate on the convergent nozzle wall [mm/s]
rdot_div_noz  = 0.10;          % Ablative regression rate on the divergent nozzle wall [mm/s]

%% NOZZLE ABLATIVE THICKNESSES
% Sized with the same approach as the chamber: erosion depth × safety factor, floored by t_min
t_noz_conv = max(rdot_conv_noz * tb * SF_insul, t_min_insul);   % Convergent nozzle ablative thickness [mm]
t_noz_div  = max(rdot_div_noz  * tb * SF_insul, t_min_insul);   % Divergent nozzle ablative thickness [mm]

%% NOZZLE ABLATIVE MASS
% Mass of each ablative layer (thickness converted from mm to m via 1e-3)
M_noz_conv = rho_insul * A_wall_conv_abl * t_noz_conv * 1e-3;   % Convergent nozzle ablative mass [kg]
M_noz_div  = rho_insul * A_wall_div_abl  * t_noz_div  * 1e-3;   % Divergent nozzle ablative mass [kg]
M_noz_tot  = M_noz_conv + M_noz_div;                             % Total nozzle ablative mass [kg]

%% PRINT NOZZLE ABLATIVE RESULTS
fprintf('===== PART E: NOZZLE ABLATIVE PROTECTION ===================\n\n');

fprintf('--- Material: %s (same as chamber) ---\n', mat_name);
fprintf('  SF_insul    = %.1f (Safety factor)\n', SF_insul);
fprintf('  t_min_insul = %.1f mm (Minimum allowable thickness)\n\n', t_min_insul);

fprintf('--- Nozzle ablative geometry ---\n');
fprintf('  Convergent: eps = 62.5 to eps = 2\n');
fprintf('    l_conv_abl_wall = %.1f mm (Convergent ablative wall slant length)\n', l_conv_abl_wall*1e3);
fprintf('    r_inlet         = %.1f mm (Convergent inlet radius)\n', r_inlet*1e3);
fprintf('    r_eps2          = %.1f mm (Radius at area ratio 2)\n', r_eps2*1e3);
fprintf('    A_wall_conv_abl = %.4f m^2 (Convergent ablative wetted wall area)\n\n', A_wall_conv_abl);

fprintf('  Divergent: eps = 2 to eps = 8.12\n');
fprintf('    l_div_abl_wall  = %.1f mm (Divergent ablative wall slant length)\n', l_div_abl_wall*1e3);
fprintf('    r_eps2          = %.1f mm (Radius at area ratio 2)\n', r_eps2*1e3);
fprintf('    r_exit          = %.1f mm (Nozzle exit radius)\n', r_exit*1e3);
fprintf('    A_wall_div_abl  = %.4f m^2 (Divergent ablative wetted wall area)\n\n', A_wall_div_abl);

fprintf('--- Ablative thicknesses ---\n');
fprintf('  Zone                     rdot [mm/s]           t_insul [mm]\n');
fprintf('  -------------------------------------------------------\n');
fprintf('  Convergent (eps>2)       rdot_conv_noz = %.2f   t_noz_conv = %.1f\n', rdot_conv_noz, t_noz_conv);
fprintf('  Divergent  (eps>2)       rdot_div_noz  = %.2f   t_noz_div  = %.1f\n\n', rdot_div_noz, t_noz_div);

fprintf('--- Nozzle ablative mass ---\n');
fprintf('  M_noz_conv = %.2f kg (Convergent nozzle ablative mass)\n', M_noz_conv);
fprintf('  M_noz_div  = %.2f kg (Divergent nozzle ablative mass)\n', M_noz_div);
fprintf('  M_noz_tot  = %.2f kg (Total nozzle ablative mass)\n\n', M_noz_tot);

fprintf('--- Total thermal protection (chamber + liner + nozzle) ---\n');
fprintf('  M_insul_tot = %.2f kg (Chamber insulation mass)\n', M_insul_tot);
fprintf('  M_liner     = %.2f kg (Liner mass)\n', M_liner);
fprintf('  M_noz_tot   = %.2f kg (Nozzle ablative mass)\n', M_noz_tot);
fprintf('  M_th_tot    = %.2f kg (Grand total thermal protection mass) (%.1f%% of Mp)\n\n', ...
         M_insul_tot + M_liner + M_noz_tot, (M_insul_tot + M_liner + M_noz_tot)/Mp*100);









