clear; clc; close all;

% This code was written and tested in MATLAB R2025b

fprintf('Please wait, the code may take a while to run\n')

%-----------------------------
% SECTION 1: inputs parameters
%-----------------------------

% GAS
T0 = 2852.3;
Mmol = 23.683;
P0 = 63.77e5;
C_star = 1542.0;
cp_0_fr = 1899.1;
Pr_0_fr = 0.6268;
mu_0_fr = 9.1736e-5;


% ENGINE
F    = 100e3;    
Isp  = 240.4;   
Itot = 2.5e6;  
g0 = 9.80665;
mdot_prop = F / (Isp * g0);                     
t_burn    = Itot / F; 

% NOZZLE GEOMETRY 
AR_in      = 2.0;    
AR_out     = 2.0;    
a_conv_deg = 45.0;   % [deg] Convergent half-angle.
a_div_deg  = 15.0;   % [deg] Divergent half-angle.
throat_curvature_ratio       = 1.5;    % Rc/Rt.

throat_area  = mdot_prop * C_star / P0;                     
throat_radius  = sqrt(throat_area / pi);                          
throat_diameter  = 2 * throat_radius;                                
throat_curvature_radius  = throat_curvature_ratio * throat_radius; 

throat_idx = 42;


% LAYERS
k_YSZ      = 1.5;        
T_lim_YSZ  = 1400+273;   

t_bond     = 0.035e-3;   
k_bond     = 20.0;   
R_bond = t_bond / k_bond;         

t_wall     = 1e-3;    
k_wall     = 320.0;      
T_lim_wall = 773;  
R_wall = t_wall  / k_wall;  

% COOLING CHANNELS
fin_th = 1.5e-3; 

% WATER
T_in   = 18 + 273.15;   
mdot_w_hp = 5;        % Initial hypotesized mdot for water          
rho_w  = 983;           
cp_w   = 4183;           
k_w    = 0.654;

% NUMERICAL PARAMETERS 
N      = 200;   % Number of axial sections.
N_sig   = 15;    % Max sigma fixed-point iterations per section.
tol_sig = 0.5;   % Convergence tolerance on T_TBC.

% CEA GAS PROPERTY TABLES
%  AR_tab: area ratio sample points [1.0 ... 2.0], shared by both sides.
%  mu in [Pa·s] (converted from mPoise: 1 mPoise = 1e-4 Pa·s).
%  Cp in [J/(kgK)] (converted from kJ/(kgK) * 1e3).

% Convergent — equilibrium model:
AR_eq_tab  = [1.00,   1.05,   1.10,   1.20,   1.50,   2.00];   
gamma_eq_tab  = [1.2076, 1.2029, 1.2014, 1.1998, 1.1977, 1.1965];        % [-]
cp_gas_eq_tab = [1875.1, 1884.7, 1887.5, 1890.7, 1894.6, 1896.8];       % [J/(kgK)]
Pr_eq_tab     = [0.6306, 0.6292, 0.6288, 0.6283, 0.6276, 0.6272];       % [-]
mu_eq_tab     = [8.5925, 8.8139, 8.8827, 8.9598, 9.0584, 9.1149] * 1e-5; % [Pa·s]
T_eq_tab      = [2598.28, 2694.19, 2724.21, 2757.98, 2801.34, 2826.29];       % [K]
P_eq_tab      = [37.704, 46.920, 50.222, 54.203, 59.761, 63.207] * 1e5;  % [Pa]

% Divergent — frozen model:
AR_fr_tab     = [1.00,   1.05,   1.15,   1.50,   2.00];
gamma_fr_tab  = [1.2076, 1.2138, 1.2187, 1.2274, 1.2333];       % [-]
cp_gas_fr_tab = [1875.1, 1861.5, 1849.0, 1820.4, 1792.1];       % [J/(kgK)]
Pr_fr_tab     = [0.6306, 0.6323, 0.6335, 0.6356, 0.6365];       % [-]
mu_fr_tab     = [8.5925, 8.2944, 8.0398, 7.5091, 7.0429] * 1e-5; % [Pa·s]
T_fr_tab      = [2598.28, 2470.88, 2363.69, 2145.11, 1958.71]; % [K]
P_fr_tab      = [37.704, 28.080, 21.804, 12.758, 7.8401] * 1e5; % [Pa]

% Mach number from CEA (avoids fzero — CEA integrates Euler with variable gamma):
M_eq_tab      = [1.000,  0.781,  0.701,  0.600,  0.439,  0.313];          
M_fr_tab      = [1.000,  1.243,  1.427,  1.774,  2.063];

% OTHER QUANTITIES

universal_gas_const = 8314.5;
specific_gas_const = universal_gas_const/Mmol;
                                                       

a_conv_rad = a_conv_deg * pi/180;                   % [rad]
a_div_rad  = a_div_deg  * pi/180;                   % [rad]
R_in   = throat_radius * sqrt(AR_in);                      % [m]
R_out  = throat_radius * sqrt(AR_out);                     % [m]
L_conv = (R_in  - throat_radius) / tan(a_conv_rad);            % [m]
L_div  = (R_out - throat_radius) / tan(a_div_rad);             % [m]
L_tot  = L_conv + L_div;                        % [m]

% Bartz prefactor — chamber (stagnation) CEA frozen properties
Bpf = (0.026 / throat_diameter^0.2) * (throat_diameter/throat_curvature_radius)^0.1 ...
    * (mu_0_fr^0.2 * cp_0_fr / Pr_0_fr^0.6) ...
    * (P0 / C_star)^0.8; 

% Setting geometry arrays
[n_conv, n_div, x_conv, x_div, x, local_area_ratio, local_diameter, dx_vec, dr_vec, ds, ds_fixed, r_avg] = nozzle_geometry_section_array_plot( ...
    N, L_conv, L_div, L_tot, throat_radius, throat_area, a_conv_rad, a_div_rad);

% Setting CEA arrays
[gamma_local, cp_local, Pr_local, mu_local, T_static, P_static, M] = cea_interpolations_plot(N, ...
    AR_eq_tab, AR_fr_tab, local_area_ratio, gamma_eq_tab, cp_gas_eq_tab, Pr_eq_tab, mu_eq_tab, T_eq_tab, P_eq_tab, M_eq_tab, gamma_fr_tab, cp_gas_fr_tab, Pr_fr_tab, mu_fr_tab, T_fr_tab, P_fr_tab, M_fr_tab, x);

%--------------------------------------------
% SECTION 2: sizing of TBC and water channels
%--------------------------------------------

% YSZ thickness computation
[hc_no_TBC, t_YSZ, R_YSZ] = tbc_thickness_calculation_plot(throat_idx, ...
    cp_w, mdot_w_hp, T_in, T_lim_wall, T0, P0, ...
    C_star, local_area_ratio, local_diameter, ds, gamma_local, ...
    throat_diameter, throat_curvature_radius, mu_0_fr, cp_local, Pr_local, t_wall, k_wall, k_YSZ, ...
    T_lim_YSZ, t_bond, k_bond, M, N);

% Water channel sizing
[a_opt, Tw_opt, Twc_opt, v_h2o, q_gas, Marg_boil, N_ch, b_ch_ottimale, mdot_fin, P_inlet_ottimale] = ...
    water_channel_sizing_dp_min_plot(N, T0, gamma_local, Pr_local, M, throat_diameter, ...
    throat_curvature_radius, mu_0_fr, cp_local, P0, C_star, local_area_ratio, T_lim_wall, ...
    t_YSZ, k_YSZ, t_bond, k_bond, t_wall, k_wall, local_diameter, ds, T_in, cp_w, k_w, rho_w, throat_idx);


%------------------------
%%% PLOTS AND OUTPUTS %%%
%------------------------

% note: not all the plots and outputs are in this section, some
% are in their specific function

% Function to plot the research space
[k] = water_channel_sizing_dp_min_SS(N, T0, gamma_local, Pr_local, M, throat_diameter, throat_curvature_radius, mu_0_fr, cp_local, P0, C_star, local_area_ratio, T_lim_wall, ...
   t_YSZ, k_YSZ, t_bond, k_bond, t_wall, local_diameter, ds, T_in, cp_w, k_w, rho_w);

%------------------------
%%% FUNCTIONS %%%
%------------------------

function [gamma_local, cp_local, Pr_local, mu_local, T_static, P_static, M] = cea_interpolations_plot( N,...
    AR_eq_tab, AR_fr_tab, local_area_ratio, gamma_eq_tab, cp_gas_eq_tab, Pr_eq_tab, mu_eq_tab, T_eq_tab, P_eq_tab, M_eq_tab, gamma_fr_tab, cp_gas_fr_tab, Pr_fr_tab, mu_fr_tab, T_fr_tab, P_fr_tab, M_fr_tab, x)

% This functions interpolates the values from CEA into arrays of N elements

% INPUTS:
% Tables with values from CEA to be interplated.
% 'eq' index for equilibrium (convergent) and 'fr' index for frozen (divergent)

% OUTPUTS:
% gamma_local
% cp_local
% Pr_local
% mu_local
% T_static
% P_static
% M for Mach 

gamma_local  = zeros(1, N);   
cp_local   = zeros(1, N);   
Pr_local   = zeros(1, N);   
mu_local   = zeros(1, N);   
T_static  = zeros(1, N);   
P_static = zeros(1, N);   
M        = zeros(1, N);   

for i = 1:N
    if x(i) <= 0
        gamma_local(i) = interp1(AR_eq_tab, gamma_eq_tab,  local_area_ratio(i), 'pchip', 'extrap');
        cp_local(i)  = interp1(AR_eq_tab, cp_gas_eq_tab,   local_area_ratio(i), 'pchip', 'extrap');
        Pr_local(i)  = interp1(AR_eq_tab, Pr_eq_tab,   local_area_ratio(i), 'pchip', 'extrap');
        mu_local(i)  = interp1(AR_eq_tab, mu_eq_tab,   local_area_ratio(i), 'pchip', 'extrap');
        T_static(i) = interp1(AR_eq_tab, T_eq_tab,    local_area_ratio(i), 'pchip', 'extrap');
        P_static(i) = interp1(AR_eq_tab, P_eq_tab,    local_area_ratio(i), 'pchip', 'extrap');
        M(i)       = interp1(AR_eq_tab, M_eq_tab, local_area_ratio(i), 'pchip', 'extrap');
    else
        gamma_local(i) = interp1(AR_fr_tab, gamma_fr_tab,  local_area_ratio(i), 'pchip', 'extrap');
        cp_local(i)  = interp1(AR_fr_tab, cp_gas_fr_tab,   local_area_ratio(i), 'pchip', 'extrap');
        Pr_local(i)  = interp1(AR_fr_tab, Pr_fr_tab,   local_area_ratio(i), 'pchip', 'extrap');
        mu_local(i)  = interp1(AR_fr_tab, mu_fr_tab,   local_area_ratio(i), 'pchip', 'extrap');
        T_static(i) = interp1(AR_fr_tab, T_fr_tab,    local_area_ratio(i), 'pchip', 'extrap');
        P_static(i) = interp1(AR_fr_tab, P_fr_tab,    local_area_ratio(i), 'pchip', 'extrap');
        M(i)       = interp1(AR_fr_tab, M_fr_tab, local_area_ratio(i), 'pchip', 'extrap');
    end
end
end

function [mu_w, Pr_w] = mu_w_to_temperature(T)
% This function computates the mu and the Pr of water depending on its
% temperature.

% INPUTS:
% T is the temperature of water in the section

% OUTPUTS:
% mu_w & Pr_w calculated at the input temperature

A = 2.414e-5;
B = 247.8;
C = 140;

cp_w   = 4183;           
k_w    = 0.654;

mu_w = A * 10^(B/(T-C));
Pr_w = mu_w * cp_w / k_w;
end

function [n_conv, n_div, x_conv, x_div, x, local_area_ratio, local_diameter, dx_vec, dr_vec, ds, ds_fixed, r_avg] = nozzle_geometry_section_array_plot( ...
    N, L_conv, L_div, L_tot, throat_radius, throat_area, a_conv_rad, a_div_rad)
% This function calculate the arrays for the sections iterations.

% INPUTS:
% N for number of sections
% L_conv & L_div are lenghts in m for convergent and divergent parts
% L_tot is total lenght in m of the nozzle (2 to 2 ratio) 

% OUTPUTS:
% n_conv is number of nodes in covergent part
% n_div is number of nodes in divergent part
% x_conv & x_div are coordinates for -L_conv | 0 | L_div
% x is total coordinates array with throat not duplicated
% dx_vec is axial spacing
% dr_vec is radial spacing
% ds is Arc-length elements (N-1)
% r_avg is average radius per element

n_conv = max(3, round(N * L_conv / L_tot));   
n_div  = N - n_conv + 1;                      
x_conv = linspace(-L_conv, 0,     n_conv);    
x_div  = linspace(0,       L_div, n_div);     
x      = [x_conv, x_div(2:end)];              

r = zeros(1, N);
for i = 1:N
    if x(i) <= 0
        r(i) = throat_radius + abs(x(i)) * tan(a_conv_rad);
    else
        r(i) = throat_radius + x(i) * tan(a_div_rad);
    end
end
local_area_ratio     = (r / throat_radius).^2;  
local_diameter = 2 * sqrt(((local_area_ratio * throat_area) / pi));

dx_vec = diff(x);                              
dr_vec = diff(r);                              
ds     = sqrt(dx_vec.^2 + dr_vec.^2);   
ds_fixed = L_tot / N;
r_avg  = (r(1:end-1) + r(2:end)) / 2;  
end

function [hc_no_TBC, t_YSZ, R_YSZ] = tbc_thickness_calculation_plot(throat_idx, ...
    cp_w, mdot_w_hp, T_in, T_lim_wall, T0, P0, ...
    C_star, local_area_ratio, local_diameter, ds, gamma_local, ...
    throat_diameter, throat_curvature_radius, mu_0_fr, cp_local, Pr_local, t_wall, k_wall, k_YSZ, ...
    T_lim_YSZ, t_bond, k_bond, M, N)

% This function calculates first the q and the water temperature in the
% throat section without the TBC. 
% Then it calculates the thickness for the TBC to be used for the next
% computations.

% INPUTS:
% throat_idx [int] Array index corresponding to the throat section
% cp_w [J/(kg K)] Specific heat of water
% mdot_w_hp [kg/s] Design mass flow rate of water in the cooling circuit
% T_in [K] Water inlet temperature
% T_lim_wall [K] Maximum allowable structural temperature for the metallic wall 
% T0 [K] Stagnation temperature of the gases in the combustion chamber
% P0 [Pa] Stagnation pressure in the combustion chamber
% C_star [m/s] Characteristic velocity of the propellant
% local_area_ratio [array 1xN] Expansion ratio along the nozzle axis
% local_diameter [array 1xN] Local nozzle diameter [m] along the axis
% ds [array 1xN-1] Length of each single axial segment [m]
% gamma_local [array 1xN] Specific heat ratio of the exhaust gases along the axis
% throat_diameter [m] Geometric diameter of the throat
% throat_curvature_radius [m] Radius of curvature of the throat region
% mu_0_fr [Pa*s] Viscosity of the hot gases calculated at film temperature
% cp_local [array 1xN] Specific heat at constant pressure of the gases [J/(kg K)]
% Pr_local [array 1xN] Local Prandtl number of the exhaust gases
% t_wall [m] Thickness of the inner metallic wall (copper)
% k_wall [W/(m K)] Thermal conductivity of the copper wall
% k_YSZ [W/(m K)] Equivalent thermal conductivity of the ceramic 
% T_lim_YSZ [K] Maximum allowable surface temperature of YSZ
% t_bond [m] Thickness of the metallic bond coat 
% k_bond [W/(m K)] Thermal conductivity of the bond coat
% M [array 1xN] Local Mach number of the gases along the axis
% N [int] Total number of nodes for spatial discretization
%
% OUTPUTS:
% hc_no_TBC [W/(m^2 K)] Water side convective heat transfer coefficient 
% t_YSZ [m] Calculated YSZ thickness 
% R_YSZ [K m^2/W] Total thermal resistance of YSZ

T_initial_no_TBC = T_in;
T_initial_TBC = T_in;

T_out_no_TBC_vec = [];
hg_no_TBC_vec = [];

T_out_TBC_vec = [];
hg_TBC_vec = [];


% Heat without TBC
for j = 1:N-1

xi = 1 + (gamma_local(j) - 1)/2 * M(j)^2;
sigma_no_TBC = ( (T_lim_wall / (2 * T0)) * xi + 0.5 )^(-0.68) * (xi)^(-0.12);

    hg_no_TBC = 0.026/(throat_diameter.^0.2) * ((throat_diameter/throat_curvature_radius).^0.1) * (((mu_0_fr.^0.2) * cp_local(j))/(Pr_local(j).^0.6)) * ((P0/C_star).^0.8) * ((local_area_ratio(j)).^-0.9) * sigma_no_TBC;

    hg_no_TBC_vec(j) = hg_no_TBC;

    R_no_TBC = (1 + ((Pr_local(j).^(1/3)) * ((gamma_local(j) - 1)/2) * (M(j).^2))) / (1 + ((gamma_local(j) - 1)/2) * (M(j).^2));
    T_adiabatic_no_TBC = R_no_TBC * T0;

    q_no_TBC = hg_no_TBC * (T_adiabatic_no_TBC - T_lim_wall);

    T_out_no_TBC = ((q_no_TBC * pi * local_diameter(j) * ds(j)) / (mdot_w_hp * cp_w)) + T_initial_no_TBC;
    T_initial_no_TBC = T_out_no_TBC; 

        T_out_no_TBC_vec(j) = T_out_no_TBC;
    
    if j == throat_idx
       q_throat_no_TBC = q_no_TBC;
       hg_throat_no_TBC = hg_no_TBC;
       T_water_throat_noTBC = T_out_no_TBC;
    end

    dT_wall_vec(j) = q_no_TBC * (t_wall / k_wall);
hc_no_TBC_vec(j) = q_no_TBC / ((T_lim_wall - dT_wall_vec(j)) - T_out_no_TBC);
end

% Thermal jump for wall
dT_wall = q_throat_no_TBC * (t_wall / k_wall);
hc_no_TBC = q_throat_no_TBC / ((T_lim_wall - dT_wall) - T_water_throat_noTBC);




% Marching with TBC
for j = 1:N-1

xi = 1 + (gamma_local(j) - 1)/2 * M(j)^2;
sigma_TBC = ( (T_lim_YSZ / (2 * T0)) * xi + 0.5 )^(-0.68) * (xi)^(-0.12);

    hg_TBC = 0.026/(throat_diameter.^0.2) * ((throat_diameter/throat_curvature_radius).^0.1) * (((mu_0_fr.^0.2) * cp_local(j))/(Pr_local(j).^0.6)) * ((P0/C_star).^0.8) * (local_area_ratio(j).^-0.9) * sigma_TBC;

    hg_TBC_vec(j) = hg_TBC;

    R_YSZ = (1 + ((Pr_local(j).^(1/3)) * ((gamma_local(j) - 1)/2) * (M(j).^2))) / (1 + ((gamma_local(j) - 1)/2) * (M(j).^2));
    T_adiabatic_TBC = R_YSZ * T0;
    
    q_TBC = hg_TBC * (T_adiabatic_TBC - T_lim_YSZ);

    T_out_TBC = ((q_TBC * pi * local_diameter(j) * ds(j)) / (mdot_w_hp * cp_w)) + T_initial_TBC;

    T_out_TBC_vec(j) = T_out_TBC;

    T_initial_TBC = T_out_TBC;
    
    if j == throat_idx
        T_water_throat_TBC = T_out_TBC;
        hg_throat_TBC = hg_TBC;
        T_aw_throat = T_adiabatic_TBC;
    end
end


% --- TBC thickness --- %
hc_max = 40000.00;
q_new = (hc_max * (T_lim_wall - T_water_throat_TBC)) / (1 + hc_max * (t_wall/k_wall));
T_surface_actual = T_aw_throat - (q_new / hg_throat_TBC);
t_YSZ = k_YSZ * (((T_surface_actual - T_lim_wall) / q_new) - (t_bond / k_bond));
R_YSZ  = (t_YSZ  / k_YSZ);   


%  PLOTS AND REPORT

fprintf('\n================================================================\n');
fprintf('       SIZING TBC      \n');
fprintf('                                     \n');
fprintf('TBC thickness computed at throat:    %.3f mm\n', t_YSZ * 1000);
fprintf('================================================================\n\n');



figure('Name', 'TBC justification', 'Color', 'w', 'Position', [100, 100, 900, 500]);

x_coords = cumsum(ds(:)');
yyaxis left
plot(x_coords(1:N-1), hg_no_TBC_vec, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'h_g without TBC'); hold on;
plot(x_coords(1:N-1), hg_TBC_vec, 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'h_g with TBC');
ylabel('h_g [W/(m^2 K)]', 'FontWeight', 'bold', 'FontSize',12);
set(gca, 'YColor', [0.8 0 0], 'FontSize',16); % Colore rosso per hg

yyaxis right
plot(x_coords(1:N-1), T_out_no_TBC_vec, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'T water without TBC');
plot(x_coords(1:N-1), T_out_TBC_vec, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'T water with TBC');
ylabel('Water temperature [K]', 'FontWeight', 'bold','FontSize',16);
set(gca, 'YColor', [0 0 0.8]); % Colore blu per T

grid on; xlabel('Axial position x [m]', 'FontSize', 16);
legend('Location', 'northeastoutside');
xline(x_coords(throat_idx), '--k', 'LabelVerticalAlignment', 'bottom','DisplayName', 'Throat','LineWidth',3);

end

function [best_a, best_T_water, best_T_wc, v_plot, q_req, Margin_plot, N_ch_ottimale, b_ch_ottimale, mdot_ottimale, P_inlet_ottimale] = ...
    water_channel_sizing_dp_min_plot(N, T0, gamma_local, Pr_local, M, throat_diameter, throat_curvature_radius, mu_0_fr, cp_local, P0, C_star, local_area_ratio, T_lim_wall, ...
   t_YSZ, k_YSZ, t_bond, k_bond, t_wall, k_wall, local_diameter, ds, T_in, cp_w, k_w, rho_w, throat_idx)

% This function compute the geometry of the water channels of the thermal
% jacket. It first computes the heat flux the water must absorb and then
% all the valid configurations in a research space. Finally it evaluates
% the optimal (dP minimum) configuartion.

% INPUTS:
% N [int] Total number of nodes for spatial discretization
% T0 [K] Stagnation temperature of the gases in the combustion chamber
% gamma_local [array 1xN] Specific heat ratio of the exhaust gases along the axis
% Pr_local [array 1xN] Local Prandtl number of the exhaust gases
% M [array 1xN] Local Mach number of the gases along the axis
% throat_diameter [m] Geometric diameter of the throat
% throat_curvature_radius [m] Radius of curvature of the throat region
% mu_0_fr [Pa*s] Viscosity of the hot gases calculated at film temperature
% cp_local [array 1xN] Specific heat at constant pressure of the gases [J/(kg K)]
% P0 [Pa] Stagnation pressure in the combustion chamber
% C_star [m/s] Characteristic velocity (c-star) of the propellant
% local_area_ratio [array 1xN] Expansion ratio (A/A_t) along the nozzle axis
% T_lim_wall [K] Maximum allowable structural temperature for the metallic wall (copper)
% t_YSZ [m] Thickness of the YSZ
% k_YSZ [W/(m K)] Thermal conductivity of the YSZ
% t_bond [m] Thickness of the metallic bond coat
% k_bond [W/(m K)] Thermal conductivity of the bond coat
% t_wall [m] Thickness of the inner metallic wall (copper)
% k_wall [W/(m K)] Thermal conductivity of the copper wall
% local_diameter [array 1xN] Local nozzle diameter [m] along the axis
% ds [array 1xN-1] Length of each single axial segment [m]
% T_in [K] Water inlet temperature
% cp_w [J/(kg K)] Specific heat of water (coolant)
% k_w [W/(m K)] Thermal conductivity of water
% rho_w [kg/m^3] Density of water
% throat_idx [int] Array index corresponding to the throat section
%
% OUTPUTS:
% best_a [m] Optimal channel width found by the sizing algorithm
% best_T_water [array 1xN] Coolant bulk temperature profile [K] along the axis for the optimal configuration
% best_T_wc [array 1xN] Cold-side wall temperature profile [K] along the axis for the optimal configuration
% v_plot [array 1xN] Coolant velocity profile [m/s] along the axis for the optimal configuration
% q_req [array 1xN] Thermal heat flux profile [W/m^2] absorbed by the coolant
% Margin_plot [array 1xN] Thermal margin profile [K] representing the difference between water saturation temperature and cold wall temperature
% N_ch_ottimale [int] Optimal number of cooling channels selected by the algorithm
% b_ch_ottimale [m] Optimal channel height selected by the algorithm
% mdot_ottimale [kg/s] Optimal total coolant mass flow rate
% P_inlet_ottimale [Pa] Calculated required inlet pressure to overcome friction losses and satisfy thermal constraints
% P_out [Pa] Resulting outlet pressure of the coolant

% Initialization outputs
q_req = zeros(1, N);
best_a = zeros(1, N);
best_T_water = zeros(1, N);
best_T_wc = zeros(1, N-1);
v_plot = zeros(1, N-1);
Margin_plot = zeros(1, N-1);
N_ch_ottimale = 0;
b_ch_ottimale = 0;
mdot_ottimale = 0;
P_inlet_ottimale = 0;
best_hc = zeros(1,N-1);

T_lim_YSZ = 1673;

% Compute the thermal flux on the gas side
for j = 1:N
    R = (1 + ((Pr_local(j).^(1/3)) * ((gamma_local(j) - 1)/2) * (M(j).^2))) / (1 + ((gamma_local(j) - 1)/2) * (M(j).^2));
    T_adiabatic = R * T0;
    T_surf_guess = 0.7 * T_adiabatic;
    err_surf = 100;
    while err_surf > 0.1
        xi = 1 + (gamma_local(j) - 1)/2 * M(j)^2;
        sigma = ( (T_surf_guess / (2 * T0)) * xi + 0.5 )^(-0.68) * (xi)^(-0.12);
        hg = 0.026/(throat_diameter.^0.2) * ((throat_diameter/throat_curvature_radius).^0.1) * (((mu_0_fr.^0.2) * cp_local(j))/(Pr_local(j).^0.6)) * ((P0/C_star).^0.8) * (local_area_ratio(j).^-0.9) * sigma;
        q_comp = (T_adiabatic - T_lim_wall) / ((1/hg) + (t_YSZ/k_YSZ) + (t_bond/k_bond));
        T_surf_new = T_adiabatic - (q_comp / hg);
        err_surf = abs(T_surf_new - T_surf_guess);
        T_surf_guess = T_surf_new; 
    end
    if T_surf_guess > T_lim_YSZ
        T_surf_fixed = T_lim_YSZ;
        sigma_lim = ( (T_surf_fixed / (2 * T0)) * xi + 0.5 )^(-0.68) * (xi)^(-0.12);
        hg_lim = 0.026/(throat_diameter.^0.2) * ((throat_diameter/throat_curvature_radius).^0.1) * (((mu_0_fr.^0.2) * cp_local(j))/(Pr_local(j).^0.6)) * ((P0/C_star).^0.8) * (local_area_ratio(j).^-0.9) * sigma_lim;
        q_comp = hg_lim * (T_adiabatic - T_surf_fixed);
    end
    q_req(j) = q_comp;
end

% Compute total heat and throat heat (Q_tot e Q_to_throat)
Q_tot = 0; Q_to_throat = 0;
for j = 1:(N-1)
    circ_avg = pi * (local_diameter(j) + local_diameter(j+1)) / 2; 
    node_area = circ_avg * ds(j); 
    q_avg = (q_req(j) + q_req(j+1)) / 2;
    Q_tot = Q_tot + (q_avg * node_area);
    if j < throat_idx; Q_to_throat = Q_to_throat + (q_avg * node_area); end
end

% Pre filtering of research space
b_test = 0.001:0.0001:0.005; 
N_ch_test = 90:1:180; 
mdot_test = 2:0.1:8; 
P_out = 50; c_ch = 0.0015;            
P_sat_tab = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150];
T_sat_tab = [180, 212, 234, 250, 264, 275, 286, 295, 303, 311, 342] + 273.15; 
T_sat = interp1(P_sat_tab, T_sat_tab, P_out, 'linear');
T_out_target = T_sat - 50;
couple_passed = [];

% Cycle for each configuration in the research space
for mdot_water = mdot_test
    T_throat_tmp = T_in + (Q_to_throat / (mdot_water * cp_w));
    T_wall_cold_throat = T_lim_wall - q_req(throat_idx) * (t_wall/k_wall);
    hc_req = q_req(throat_idx) / (T_wall_cold_throat - T_throat_tmp);
    [mu_w_tmp, Pr_w_tmp] = mu_w_to_temperature(T_throat_tmp);
    for b_ch = b_test
        for N_ch = N_ch_test
            a_ch = 2 * pi * (throat_diameter/2 + t_wall + t_YSZ + t_bond) / (N_ch) - c_ch;
            if a_ch < 0.0008; continue; end
            v_w = (mdot_water / N_ch) / (rho_w * (a_ch * b_ch)); 
            D_h = 4 * (a_ch * b_ch) / (2 * (a_ch + b_ch));
            Re_w = (rho_w * v_w * D_h) / mu_w_tmp;
            hc_calc = (0.023 * (Re_w^0.8) * (Pr_w_tmp^0.4) * k_w) / D_h;

            if hc_calc >= hc_req
                couple_passed = [couple_passed; mdot_water, N_ch, b_ch]; 
            end 
        end
    end
end

% Optimization cycle
num_configurazioni = size(couple_passed, 1);
best_dP_opt = inf; 
best_v_throat = 0;    
config_valide = 0;
best_Margin = inf;
Re_water = [];

% *Experimental: factors to scale inlet area of the channels to have a
% variable number of channels and a fixed area (fattore = 1) or a varaible
% area from 1 (throat) to 1.2 (inlet)
%fattore_ingresso = 1; 
fattore_uscita = 1;
%fattori_ingresso = linspace(fattore_ingresso, 1.0, throat_idx);
%fattori_uscita = linspace(1, fattore_uscita, N-throat_idx);

% Testing each valid configurations for the boiling margin
for j = 1:num_configurazioni
    mdot_water = couple_passed(j, 1); N_ch = couple_passed(j, 2); b_ch = couple_passed(j, 3);
    a_local = zeros(1, N); dP_local = zeros(1, N-1);
    T_water_local = zeros(1, N); T_water_local(1) = T_in;
    T_wc_local = zeros(1, N-1); v_w_temp = zeros(1, N-1); Margin_temp = zeros(1, N-1);
    Margin_min = inf; fail_flag = false; 
    hc_throat_current = 0; v_throat_current = 0;

    for i = 1:N-1

        % *Experimental
        % if i <= throat_idx
        %     a_th = 2 * pi * (throat_diameter/2 + t_wall + t_TBC + t_bond) / (N_ch) - c_ch;
        %     a_local(i) = 2 * pi * (throat_diameter/2 + t_wall + t_TBC + t_bond) / (N_ch) - c_ch;
        %     a_local(i) = a_local(i) * fattori_ingresso(i);
        %     R_iniziale = local_diameter(1)/2 + t_wall + t_TBC + t_bond;
        %     N_iniziale = 2*pi*R_iniziale / (a_th + c_ch);
        %     N_iniziale = round(N_iniziale);   % perché deve essere intero
        %     if N_iniziale >= 210
        %         fail_flag = true;
        %         break;
        %     end
        % 
        % else        
        %     a_local(i) = 2 * pi * (local_diameter(i)/2 + t_wall + t_TBC + t_bond) / (N_ch) - c_ch;
        %     a_local(i) = a_local(i) * fattori_uscita(i - throat_idx);
        % end
        
        a_local(i) = 2 * pi * (local_diameter(i)/2 + t_wall + t_YSZ + t_bond) / (N_ch) - c_ch;
        ch_area = a_local(i) * b_ch; D_h = 4 * ch_area / (2 * (a_local(i) + b_ch));
        [mu_w_l, Pr_w_l] = mu_w_to_temperature(T_water_local(i));
        v_w = (mdot_water / N_ch) / (rho_w * ch_area); 
        v_w_temp(i) = v_w;

        Re_w = (rho_w * v_w * D_h) / mu_w_l;
        Re_water(i) = Re_w;
        hc_local = (0.023 * (Re_w^0.8) * (Pr_w_l^0.4) * k_w) / D_h;
        hc_vec(i) = hc_local;
        
        if i == throat_idx
            hc_throat_current = hc_local; 
            v_throat_current = v_w; 
        end

        T_water_local(i+1) = ((q_req(i) * pi * local_diameter(i) * ds(i)) / (mdot_water * cp_w)) + T_water_local(i);
        T_wc_vera = T_water_local(i+1) + (q_req(i) / hc_local);
        T_wc_local(i) = T_wc_vera;
        Margin_val = T_sat - T_wc_vera;
        Margin_temp(i) = Margin_val;
        
        if Margin_val < Margin_min; Margin_min = Margin_val; end
        if Margin_val <= 15; fail_flag = true; break; end 
        
        dP_local(i) = (0.316 * (Re_w^-0.25)) * (ds(i) / D_h) * 0.5 * rho_w * v_w^2;
    end

    a_local(N) = 2 * pi * (throat_diameter/2 + t_wall + t_YSZ + t_bond) / (N_ch) - c_ch;
    a_local(N) = a_local(N) * fattore_uscita;
    
    % Validating configurations
    if ~fail_flag && (T_out_target - T_water_local(end) > 50)
        config_valide = config_valide + 1;
        current_dP_total = sum(dP_local) / 1e5; 
        
        % If the one computed is the best, change it
        if (current_dP_total < best_dP_opt)
            best_dP_opt = current_dP_total;
            best_v_throat = v_throat_current; 
            best_hc_throat = hc_throat_current;
            
            best_a = a_local;
            best_T_water = T_water_local;
            best_T_wc = T_wc_local;
            v_plot = v_w_temp;
            Margin_plot = Margin_temp;
            mdot_ottimale = mdot_water;
            N_ch_ottimale = N_ch;
            b_ch_ottimale = b_ch;
            best_dP = best_dP_opt;
            P_inlet_ottimale = P_out + best_dP;
            best_Margin = Margin_min;
            min_Re = min(Re_water); 
            best_hc = best_hc_throat; 
            final_min_Re = min_Re;
        end
    end    
end

% --- REPORT FINALE A VIDEO ---
fprintf('\n======================================================================\n');
fprintf('             REPORT WATER CHANNELS          \n');
fprintf('             parameters used for the verification labeled with "!"\n other numbers are used for controlling and debugging\n          \n');

fprintf('======================================================================\n');
if config_valide > 0
    fprintf(' [SUCCES]\n');
    fprintf(' (Q_tot): %8.2f MW\n', Q_tot / 1e6);
    fprintf('!  (mdot_w): %8.3f kg/s\n', mdot_ottimale);
    fprintf('!  (N_ch): %8d\n', N_ch_ottimale); 

    fprintf('!  channel high (b_ch): %8.3f mm\n', b_ch_ottimale * 1000);
    fprintf('  Speed at throat (v_throat): %8.2f m/s\n', best_v_throat);
    fprintf('  total dP: %8.2f bar\n', best_dP);
    fprintf('!  inlet P: %8.2f bar\n', P_inlet_ottimale);
    fprintf('  minimum boiling margin: %8.1f °C\n', best_Margin);
    
    if final_min_Re < 10000
        fprintf('\n  [WARNING]: minimum Reynolds = %.0f. out of limit for Dittus-Boelter (>10000)!\n', final_min_Re);
    else
        fprintf('\n  Fluidodynamic: OK (minimum Reynolds = %.0f > 10000)\n', final_min_Re);
    end

else
    fprintf(' [FAILURE] - No acceptable configurations were found.\n');
end
fprintf('======================================================================\n\n');

end

function [k] = water_channel_sizing_dp_min_SS(N, T0, gamma_local, Pr_local, M, throat_diameter, throat_curvature_radius, mu_0_fr, cp_local, P0, C_star, local_area_ratio, T_lim_wall, ...
   t_TBC, k_TBC, t_bond, k_bond, t_wall, local_diameter, ds, T_in, cp_w, k_w, rho_w)

% This function repeats the computations of the sizing of the water
% channels in water_channel_sizing_dp_min_plot with a slightly different code to save some variables and arrays. The scope of this function
% is to generate the plot for the research space. This function was included to keep the main sizing function clear.

q_req = zeros(1, N);
T_lim_TBC = 1673;


for j = 1:N
    R = (1 + ((Pr_local(j).^(1/3)) * ((gamma_local(j) - 1)/2) * (M(j).^2))) / (1 + ((gamma_local(j) - 1)/2) * (M(j).^2));
    T_adiabatic = R * T0;
    T_surf_guess = 0.7 * T_adiabatic;
    err_surf = 100;
    
    while err_surf > 0.1
        xi = 1 + (gamma_local(j) - 1)/2 * M(j)^2;
        sigma = ( (T_surf_guess / (2 * T0)) * xi + 0.5 )^(-0.68) * (xi)^(-0.12);
        hg = 0.026/(throat_diameter.^0.2) * ((throat_diameter/throat_curvature_radius).^0.1) * (((mu_0_fr.^0.2) * cp_local(j))/(Pr_local(j).^0.6)) * ((P0/C_star).^0.8) * (local_area_ratio(j).^-0.9) * sigma;
        q_comp = (T_adiabatic - T_lim_wall) / ((1/hg) + (t_TBC/k_TBC) + (t_bond/k_bond));
        T_surf_new = T_adiabatic - (q_comp / hg);
        err_surf = abs(T_surf_new - T_surf_guess);
        T_surf_guess = T_surf_new; 
    end
    
    if T_surf_guess > T_lim_TBC
        T_surf_fixed = T_lim_TBC;
        sigma_lim = ( (T_surf_fixed / (2 * T0)) * xi + 0.5 )^(-0.68) * (xi)^(-0.12);
        hg_lim = 0.026/(throat_diameter.^0.2) * ((throat_diameter/throat_curvature_radius).^0.1) * (((mu_0_fr.^0.2) * cp_local(j))/(Pr_local(j).^0.6)) * ((P0/C_star).^0.8) * (local_area_ratio(j).^-0.9) * sigma_lim;
        q_comp = hg_lim * (T_adiabatic - T_surf_fixed);
    end
    
    q_req(j) = q_comp;
end


b_ch = 0.0016; 
N_ch_test = 90:1:180; 
mdot_test = 0.1:0.1:8.0; 
P_out = 50;                         
c_ch = 0.0015;            
P_sat_tab = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150];
T_sat_tab = [180, 212, 234, 250, 264, 275, 286, 295, 303, 311, 342] + 273.15; 
T_sat = interp1(P_sat_tab, T_sat_tab, P_out, 'linear');
T_out_target = T_sat - 50;

DP_MAP = zeros(length(mdot_test), length(N_ch_test)); 
FEASIBLE_MAP_THERMAL = zeros(length(mdot_test), length(N_ch_test)); 
FEASIBLE_MAP_GEOM = zeros(length(mdot_test), length(N_ch_test)); 

best_dP_opt = inf; 
config_valide = 0;

for m_idx = 1:length(mdot_test)
    mdot_water = mdot_test(m_idx);
    
    for n_idx = 1:length(N_ch_test)
        N_ch = N_ch_test(n_idx);
        
        dP_local = zeros(1, N-1);
        T_water_local = zeros(1, N);
        T_water_local(1) = T_in;
        Margin_min = inf;
        geom_valid = true;
        
        for i = 1:N-1
            a_local_i = 2 * pi * (local_diameter(i)/2 + t_wall + t_TBC + t_bond) / (N_ch) - c_ch;
            
         
            if a_local_i < 0.0005
                geom_valid = false;
                a_local_i = 0.0005; 
            end
            
            ch_area_local = a_local_i * b_ch;
            hydraulic_diameter = 4 * ch_area_local / (2 * (a_local_i + b_ch));
            
            [mu_w_local, Pr_w_local] = mu_w_to_temperature(T_water_local(i));
            v_w = (mdot_water / N_ch) / (rho_w * ch_area_local); 
            
            Re_w = (rho_w * v_w * hydraulic_diameter) / mu_w_local;
            Nu_w = 0.023 * (Re_w^0.8) * (Pr_w_local^0.4);
            hc_local = (Nu_w * k_w) / hydraulic_diameter;
            
            T_out_calc = ((q_req(i) * pi * local_diameter(i) * ds(i)) / (mdot_water * cp_w)) + T_water_local(i);
            T_water_local(i+1) = T_out_calc;
            
            T_wc_vera = T_out_calc + (q_req(i) / hc_local);
            Margin_tmp = T_sat - T_wc_vera;
            
            if Margin_tmp < Margin_min
                Margin_min = Margin_tmp;
            end
            
            f_fric = 0.316 * (Re_w^(-0.25));
            dP_local(i) = f_fric * (ds(i) / hydraulic_diameter) * 0.5 * rho_w * v_w^2;
        end
        
        current_dP_total = sum(dP_local) / 1e5; 
        DP_MAP(m_idx, n_idx) = current_dP_total;
        
    
        if geom_valid
            FEASIBLE_MAP_GEOM(m_idx, n_idx) = 1;
        end
        
        if Margin_min > 15 && (T_out_target - T_water_local(end) > 50)
            FEASIBLE_MAP_THERMAL(m_idx, n_idx) = 1; 
            
            if geom_valid && (current_dP_total < best_dP_opt)
                best_dP_opt = current_dP_total;
                mdot_ottimale = mdot_water;
                N_ch_ottimale = N_ch;
                b_ch_ottimale = b_ch;
                config_valide = config_valide + 1;
            end
        end
    end    
end
 
% =========================================================================
% PLOT
% =========================================================================
if config_valide > 0
    figure('Name', 'Design Search Space', 'Color', 'w', 'Position', [100, 100, 900, 600]);
    hold on; box on;
   
    [~, h_contour] = contourf(N_ch_test, mdot_test, DP_MAP, 100, 'LineStyle', 'none');
    colormap(turbo); 
    
    cb = colorbar;
    ylabel(cb, '\Delta P [bar]',  'Interpreter', 'tex','FontSize', 18, 'FontWeight', 'bold');
    
    [~, h_bound_th] = contour(N_ch_test, mdot_test, FEASIBLE_MAP_THERMAL, [0.5 0.5], 'r-', 'LineWidth', 3);
    
    [~, h_bound_geom] = contour(N_ch_test, mdot_test, FEASIBLE_MAP_GEOM, [0.5 0.5], 'k--', 'LineWidth', 2.5);
    
    h_opt = plot(N_ch_ottimale, mdot_ottimale, 'p', ...
             'MarkerSize', 15, ...
             'MarkerFaceColor', 'y', ...
             'MarkerEdgeColor', 'k', ...
             'LineWidth', 1);
             
    xlabel('Number of Channels (N_{ch})', 'Interpreter', 'tex', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Coolant Mass Flow Rate (m_{dot}) [kg/s]', 'Interpreter', 'tex', 'FontSize', 18, 'FontWeight', 'bold');
    title('Design Search Space', 'Interpreter', 'tex', 'FontSize', 18);
    
    legend([h_bound_th, h_bound_geom, h_opt], ...
           {'Thermal Limit Set', ...
            'Manufacturing Limit (Channel Width > 0.5 mm)', ...
            'Optimal Design Point'}, ...
           'Location', 'southwest', 'FontSize', 18);
           
    grid on;
    set(gca, 'GridAlpha', 0.3, 'Layer', 'top', 'FontSize', 18); 
    hold off;
end

k = 1;
end


