%% ========================================================================
%  UNCERTAINTY ANALYSIS
%  ========================================================================
    
    Ncases = 5000; %Number of analyzed cases in the uncertainty analysis
    c_star = 1542; % This quantity has been calculated using NASA CEA and it's assumed to be constant
    a_nom = 2.0006e-05; % Burn rate coefficient [m/s/Pa^n]
    n_nom = 0.38155; % Pressure exponent
    pc0_nom = 63.7729e5; % Initial pressure chamber pressure [Pa]
    L0_nom = 1.5706; % Initial length of the grain [m]
    R_nom = 0.4432; % External radius of the grain [m]
    Dt_nom = 0.1143; % Throat diameter [m]
    r0_nom = 0.241; % Internal radius dimension [m]
    wt_AP_nom = 84.18; % Weight percentage of ammonium perchlorate
    wt_HTPB_nom = 100-wt_AP_nom; % Weight percentage of HTPB
    tb = 1000; % Burning time [s] (Burning time is higher than usual because here it also important to see when the grain stops burning)
    dt = 0.01; % Integration step [s]
    AP_rho_nom = 1950; % AP density [kg/m^3]
    HTPB_rho_nom = 915; % HTPB density [kg/m^3]
    rho_nom = 1/(wt_AP_nom/(100*AP_rho_nom) + wt_HTPB_nom/(100*HTPB_rho_nom)); % Nominal density of propellant grain [kg/m^3]

    st_rho = 0.0033*rho_nom; % Density standard deviation
    st_L0 = 0.0047*L0_nom; % Length of the grain standard deviation
    st_Dt = 0.008*Dt_nom; % Throat diameter standard deviantion


for x = 1:Ncases
    % NOW WE CALCULATE THE PERTURBED QUANTITIES 
    [a,n] = a_n_casual_generator(a_nom,n_nom);
    pc0 = pc0_nom; % Initial pressure chamber pressure (unperturbed) [Pa]
    L0 = L0_nom + st_L0*randn; % Initial length of the grain (perturbed) [m]
    R = R_nom; % External radius of the grain (unperturbed)
    Dt = Dt_nom + st_Dt*randn; % Throat diameter (perturbed)
    At = pi*(Dt/2)^2; % Throat area (perturbed)
    rho = rho_nom + st_rho*randn; % Grain density (perturbed)
    r0 = r0_nom; % Internal radius of the grain (unperturbed)
    wt_AP = wt_AP_nom; % Initial weight percentage of AP (unperturbed)
    wt_HTPB = 100-wt_AP; % Initial weight percentage of HTPB (unperturbed)
    AP_rho = AP_rho_nom; % Initial density of AP (unperturbed)
    HTPB_rho = HTPB_rho_nom; % Initial density of HTPB (unperturbed)

    Vol = pi*(R^2 - r0^2)*L0; % Volume of the grain [m^3]
    Mass = rho*Vol; % Mass of the grain [kg]

    Ru = 8.314; % Universal gas constant [J/(K*mol)]
    rb0 = a*pc0^n; % Initial burning rate [m/s]
    Ab0 = 2*pi*r0*L0 + 2*pi*(R^2 - r0^2); % Initial burning surface [m^2]
    mp0 = rho*rb0*Ab0; % Initial mass flow rate [kg/s]
  
    % Vector preallocation
    N = floor(tb/dt) + 1;
    rb = zeros(1,N);
    Ab = zeros(1,N);
    L  = zeros(1,N);
    r  = zeros(1,N);
    pc = zeros(1,N);
    mp = zeros(1,N);
    time = zeros(1,N);

    i = 2;
    rb(1) = rb0;
    Ab(1) = Ab0;
    L(1) = L0;
    r(1) = r0;
    pc(1) = pc0;
    mp(1) = mp0;
    time(1) = 0;

    % The following for cycle performs the chamber pressure forward
    % integration using the quasi-steady model
    for t = dt:dt:tb
        time(i) = t;
        L(i) = L(i-1) - 2*rb(i-1)*dt; % Updated length of the grain [m]
        r(i) = r(i-1) + rb(i-1)*dt; % Updated internal radius of the grain [m]
        if L(i)<=0 || r(i)>=R % If this condition is satisfied, it means the grain stopped burning
            break 
        end
        Ab(i) = 2*pi*r(i)*L(i) + 2*pi*(R^2 - r(i)^2); % Updated burning surface [m^2]
        % We now assume coefficient c_star to be constant
        pc(i) = (rho*a*Ab(i)*c_star/At)^(1/(1-n)); % Updated chamber pressure [Pa]
        rb(i) = a*pc(i)^n; % Updated burning rate [m/s]
        mp(i) = rho*rb(i)*Ab(i); % Updated mass flow rate [kg/s]
        i = i+1;
    end

    %Cut all vector components that come after the time of the end of combustion of the grain
    time = time(1:i-1);
    rb   = rb(1:i-1);
    Ab   = Ab(1:i-1);
    L    = L(1:i-1);
    r    = r(1:i-1);
    pc   = pc(1:i-1);
    mp   = mp(1:i-1);

    pc_max_vector(x) = max(pc);
    tb_max_vector(x) = max(time);
end

pc_max = max(pc_max_vector); % Maximum chamber pressure computed [Pa]
tb_max = max(tb_max_vector); % Maximum burning time computed [s]
tb_min = min(tb_max_vector); % Minimum burning time computed [s]

mu_pc  = mean(pc_max_vector); % Maximum chamber pressure average value [Pa]
sig_pc = std(pc_max_vector); % Standard deviation of maximum chamber pressure [Pa]
rel_st_pc = 100*sig_pc/mu_pc; % Relative standard deviation of maximum chamber pressure
mu_tb  = mean(tb_max_vector); % Burning time average value [s]
sig_tb = std(tb_max_vector); % Standard deviation of burning time [s]
rel_st_tb = 100*sig_tb/mu_tb; % Relative standard deviation of burning time

fprintf('===== UNCERTAINTY ANALYSIS ============================\n\n');

fprintf('--- Chamber pressure statistics ---\n');
fprintf('  Pc_max      = %.4e Pa  = %.4f bar (Maximum pressure reached in chamber)\n', pc_max, pc_max/1e5);
fprintf('  mean(Pcmax) = %.4e Pa  = %.4f bar (Average maximum chamber pressure)\n', mu_pc, mu_pc/1e5);
fprintf('  std(Pcmax)  = %.4e Pa  = %.4f bar (Maximum chamber pressure standard deviation)\n', sig_pc, sig_pc/1e5);
fprintf('  rel_std(Pc) = %.4f %% (Maximum chamber pressure relative standard deviation)\n\n', rel_st_pc);

fprintf('--- Burning time statistics ---\n');
fprintf('  tb_max      = %.4f s (Maximum burning time)\n', tb_max);
fprintf('  tb_min      = %.4f s (Minimum burning time)\n', tb_min);
fprintf('  mean(tb)    = %.4f s (Average burning time)\n', mu_tb);
fprintf('  std(tb)     = %.4f s (Burning time standard deviation)\n', sig_tb);
fprintf('  rel_std(tb) = %.4f %% (Burning time relative standard deviation)\n\n', rel_st_tb);

fprintf('=============================================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================================================================
%  CUMULATIVE STATISTICS PLOTS
%  ========================================================================

Nmc = length(pc_max_vector); % Number of Monte Carlo simulations
iterations = 1:Nmc; % Iteration vector

% Vector preallocation
pc_cum_mean    = zeros(1,Nmc);
pc_cum_std     = zeros(1,Nmc);
pc_cum_rel_std = zeros(1,Nmc);

tb_cum_mean    = zeros(1,Nmc);
tb_cum_std     = zeros(1,Nmc);
tb_cum_rel_std = zeros(1,Nmc);

% Cumulative statistics calculation
for k = 1:Nmc

    % Cumulative mean and standard deviation of maximum chamber pressure
    pc_cum_mean(k) = mean(pc_max_vector(1:k));
    pc_cum_std(k)  = std(pc_max_vector(1:k));

    % Cumulative relative standard deviation of maximum chamber pressure
    pc_cum_rel_std(k) = 100*pc_cum_std(k)/pc_cum_mean(k);

    % Cumulative mean and standard deviation of burning time
    tb_cum_mean(k) = mean(tb_max_vector(1:k));
    tb_cum_std(k)  = std(tb_max_vector(1:k));

    % Cumulative relative standard deviation of burning time
    tb_cum_rel_std(k) = 100*tb_cum_std(k)/tb_cum_mean(k);

end

%% ========================================================================
%  CONVERGENCE MONITORING PLOTS: MEAN AND STANDARD DEVIATION
%  ========================================================================

figure

subplot(2,2,1)
plot(iterations, pc_cum_mean/1e5)
grid on
title("Convergence of Pc_{max} mean")
xlabel("Number of simulations")
ylabel("Cumulative mean [bar]")

subplot(2,2,2)
plot(iterations, pc_cum_std/1e5)
grid on
title("Convergence of Pc_{max} standard deviation")
xlabel("Number of simulations")
ylabel("Cumulative standard deviation [bar]")

subplot(2,2,3)
plot(iterations, tb_cum_mean)
grid on
title("Convergence of burning time mean")
xlabel("Number of simulations")
ylabel("Cumulative mean [s]")

subplot(2,2,4)
plot(iterations, tb_cum_std)
grid on
title("Convergence of burning time standard deviation")
xlabel("Number of simulations")
ylabel("Cumulative standard deviation [s]")

%% ========================================================================
%  ADDITIONAL RELATIVE STANDARD DEVIATION PLOTS
%  ========================================================================

figure

subplot(1,2,1)
plot(iterations, pc_cum_rel_std)
grid on
title("Convergence of Pc_{max} relative standard deviation")
xlabel("Number of simulations")
ylabel("Cumulative relative standard deviation [%]")

subplot(1,2,2)
plot(iterations, tb_cum_rel_std)
grid on
title("Convergence of burning time relative standard deviation")
xlabel("Number of simulations")
ylabel("Cumulative relative standard deviation [%]")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,2,1)
h1 = plot([1:Ncases], pc_max_vector/10^5);
hold on
h2 = yline(mu_pc/10^5,'r','LineWidth',2);
grid on
title("Monte Carlo analysis of chamber pressure")
xlabel("Number of simulations")
ylabel("Maximum pressure [Bar] for each iteration")
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend(h2,' Pc mean value')
subplot(1,2,2)
h3 = plot([1:Ncases], tb_max_vector);
grid on
h4 = yline(mu_tb,'r','LineWidth',2);
title("Monte Carlo analysis of burning time")
xlabel("Number of simulations")
ylabel("Burning time [s] for each iteration")
h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend(h4,' tb mean value')





function [a,n] = a_n_casual_generator(a_nom,n_nom)

p  = (10^5)*[10.1 9.7 10.3 30.2 31.0 29.8 50.2 51.0 50.3 70.2 69.0 69.8 91.2 89.1 89.0]'; %Chamber pressure dataset [Pa]
rb = (10^-3)*[4.0  3.8 4.1  5.6  6.0  5.7  7.0  7.2  7.1  8.4  8.3  8.6  8.8  9.0  9.2 ]'; %Burning rate dataset [m/s]

% Linearization: log(rb) = log(a) + n*log(p)
x = log(p); % Independent variable of the linear regression
y = log(rb); % Dependent variable of the linear regression

% Matrix of the model: y = b + n*x
X = [ones(length(x),1), x]; % Regression matrix
 
% Estimation of [b; n] parameters
beta_hat = X \ y; % Solves the overdetermined linear system y = X*beta in the least-squares sense.

b_hat = beta_hat(1); % Extracts the first estimated parameter b_hat = log(a_hat)
n_hat = beta_hat(2); % Extracts the second estimated parameter n_hat

a_hat = exp(b_hat); % Since b_hat = log(a_hat), the coefficient a_hat is obtained by exponentiation

y_hat = X*beta_hat; % Computes the values of log(rb) predicted by the linear regression model
res   = y - y_hat; % Computes the difference between the experimental data and the values predicted by the regression

N = length(y);
m = 2;  % Parameters number: b e n
sigma2_hat = (res' * res) / (N - m); % Computes the unbiased estimate of the residual variance


Cov_beta = sigma2_hat * inv(X' * X); % Covariance matrix of the estimated parameters [b, n]

var_b  = Cov_beta(1,1); % Variance of the parameter b = log(a)
var_n  = Cov_beta(2,2); % Variance of the parameter n
cov_bn = Cov_beta(1,2); % Covariance between b = log(a) and n

sigma_b = sqrt(var_b); % Standard deviation of b = log(a)
sigma_n = sqrt(var_n); % Standard deviation of n


b_nom = log(a_nom);
mu = [b_nom, n_nom]; % Mean vector of the multivariate normal distribution
x = mvnrnd(mu, Cov_beta); % Generates one random sample [b_random, n_random] from a multivariate normal distribution.
a = exp(x(1)); % Converts the first randomly generated parameter from the logarithmic domain back to the physical domain
n = x(2); % Assigns to n the second randomly generated parameter
end
