%% INPUT-OUTPUT MODEL ESTIMATION: C-PEPTIDE KINETICS IN 7 SUBJECTS 
% This script evaluates single-, bi-, and tri-compartmental pharmacokinetic 
% models using Ordinary Least Squares (LS) and Weighted Least Squares (WLS), 
% followed by population analysis (NAD, NPD, STS).

clc; clear; close all; 

%% ========================================================================= 
% 1. DATA IMPORTATION AND PREPROCESSING 
% =========================================================================
[t1, conc1] = textread('DatiCPsog1.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t2, conc2] = textread('DatiCPsog2.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t3, conc3] = textread('DatiCPsog3.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t4, conc4] = textread('DatiCPsog4.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t5, conc5] = textread('DatiCPsog5.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t6, conc6] = textread('DatiCPsog6.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t7, conc7] = textread('DatiCPsog7.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');

d = 49650; % Administered dose of synthetic C-peptide [pmol]
CV = 0.04; % Coefficient of variation of the measurement error

% Ask the user to select the subject to be analyzed.
soggetto = input('Enter patient number to process data (1-7): ');
num = num2str(soggetto);
file = strcat('DatiCPsog', num, '.dat');
disp(['Selected file ', file, '.'])

[t, conc] = textread(file, '%d%f', 'headerlines', 4, 'delimiter', 't');

% Subtract the baseline concentration from the selected subject. 
% The baseline corresponds to the first concentration measurement.
basale = conc(1);
conc = conc - basale;

% I assume that the baseline corresponds to the steady-state condition
% of the system before the experiment.
% The baseline can be subtracted because the model is linear and therefore
% the superposition principle applies.
% By subtracting the baseline, I obtain a difference model that describes
% the kinetics of the synthetic C-peptide only. Therefore, I focus on the
% system response to the experimental perturbation rather than on its
% equilibrium condition, where the baseline secretion is still present.

% Remove the measurement taken at minute 1, as it is considered unreliable
% according to the instructions.
conc(1:2) = [];
t(1:2) = [];

% Plot the measurements of the selected subject together with those of the
% other 6 subjects.
figure(1);
plot(t1, conc1); hold on;
plot(t2, conc2); hold on;
plot(t3, conc3); hold on;
plot(t4, conc4); hold on;
plot(t5, conc5); hold on;
plot(t6, conc6); hold on;
plot(t7, conc7); hold on;
plot(t, conc, 'r', 'LineWidth', 2); 
title('Dataset ($z_i$ measurements)');
xlabel('Time [min]');
ylabel('Concentration [pmol/mL]');

%% ========================================================================= 
% 2. ERROR VARIANCE AND INITIAL PARAMETER SETUP 
% =========================================================================
% C-peptide kinetics can be well-approximated by a mono-, bi-, or 
% tri-compartmental model. However, knowing the system model, we leverage 
% the analytical responses (inverse transforms) identified in Exercise 1 
% and consider the outputs as model predictions (y_i).

% Identify the variance, standard deviation, and covariance matrix 
% of the measurement error, given its known CV.
sigma = CV * (conc + basale); % Note: The CV applies to the entire measurement, baseline included!
sigma2 = sigma .^ 2;
SIGMA_V = diag(sigma2);        % Measurement error covariance matrix (v)

% Error measurements are typically uncorrelated, so the covariance matrix 
% is diagonal: v has zero mean and independent samples (uncorrelated measurements).

% For both LS and WLS, initial guesses must be defined for the unknown 
% parameters to be estimated (k01, k12, k21, V1). They can be chosen 
% arbitrarily, but they must match the physical units of the target parameters 
% to prevent optimization failures and avoid falling into local minima.

% Since these unknown parameters have known physiological values from 
% the literature, we set the initial guesses to those provided in Exercise 1:
p0 = [0.0654, 0.0568, 0.0716, 3290]; % [k01, k12, k21, V1] - Note: unit consistency!
options = optimset('Display', 'Off');

%% ========================================================================= 
% 3. INITIAL PARAMETER CONSISTENCY CHECK (P0 MODEL CHECK)
% =========================================================================
A0 = [-p0(1)-p0(3), p0(2); p0(3), -p0(2)];
B0 = [1; 0];
C0 = [1/p0(4), 0];
D0 = [0];
SDL0 = ss(A0, B0, C0, D0);
[y0, t0] = impulse(d * SDL0);

figure(2);
subplot(2, 1, 1);
plot(t0, log(y0));
hold on;
plot(t, log(conc), '*');
title('Model vs Data Comparison (Semilog Scale)');
legend('y', 'z');

subplot(2, 1, 2);
plot(t0, y0);
hold on;
plot(t, conc, '*');
legend('y', 'z');
title('Model vs Data Comparison (Natural Scale)');
xlabel('Time [min]');
ylabel('Concentration [pmol/L]');

% We can see that overall, the initial values are reasonable since they 
% are derived from literature. This assists the optimization algorithm 
% (lsqnonlin) and decreases the probability of getting trapped in a local minimum.

%% ========================================================================= 
% 4. LEAST SQUARES (LS) AND WEIGHTED LEAST SQUARES (WLS) ESTIMATION
% =========================================================================
% Ordinary Least Squares (LS) Estimation:
w = 0;
[phat, resnorm, ~, exitflag, output, lambda, S] = lsqnonlin(@(p) residui_new(p, conc, t, w), p0, [], [], options);
[residual, y] = residui_new(phat, conc, t, w);

figure(2);
subplot(2, 3, 1);
plot(t, conc, '.');
hold on;
plot(t, y, 'r');
title('LS Model');
legend('z', 'phat');

subplot(2, 3, 2);
plot(t, log(conc), '*');
hold on;
plot(t, log(y), 'r');
title('LS Model (Semilog Scale)');
legend('z', 'phat');

subplot(2, 3, 3);
plot(t, residual, '-ob');
yline(0, '--r');
xlabel('Time [h]');
title('LS Model Residuals');

% Evaluate LS goodness of fit:
SIGMA_p_ls = full(inv((S' * inv(SIGMA_V) * S)));
sigma_p_ls = sqrt(diag(SIGMA_p_ls));
CV_p_ls = 100 * (sigma_p_ls ./ phat');

% Weighted Least Squares (WLS) Estimation:
W = SIGMA_V; % A good choice for W is the inverse of the measurement error covariance matrix
w = sigma;   % Note: Since residui_new.m returns the residual rather than its squared value, 
             % we pass the standard deviation as weights. lsqnonlin squares them internally, 
             % correctly using the variance matrix instead of standard deviations!
             
[phat_pes, resnorm_pes, residual_pes, exitflag_pes, output_pes, lambda_pes, S_pes] = lsqnonlin(@(p) residui_new(p, conc, t, w), p0, [], [], options);
[residual_pes, y_pes] = residui_new(phat_pes, conc, t, w);

subplot(2, 3, 4);
plot(t, conc, '.');
hold on;
plot(t, y_pes, 'r');
title('WLS Model');
legend('z', 'phat');

subplot(2, 3, 5);
plot(t, log(conc), '*');
hold on;
plot(t, log(y_pes), 'r');
title('WLS Model (Semilog Scale)');
legend('z', 'phat');

subplot(2, 3, 6);
plot(t, residual_pes, '-ob');
yline(0, '--r');
xlabel('Time [h]');
title('WLS Model Residuals');

% Evaluate WLS goodness of fit:
SIGMA_p_wls = full(inv((S_pes' * inv(SIGMA_V) * S_pes)));
sigma_p_wls = sqrt(diag(SIGMA_p_wls));
CV_p_wls = 100 * (sigma_p_wls ./ phat_pes');

% Which approach to choose?
% We should use Weighted Least Squares (WLS) because the defined model has a 
% constant CV (CV = 4%). Analyzing the measurement variances, we notice they 
% are all different (it is not a constant standard deviation model), so we must 
% downweight data more heavily affected by noise -> this can only be achieved 
% using WLS! This is also evident from the fact that the parameter CVs calculated 
% via WLS are smaller than those calculated via LS!

%% ========================================================================= 
% 5. COMPARTMENTAL VS BI-EXPONENTIAL MODEL COMPARISON
% =========================================================================
% Now that the unknown parameters are estimated, define the system:
k01hat = phat_pes(1);
k12hat = phat_pes(2);
k21hat = phat_pes(3);
Vhat = phat_pes(4);
d = 49650; % Bolus of 49650 pmol defined in Exercise 1

A = [-k01hat - k21hat, k12hat; k21hat, -k12hat];
B = [1; 0];
C = [1/Vhat, 0];
D = [0];

SDLhat = ss(A, B * d, C, D);
FDThat = tf(SDLhat);
[numeratore, denominatore] = tfdata(FDThat); % tfdata extracts the numerator and denominator of the transfer function
[zeri, poli, k] = residue(numeratore{:}, denominatore{:});

% Using partial fraction expansion, we obtain a transfer function of the form:
% X1/(s + P1) + X2/(s + P2)
% Inverse transforming yields: Y(t) = X1*exp(-P1*t) + X2*exp(-P2*t)

% We can therefore numerically compare the exponential model outputs with 
% those of the compartmental model:
% Y(t) = X1*exp(-P1*t) + X2*exp(-P2*t) = A*exp(-alfa*t) + B*exp(-beta*t)
A = zeri(1);
alfa = poli(1);
B = zeri(2);
beta = poli(2);

% Compare predictions obtained from the bi-exponential model and the compartmental model
figure(4);
subplot(2, 1, 1);
[y_imp, t_imp] = impulse(SDLhat);
plot(t_imp, y_imp);
title('Compartmental Model Impulse Response');
xlabel('Time [min]');
ylabel('Concentration [pmol/L]');

subplot(2, 1, 2);
plot(t, A * exp(alfa * t) + B * exp(beta * t));
title('Bi-exponential Model Impulse Response');
xlabel('Time [min]');
ylabel('Concentration [pmol/L]');

disp(['Bi-exponential Inverse Transform: y = ', num2str(A), '*exp(', num2str(alfa), 't) + ', num2str(B), '*exp(', num2str(beta), 't)']);

%% =========================================================================
% 6. MONTE CARLO SIMULATION METHOD
% =========================================================================
% - Assuming the probability density function of error v is known:
%   knowing that the measurement error is typically an uncorrelated Gaussian noise 
%   with zero mean and variance sigma^2.
% - Synthetically generate a new dataset z = y(phat) + v.
% - Repeat the WLS estimation M times to estimate the N parameters multiple times 
%   and then compute their means, variances, and CVs.

for i = 1:50
    % z_mc = y_pes + normrnd(0, sigma);
    z_mc(:, i) = y_pes + normrnd(0, (CV * (conc + basale))); 
    [phat_MC, resnorm_mc, residual_mc, exitflag_mc, output_mc, lambda_mc, S_mc] = lsqnonlin(@(p) residui_new(p, z_mc(:, i), t, w), p0, [], [], options);
    phat_mc(i, :) = phat_MC;
end

medie_mc = mean(phat_mc);
varianze_mc = var(phat_mc);

% We notice that the CVs of the estimated parameters in the compartmental model case 
% are higher compared to the bi-exponential model case -> USING THE EXPONENTIAL MODEL 
% PARAMETRIZATION YIELDS LOWER CVs THAN THE COMPARTMENTAL MODEL PARAMETRIZATION.

SIGMAP_mc = cov(phat_mc); % Note: This is the exact SIGMA_P: represents the uncertainty of parameters estimated via Monte Carlo
CV_mc = 100 * (diag(sqrt(SIGMAP_mc)) ./ phat_MC');

figure(3);
subplot(2, 2, 1);
histogram(phat_mc(:, 1));
title('MC Histogram k01');
subplot(2, 2, 2);
histogram(phat_mc(:, 2));
title('MC Histogram k12');
subplot(2, 2, 3);
histogram(phat_mc(:, 3));
title('MC Histogram k21');
subplot(2, 2, 4);
histogram(phat_mc(:, 4));
title('MC Histogram V1');

% Since we have few data points (30), it's hard to clearly tell from individual 
% histograms whether the estimated parameter distributions are normal. 
% We apply the Kolmogorov-Smirnov test:
for i = 1:4
    pd = makedist('Normal', 'mu', medie_mc(i), 'sigma', sqrt(varianze_mc(i)));
    [h(i), p(i)] = kstest(phat_mc(:, i), pd); % p-value > 0.05 -> accept H0 -> populations are normally distributed
end

% Find Confidence Intervals
% From histograms and the K-S test, we note a normal distribution of estimated parameters.
% Therefore, find the 95% CIs as xbar +/- 2*std:
IC1 = medie_mc - 2 * sqrt(varianze_mc);
IC2 = medie_mc + 2 * sqrt(varianze_mc);

%% =========================================================================
% 7. BOOTSTRAP METHOD
% =========================================================================
% - We should randomly extract N sets of 30 data points each and compute 
%   a pointwise mean of the residuals. These will then be sampled with 
%   replacement to estimate z.
% - Normalize the weighted residuals because we have a constant CV model 
%   (and not a constant standard deviation model!):
residui_normalizzati = residual_pes; % Normalize residuals -> Note: lsqnonlin already outputs normalized residuals

% Since in the function passed to lsqnonlin we already defined normalized residuals as:
% res = (conc - y) ./ w, where w = sigma_v.

% - Shuffle the normalized residuals (randomly sampling with replacement) 
%   performing a total of 50 experiments:
for i = 1:50
    n = randi(length(conc), [length(conc), 1]);
    residui_estratti_normalizzati = residui_normalizzati(n); % Draw 30 random residuals with replacement
    
    % - Reconstruct measurements z:
    % Each residual is multiplied by its standard deviation (sigma) to extract 
    % a denormalized sample -> this is the extracted error which is then added to measurement y.
    z_bootstrap(:, i) = y_pes + sigma .* residui_estratti_normalizzati; % Denormalize weighted residuals
    
    [phat_BOOT, resnorm_boot, residual_boot, exitflag_boot, output_boot, lambda_boot, S_boot] = lsqnonlin(@(p) residui_new(p, z_bootstrap(:, i), t, w), p0, [], [], options);
    phat_boot(i, :) = phat_BOOT;   
end

SIGMAP_boot = cov(phat_boot); % Note: This is the exact SIGMA_P: represents the uncertainty of parameters estimated via Bootstrap

% Now we summarize the impulse responses found with the 2 methods (Monte Carlo and Bootstrap) 
% and the 2 adopted models (Bi-exponential and Compartmental)
figure(5);
subplot(2, 2, 1);
histogram(phat_boot(:, 1));
title('Bootstrap Histogram k01');
subplot(2, 2, 2);
histogram(phat_boot(:, 2));
title('Bootstrap Histogram k12');
subplot(2, 2, 3);
histogram(phat_boot(:, 3));
title('Bootstrap Histogram k21');
subplot(2, 2, 4);
histogram(phat_boot(:, 4));
title('Bootstrap Histogram V1');

figure(6);
plot(t, z_bootstrap(:, 1));
hold on;
plot(t, z_mc(:, 1));
xlabel('Time [min]');
ylabel('Concentration [pmol/L]');
title('Impulse Response: Monte Carlo vs Bootstrap Methods');
legend('Bootstrap', 'Monte Carlo');

CV_boot = (diag(sqrt(SIGMAP_boot)) ./ phat_BOOT') * 100;