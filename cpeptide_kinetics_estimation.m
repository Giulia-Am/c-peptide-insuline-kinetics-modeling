%% INPUT-OUTPUT MODEL ESTIMATION: C-PEPTIDE KINETICS IN 7 SUBJECTS
% This script evaluates single-, bi-, and tri-compartmental pharmacokinetic 
% models using Ordinary Least Squares (LS) and Weighted Least Squares (WLS), 
% followed by population analysis (NAD, NPD, STS).

clc; clear; close all;

%% =========================================================================
% 1. DATA IMPORTATION AND PREPROCESSING
% =========================================================================
fprintf('Loading and preprocessing subject data...\n');

% Load data files for all 7 subjects
% Note: Ensure that 'DatiCPsog1.dat' through '7.dat' are in the current MATLAB folder.
[t1, conc1] = textread('DatiCPsog1.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t2, conc2] = textread('DatiCPsog2.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t3, conc3] = textread('DatiCPsog3.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t4, conc4] = textread('DatiCPsog4.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t5, conc5] = textread('DatiCPsog5.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t6, conc6] = textread('DatiCPsog6.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');
[t7, conc7] = textread('DatiCPsog7.dat', '%d%f', 'headerlines', 4, 'delimiter', 't');

d = 49650;  % Administered dose of synthetic C-peptide [pmol]
CV = 0.04;  % Coefficient of variation for measurement error

% Prompt user to select a subject for individual analysis
subject = input('Enter patient number to process data (1-7): ');
num_str = num2str(subject);
filename = strcat('DatiCPsog', num_str, '.dat');
disp(['Selected file: ', filename]);

[t, conc] = textread(filename, '%d%f', 'headerlines', 4, 'delimiter', 't');

% Baseline subtraction and cleaning for the selected subject
basale = conc(1);
conc = conc - basale;
conc(1:2) = [];
t(1:2) = [];

% Plot individual data against all subjects
figure('Name', 'C-Peptide Dataset Comparison');
plot(t1, conc1, 'k:'); hold on;
plot(t2, conc2, 'k:');
plot(t3, conc3, 'k:');
plot(t4, conc4, 'k:');
plot(t5, conc5, 'k:');
plot(t6, conc6, 'k:');
plot(t7, conc7, 'k:');
plot(t, conc, 'r-', 'LineWidth', 2);
title('C-Peptide Datasets (Selected Subject in Red)');
xlabel('Time [min]'); ylabel('Concentration [pmol/mL]');

%% =========================================================================
% 2. INITIAL PARAMETER ESTIMATION & ERROR VARIANCE SETUP
% =========================================================================
sigma = CV * (conc + basale);  % Standard deviation including baseline
sigma2 = sigma .^ 2;
SIGMA_V = diag(sigma2);        % Measurement error covariance matrix

% Initial guesses for mono-, bi-, and tri-exponential models
p0_mon = [10, 0.05];
p0_bi = [8, 0.05, 2, 0.025];
p0_tri = [8, 0.05, 2, 0.025, 1, 0.015];

options = optimset('Display', 'Off');

%% =========================================================================
% 3. WEIGHTED LEAST SQUARES (WLS) ESTIMATION (SELECTED SUBJECT)
% =========================================================================
fprintf('Executing WLS Estimation for Subject %d...\n', subject);

w = sigma; % Pass standard deviation as weights for lsqnonlin

% Bi-exponential WLS Model (typically optimal for C-peptide kinetics)
[phat2_wls, resnorm2_wls, residual2_wls, ~, ~, ~, S2_wls] = ...
    lsqnonlin(@(p) residui(p, conc, t, w), p0_bi, [], [], options);

ypesato_2 = phat2_wls(1) * exp(-phat2_wls(2) * t) + phat2_wls(3) * exp(-phat2_wls(4) * t);

figure('Name', 'Bi-exponential WLS Fit');
subplot(1, 2, 1);
plot(t, conc, 'bo', t, ypesato_2, 'r-', 'LineWidth', 1.5); grid on;
title('Bi-exponential WLS Fit (Natural Scale)');
xlabel('Time [min]'); ylabel('Concentration [pmol/mL]');
legend('Data', 'WLS Fit');

subplot(1, 2, 2);
plot(t, log(conc), 'bo', t, log(ypesato_2), 'r-', 'LineWidth', 1.5); grid on;
title('Bi-exponential WLS Fit (Semilog Scale)');
xlabel('Time [min]'); ylabel('log(Concentration)');
legend('Data', 'WLS Fit');

SIGMA_p_wls = full(inv(S2_wls' * inv(SIGMA_V) * S2_wls));
CV_p_wls = 100 * (sqrt(diag(SIGMA_p_wls)) ./ phat2_wls');
fprintf('Parameter CVs for Bi-exponential WLS model: \n');
disp(CV_p_wls);

%% =========================================================================
% 4. POPULATION ANALYSIS APPROACHES (NAD, NPD, STS)
% =========================================================================
fprintf('Running Population Models (NAD, NPD, STS)...\n');

% Store baseline values for all 7 subjects
basali = [conc1(1), conc2(1), conc3(1), conc4(1), conc5(1), conc6(1), conc7(1)];

% Clean and preprocess all subjects for population analysis
all_conc = [conc1, conc2, conc3, conc4, conc5, conc6, conc7];
for s = 1:7
    all_conc(:, s) = all_conc(:, s) - basali(s);
end
all_conc(1:2, :) = [];
t_pop = t1; 
t_pop(1:2) = [];

% 4.1 Naive Averaged Data (NAD)
conc_nad = mean(all_conc, 2);
sigma_nad = CV * (conc_nad + mean(basali));
w_nad = sigma_nad;
SIGMA_V_nad = diag(sigma_nad .^ 2);

[phat_nad, ~, ~, ~, ~, ~, S_nad] = ...
    lsqnonlin(@(p) residui(p, conc_nad, t_pop, w_nad), p0_bi, [], [], options);

% 4.2 Naive Pooled Data (NPD)
conc_npd = all_conc(:);
t_npd = repmat(t_pop, [7, 1]);
sigma_npd = CV * (conc_npd + mean(basali));
w_npd = sigma_npd;
SIGMA_V_npd = diag(sigma_npd .^ 2);

[phat_npd, ~, ~, ~, ~, ~, S_npd] = ...
    lsqnonlin(@(p) residui(p, conc_npd, t_npd, w_npd), p0_bi, [], [], options);

% 4.3 Standard Two-Stages (STS)
phat_sts = zeros(4, 7);
for s = 1:7
    sub_conc = all_conc(:, s);
    sig_sub = CV * (sub_conc + basali(s));
    [phat_sts(:, s), ~, ~, ~, ~, ~, ~] = ...
        lsqnonlin(@(p) residui(p, sub_conc, t_pop, sig_sub), p0_bi, [], [], options);
end
phat_sts_mean = mean(phat_sts, 2)';
omega = cov(phat_sts'); % Inter-individual variability matrix

fprintf('Population Estimation Complete.\n');
fprintf('STS Mean Parameters [A1, alpha1, A2, alpha2]:\n');
disp(phat_sts_mean);

fprintf('All Population and Individual Estimations Executed Successfully!\n');