%% PHARMACOKINETIC MODELS ANALYSIS AND SIMULATION (LTI SYSTEMS)
% This script integrates 5 pharmacokinetic models modeled as LTI systems:
% 1. One-Compartment Model (IV Bolus)
% 2. One-Compartment Model with Extravascular Absorption
% 3. Two-Compartment Model (C-Peptide) with Optimal Experiment Design
% 4. Three-Compartment Model (Pole-Zero Cancellation Analysis)
% 5. Repeated Dosing Regimen (Pulse Train Simulation)

clc; clear; close all;

%% =========================================================================
% 1. ONE-COMPARTMENT MODEL (IV BOLUS)
% =========================================================================
fprintf('Running Section 1: One-Compartment Model...\n');

V1 = 5;      % Volume [L]
k01 = 1.2;   % Elimination rate constant [1/h]
d = 500;     % Dose [mg]

% LTI System Matrices
A1 = -k01; B1 = 1; C1 = 1/V1; D1 = 0;
sys1 = ss(A1, B1, C1, D1);
G1 = tf(sys1);

% Partial Fraction Expansion (Residue)
[r1, p1, ~] = residue(G1.Numerator{1}, G1.Denominator{1});

% Bode Diagram
figure('Name', 'Bode Diagram - One-Compartment Model');
bode(G1); grid on;
title('Bode Diagram - One-Compartment Model');

% Impulsive Response
figure('Name', 'Impulse Response - One-Compartment Model');
subplot(2,1,1);
[Y1, T1] = impulse(sys1);
plot(T1, Y1 * d, 'LineWidth', 1.5); grid on;
title('Impulse Response (Concentration vs Time)');
xlabel('Time [h]'); ylabel('Concentration [mg/L]');

subplot(2,1,2);
semilogy(T1, Y1 * d, 'LineWidth', 1.5); grid on;
title('Impulse Response (Semilogarithmic Scale)');
xlabel('Time [h]'); ylabel('Concentration [mg/L]');

% Analytical vs Numerical Solution
y1_analytical = (d / V1) * exp(-k01 * T1);
figure('Name', 'Analytical vs MATLAB Solution');
plot(T1, y1_analytical, 'b-', 'LineWidth', 2); hold on;
plot(T1, d * Y1, 'r--', 'LineWidth', 1.5); grid on;
title('Analytical vs MATLAB Impulse Response');
xlabel('Time [h]'); ylabel('Concentration [mg/L]');
legend('Analytical', 'MATLAB Numerical');

% Half-Life and Parameter Calculations
half_life1 = log(2) / k01;
xline(half_life1, 'k--', 'Half-Life', 'LineWidth', 1.2);

[max_conc1, idx_max1] = max(Y1 * d);
t_max_conc1 = T1(idx_max1);
tau1 = 1 / k01;
t_elim1 = 5 * tau1; % Elimination time (~5*tau)

% Non-Compartmental Parameters (AUC, AUMC, MRT, Clearance)
AUC1 = trapz(T1, d * Y1);
AUMC1 = trapz(T1, T1 .* (d * Y1));
MRT1 = AUMC1 / AUC1;
CL1 = d / AUC1;

% Sensitivity Analysis: Varying k01 and V1
k_coeffs = [0.01, 0.1, 0.5, 1, 2, 5, 10];
figure('Name', 'Impulse Response vs k01');
for i = 1:length(k_coeffs)
    sys_k = ss(-k01 * k_coeffs(i), 1, 1/V1, 0);
    [Y_k, T_k] = impulse(d * sys_k);
    plot(T_k, Y_k, 'LineWidth', 1.2); hold on;
end
grid on; title('Impulse Response for Different k01 Values');
xlabel('Time [h]'); ylabel('Concentration [mg/L]');

%% =========================================================================
% 2. ONE-COMPARTMENT MODEL WITH EXTRAVASCULAR ABSORPTION
% =========================================================================
fprintf('Running Section 2: One-Compartment Model with Absorption...\n');

V2 = 5; k01_abs = 1.2; k02 = 1.2; k21 = 2.2; d = 500;
ka = k01_abs + k21;     % Absorption rate constant
F = k21 / ka;           % Bioavailable fraction

% State-Space Model
A2 = [-k01_abs-k21, 0; k21, -k02];
B2_oral = [1; 0];
C2 = [0, 1/V2];
D2 = 0;

sys2_oral = ss(A2, B2_oral, C2, D2);
[Y2_oral, T2_oral] = impulse(sys2_oral);

% IV Reference Model (for Non-Compartmental Bioavailability)
B2_iv = [0; 1];
sys2_iv = ss(A2, B2_iv, C2, D2);
[Y2_iv, T2_iv] = impulse(sys2_iv, T2_oral);

% AUC and Bioavailability Calculation
AUC_oral = trapz(T2_oral, d * Y2_oral);
AUC_iv = trapz(T2_iv, d * Y2_iv);
F_noncomp = AUC_oral / AUC_iv;

figure('Name', 'Oral vs IV Impulse Response');
plot(T2_oral, d * Y2_oral, 'b-', 'LineWidth', 1.5); hold on;
plot(T2_iv, d * Y2_iv, 'r--', 'LineWidth', 1.5); grid on;
title('Oral (Extravascular) vs IV Impulse Response');
xlabel('Time [h]'); ylabel('Concentration [mg/L]');
legend('Oral Administration', 'IV Administration');

%% =========================================================================
% 3. TWO-COMPARTMENT MODEL (C-PEPTIDE) AND OPTIMAL EXPERIMENT DESIGN
% =========================================================================
fprintf('Running Section 3: Two-Compartment Model (C-Peptide)...\n');

V1_c = 3.29; k01_c = 6.54e-2; k12_c = 5.68e-2; k21_c = 7.16e-2; d_c = 49650;

A3 = [-k01_c-k21_c, k12_c; k21_c, -k12_c];
B3 = [1; 0];
C3 = [1/V1_c, 0];
D3 = 0;

sys3 = ss(A3, B3, C3, D3);
[Y3, T3] = impulse(sys3, 0:0.5:180);

% Eigenvalues (System Constants)
eigenvalues3 = abs(eig(A3));
alpha = max(eigenvalues3);
beta = min(eigenvalues3);

% Macro-constants for Optimal Sampling Design
A_opt = d_c * (alpha - k12_c) / ((alpha - beta) * V1_c);
B_opt = d_c * (beta - k12_c) / ((beta - alpha) * V1_c);

% Optimal Sampling Times via fmincon
n_samples = 5;
t0 = linspace(1, 180, n_samples)';
lb = zeros(n_samples, 1);
ub = 180 * ones(n_samples, 1);

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
optimal_times = fmincon(@(t) optimal_sampling_obj(t, A_opt, B_opt, alpha, beta), ...
                        t0, [], [], [], [], lb, ub, [], options);

opt_conc = A_opt * exp(-alpha * optimal_times) + B_opt * exp(-beta * optimal_times);

figure('Name', 'Two-Compartment Model & Optimal Sampling Times');
plot(T3, d_c * Y3, 'b-', 'LineWidth', 1.5); hold on;
scatter(optimal_times, opt_conc, 60, 'r', 'filled'); grid on;
title('C-Peptide Response with Optimal Sampling Schedule');
xlabel('Time [min]'); ylabel('Concentration [pmol/L]');
legend('Continuous Response', 'Optimal Sampling Instants');

%% =========================================================================
% 4. THREE-COMPARTMENT MODEL (POLE-ZERO CANCELLATION ANALYSIS)
% =========================================================================
fprintf('Running Section 4: Three-Compartment Model...\n');

V1_3c = 5; k21_3c = 2.22; k12_3c = 0.859; 
k31_3c = 0.031; k13_3c = 0.008; k01_3c = 1.2; d_3c = 500;

A4 = [-k21_3c-k31_3c-k01_3c, k12_3c, k13_3c; 
       k21_3c,             -k12_3c, 0; 
       k31_3c,              0,     -k13_3c];
B4 = [1; 0; 0];
C4 = [1/V1_3c, 0, 0];
D4 = 0;

sys4 = ss(A4, B4, C4, D4);
G4 = tf(sys4);

% Residue Analysis to Detect Pole-Zero Cancellation
[r4, p4, ~] = residue(G4.Numerator{1} * d_3c, G4.Denominator{1});

figure('Name', 'Three-Compartment Impulse Response');
[Y4, T4] = impulse(sys4, 0:0.1:20);
plot(T4, d_3c * Y4, 'm-', 'LineWidth', 1.5); grid on;
title('Three-Compartment Model Impulse Response');
xlabel('Time [h]'); ylabel('Concentration [mg/L]');

%% =========================================================================
% 5. REPEATED DOSING REGIMEN (PULSE TRAIN SIMULATION)
% =========================================================================
fprintf('Running Section 5: Repeated Dosing Regimen...\n');

tau_vec = [12, 8, 6, 4]; % Dosing intervals [hours]
t_sim = 0:0.05:72;       % 72-hour simulation timeline

figure('Name', 'Repeated Dosing Accumulation');
for k = 1:length(tau_vec)
    tau_i = tau_vec(k);
    u = zeros(size(t_sim));
    
    % Generate pulse train for repeated bolus administration
    pulse_indices = 1:round(tau_i / 0.05):length(t_sim);
    u(pulse_indices) = d;
    
    Y_rep = lsim(sys1, u, t_sim);
    
    subplot(2, 2, k);
    plot(t_sim, Y_rep, 'LineWidth', 1.2); grid on;
    title(sprintf('Dosing Interval \\tau = %d hours', tau_i));
    xlabel('Time [h]'); ylabel('Concentration [mg/L]');
end

fprintf('All Pharmacokinetic Models Executed Successfully!\n');

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function J = optimal_sampling_obj(t, A, B, alpha, beta)
    % Computes the D-optimality criterion (determinant of Fisher Info Matrix)
    % for parameter estimation in a 2-exponential decay model.
    t = t(:);
    
    % Sensitivity matrix entries
    dy_dA = exp(-alpha * t);
    dy_dB = exp(-beta * t);
    dy_dalpha = -A * t .* exp(-alpha * t);
    dy_dbeta = -B * t .* exp(-beta * t);
    
    FIM = [dy_dA, dy_dB, dy_dalpha, dy_dbeta]' * [dy_dA, dy_dB, dy_dalpha, dy_dbeta];
    
    % Minimize negative log-determinant (D-optimal design)
    J = -log(det(FIM) + 1e-8);
end