%% NONLINEAR COMPARTMENTAL MODELS ANALYSIS (SINGLE AND REPEATED DOSING)
% This script simulates and compares linear vs. nonlinear (Michaelis-Menten 
% and Hill kinetics) compartmental pharmacokinetic models with absorption.

clc; clear; close all;

%% =========================================================================
% 1. SINGLE-DOSE NONLINEAR MODEL WITH MICHAELIS-MENTEN ABSORPTION
% =========================================================================
fprintf('Running Section 1: Single-Dose Michaelis-Menten Absorption...\n');

V2 = 5;          % Volume [L]
k01 = 1.2;       % Transfer rate constant [1/h]
k02 = 1.2;       % Elimination rate constant [1/h]
Vmax = 110;      % Maximum velocity [mg/h]
km = 50;         % Michaelis constant [mg]
dose0 = 500;     % Initial dose [mg]

t0 = (0:0.1:6)'; % Time vector for elimination (~5/k01)
q0 = [dose0; 0]; % Initial quantities in compartments 1 and 2

% Solve MM Nonlinear System via ode45
[t_MM, q_MM] = ode45(@(t, q) fun_MM(t, q, k01, k02, Vmax, km), t0, q0);

k21_MM = Vmax ./ (km + q_MM(:, 1));          % Nonlinear absorption rate [1/h]
F21_MM = k21_MM .* q_MM(:, 1);               % Input flux into compartment 2
c2_MM = q_MM(:, 2) / V2;                     % Concentration in accessible compartment 2
ka_MM = k01 + k21_MM;                        % Apparent absorption constant
F_fraction_MM = k21_MM ./ ka_MM;             % Bioavailability fraction

% Plot MM Kinetics Results
figure('Name', 'Michaelis-Menten Absorption Kinetics');
subplot(2,2,1);
plot(t_MM, q_MM(:, 1), 'c-', 'LineWidth', 1.5); grid on; hold on;
plot(t_MM, q_MM(:, 2), 'b-', 'LineWidth', 1.5);
xlabel('Time [h]'); ylabel('Quantity [mg]');
legend('q_1(t)', 'q_2(t)'); title('Drug Quantities in Compartments');

subplot(2,2,2);
plot(t_MM, c2_MM, 'b-', 'LineWidth', 1.5); grid on;
xlabel('Time [h]'); ylabel('Concentration [mg/L]');
title('Concentration in Compartment 2');

subplot(2,2,3);
plot(t_MM, k21_MM, 'b-', 'LineWidth', 1.5); grid on;
xlabel('Time [h]'); ylabel('k_{21} [1/h]');
title('Absorption Velocity (Saturation Curve)');

subplot(2,2,4);
plot(q_MM(:, 1), F21_MM, 'b-', 'LineWidth', 1.5); grid on;
xlabel('q_1 [mg]'); ylabel('F_{21}(t) [mg/h]');
title('Input Flux vs. Quantity q_1');

%% =========================================================================
% 2. COMPARISON: NONLINEAR (MM) VS. LINEAR SYSTEM
% =========================================================================
fprintf('Running Section 2: MM vs. Linear Model Comparison...\n');

k21_lin = 2.2;   % Linear absorption rate [1/h]
A_lin = [-(k01 + k21_lin), 0; k21_lin, -k02];
B_lin = [1; 0];
C_lin = [0, 1/V2];
D_lin = 0;

sys_lin = ss(A_lin, B_lin, C_lin, D_lin);
[c2_lin, t_lin, q_lin] = impulse(sys_lin * dose0, t0);
F21_lin = k21_lin * q_lin(:, 1);

figure('Name', 'MM vs Linear Kinetics Comparison');
subplot(2,2,1);
plot(t_MM, q_MM(:, 1), 'c-', t_MM, q_MM(:, 2), 'b-', ...
     t_lin, q_lin(:, 1), 'm--', t_lin, q_lin(:, 2), 'r--', 'LineWidth', 1.5); 
grid on; xlabel('Time [h]'); ylabel('Quantity [mg]');
legend('q_{MM1}', 'q_{MM2}', 'q_{lin1}', 'q_{lin2}');
title('Quantity Comparison');

subplot(2,2,2);
plot(t_MM, c2_MM, 'b-', t_MM, c2_lin, 'r--', 'LineWidth', 1.5); grid on;
xlabel('Time [h]'); ylabel('Concentration [mg/L]');
legend('c_{MM2}', 'c_{lin2}'); title('Concentration Comparison');

subplot(2,2,3);
plot(t_MM, k21_MM, 'b-', 'LineWidth', 1.5); hold on;
yline(k21_lin, 'r--', 'LineWidth', 1.5); grid on;
xlabel('Time [h]'); ylabel('k_{21} [1/h]');
legend('k_{21,MM}', 'k_{21,lin}'); title('Absorption Rate Constant');

%% =========================================================================
% 3. PARAMETER VARIATION: 2*Vmax AND km/2 SENSITIVITY ANALYSIS
% =========================================================================
fprintf('Running Section 3: Parameter Sensitivity Analysis...\n');

Vmax_mod = Vmax * 2;
km_mod = km / 2;

[t_MM1, q_MM1] = ode45(@(t, q) fun_MM(t, q, k01, k02, Vmax_mod, km), t0, q0);
[t_MM2, q_MM2] = ode45(@(t, q) fun_MM(t, q, k01, k02, Vmax, km_mod), t0, q0);

c2_MM1 = q_MM1(:, 2) / V2;
c2_MM2 = q_MM2(:, 2) / V2;

figure('Name', 'Parameter Sensitivity: 2*Vmax & km/2');
subplot(2,1,1);
plot(t_MM, c2_MM, 'b-', t_MM, c2_MM1, 'r-', t_MM, c2_MM2, 'k-', 'LineWidth', 1.5);
grid on; xlabel('Time [h]'); ylabel('Concentration [mg/L]');
legend('Standard MM', '2*Vmax', 'km/2');
title('Concentration Response to Parameter Modifications');

subplot(2,1,2);
semilogy(t_MM, c2_MM, 'b-', t_MM, c2_MM1, 'r-', t_MM, c2_MM2, 'k-', 'LineWidth', 1.5);
grid on; xlabel('Time [h]'); ylabel('log(c) [mg/L]');
title('Semilogarithmic Concentration Response');

%% =========================================================================
% 4. REPEATED DOSING REGIMEN: LINEAR VS. NONLINEAR
% =========================================================================
fprintf('Running Section 4: Repeated Dosing Simulation...\n');

dose0 = 500;            % Dose [mg] (allineata con la Sezione 1)
tau = 4;                % Dosing interval [h]
n_pulses = round(72 / tau);

% Linear Repeated Dosing
X1_lin = []; X2_lin = []; T_lin = [];
[yf, tf, xf] = impulse(sys_lin, tau);
X1_lin = dose0 * xf(:, 1); 
X2_lin = dose0 * xf(:, 2); 
T_lin = tf;

for i = 1:n_pulses
    [yl, tl, xl] = initial(sys_lin, [X1_lin(end) + dose0, X2_lin(end)], tau);
    X1_lin = [X1_lin; xl(:, 1)];
    X2_lin = [X2_lin; xl(:, 2)];
    T_lin = [T_lin; tl + tau * i];
end

% Nonlinear Repeated Dosing
q0_rep = [dose0; 0];
X1_nl = []; X2_nl = []; T_nl = [];
[tf_nl, xf_nl] = ode45(@(t, q) fun_MM(t, q, k01, k02, Vmax, km), [0 tau], q0_rep);
X1_nl = xf_nl(:, 1); X2_nl = xf_nl(:, 2); T_nl = tf_nl;

for i = 1:n_pulses
    [tl_nl, xl_nl] = ode45(@(t, q) fun_MM(t, q, k01, k02, Vmax, km), [0 tau], [X1_nl(end) + dose0, X2_nl(end)]);
    X1_nl = [X1_nl; xl_nl(:, 1)];
    X2_nl = [X2_nl; xl_nl(:, 2)];
    T_nl = [T_nl; tl_nl + tau * i];
end

figure('Name', 'Repeated Dosing: Linear vs Nonlinear');
subplot(2,1,1);
plot(T_lin, X2_lin / V2, 'r-', 'LineWidth', 1.2); grid on;
title('Linear Repeated Dosing Concentration');
xlabel('Time [h]'); ylabel('Concentration [mg/L]');

subplot(2,1,2);
plot(T_nl, X2_nl / V2, 'b-', 'LineWidth', 1.2); grid on;
title('Nonlinear Repeated Dosing Concentration (MM Kinetics)');
xlabel('Time [h]'); ylabel('Concentration [mg/L]');

fprintf('All Nonlinear Simulations Executed Successfully!\n');

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function dq = fun_MM(~, q, k01, k02, Vmax, km)
    % Michaelis-Menten absorption kinetics differential equations
    q1 = q(1);
    q2 = q(2);
    
    k21 = Vmax / (km + q1);
    
    dq1 = -(k01 + k21) * q1;
    dq2 = k21 * q1 - k02 * q2;
    
    dq = [dq1; dq2];
end
