%% ENZYMATIC REACTIONS SIMULATION: CINETICS AND INHIBITION ANALYSIS
% This script evaluates various enzymatic reaction scenarios (enzyme/substrate 
% ratios, carbonic anhydrase kinetics without recombination, and different types 
% of enzyme inhibition: competitive, uncompetitive, and non-competitive).

clear all
close all
clc

%% =========================================================================
% 1. CASE 1: ENZYME CONCENTRATION MUCH GREATER THAN SUBSTRATE
% =========================================================================
% - The substrate is completely consumed, whereas the enzyme is not fully consumed.

kpiu_1 = 1;     % mM^-1 sec^-1
kmeno_1 = 10;   % sec^-1
kpiu_2 = 1;     % sec^-1
kmeno_2 = 10;   % mM^-1 sec^-1
s0 = 1;         % mM
e0 = 100;       % mM
c0 = 0;         % mM
p0 = 0;         % mM

% Simulate the time course of substrate, enzyme, complex, and product concentrations over 0.1 seconds
[T, X] = ode23(@(t,x) fun_enzimi(t, x, kpiu_1, kmeno_1, kpiu_2, kmeno_2), [0 10e-2], [s0 e0 c0 p0]);
% x = [s, e, c, p]

figure(1);
subplot(3, 2, 1);
plot(T, X(:, 1));
title('Substrate Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 2);
plot(T, X(:, 2));
title('Enzyme Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 3);
plot(T, X(:, 3));
title('Complex Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 4);
plot(T, X(:, 4));
title('Product Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 5);
semilogy(T, X(:, 4));
title('Product Concentration vs Time (Semilog Scale)');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 6);
semilogy(T, X(:, 1));
title('Substrate Concentration vs Time (Semilog Scale)');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

% Free and bound enzyme at steady-state:
enzima_libero_regime = X(end, 2) / max(X(:, 2)) * 100; % Free enzyme at steady-state: ~99%
enzima_legato_regime = 100 - enzima_libero_regime;      % Bound enzyme at steady-state: ~1%

% Product formation rate
figure(2);
plot(T, gradient(X(:, 4), T));
title('Product Formation Rate over Time (dp/dt)');
xlabel('Time [sec]');
ylabel('Velocity [mM*sec^-1]');

% QUESTION 1:
% The concentration profile of P is neither exponential nor linear (observed from the log scale).
% The concentration profile of S is exponential (observed from the log scale) because 
% of the high enzyme concentration, causing S to decrease rapidly.
% QUESTION 2:
% The steady-state value is reached when the product formation rate vanishes, i.e., dp/dt = 0.
% From the plot, the product concentration reaches its maximum in approximately 0.1 seconds, 
% and the velocity drops to zero at the exact same point.
% QUESTION 3:
% At steady-state, the free enzyme is 99% of the total, and the bound enzyme is 1%, 
% given that s0 << e0. Furthermore, not all enzyme is free because k2meno is non-zero, 
% allowing the free enzyme to rebind to the product.


%% =========================================================================
% 2. CASE 2: ENZYME CONCENTRATION MUCH LESS THAN SUBSTRATE
% =========================================================================
% - Substrate is not fully consumed, while the enzyme is consumed instantaneously.

kpiu_1 = 1;     % mM^-1 sec^-1
kmeno_1 = 1;    % sec^-1
kpiu_2 = 1;     % sec^-1
kmeno_2 = 1;    % mM^-1 sec^-1
s0 = 100;       % mM
e0 = 1;         % mM
c0 = 0;         % mM
p0 = 0;         % mM

% RAPID KINETICS: Simulate over 0.15 seconds
[Tr, Xr] = ode23(@(t,x) fun_enzimi(t, x, kpiu_1, kmeno_1, kpiu_2, kmeno_2), [0 0.15], [s0 e0 c0 p0]);

figure('Name', 'Rapid Kinetics');
subplot(3, 2, 1);
plot(Tr, Xr(:, 1));
title('Substrate Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 2);
plot(Tr, Xr(:, 2));
title('Enzyme Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 3);
plot(Tr, Xr(:, 3));
title('Complex Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 4);
plot(Tr, Xr(:, 4));
title('Product Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 5);
semilogy(Tr, Xr(:, 4));
title('Product Concentration vs Time (Semilog Scale)');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(3, 2, 6);
semilogy(Tr, Xr(:, 1));
title('Substrate Concentration vs Time (Semilog Scale)');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

% SLOW KINETICS: Simulate over 300 seconds
[Tl, Xl] = ode23(@(t,x) fun_enzimi(t, x, kpiu_1, kmeno_1, kpiu_2, kmeno_2), [0 300], [s0 e0 c0 p0]);

figure('Name', 'Slow Kinetics');
subplot(2, 2, 1);
plot(Tl, Xl(:, 1));
title('Substrate Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(2, 2, 2);
plot(Tl, Xl(:, 2));
title('Enzyme Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(2, 2, 3);
plot(Tl(1:80), Xl(1:80, 3));
title('Complex Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

subplot(2, 2, 4);
plot(Tl, Xl(:, 4));
title('Product Concentration vs Time');
xlabel('Time [sec]');
ylabel('Concentration [mM]');

% Free and bound enzyme at steady-state:
enzima_libero_regime = Xl(end, 2) / max(Xl(:, 2)) * 100; % Free enzyme at steady-state: ~99%
enzima_legato_regime = 100 - enzima_libero_regime;      % Bound enzyme at steady-state: ~1%

% Product formation rate
figure('Name', 'Product Formation Rate');
plot(Tl, gradient(Xl(:, 4), Tl));
title('Product Formation Rate over Time (dp/dt)');
xlabel('Time [sec]');
ylabel('Velocity [mM*sec^-1]');


%% =========================================================================
% 3. CASE 3: CARBONIC ANHYDRASE (NO RECOMBINATION: k-2 = 0)
% =========================================================================

kpiu_1 = 75e3;  % mM^-1 sec^-1
kmeno_1 = 75;   % sec^-1
kpiu_2 = 600e3; % sec^-1
s0 = 100;       % mM
e0 = 1;         % mM
c0 = 0;         % mM
p0 = 0;         % mM

% Simulate concentration profiles over 100 nanoseconds (Rapid Kinetics)
[T, X] = ode45(@(t,x) fun_enzimi_co2(t, x, kpiu_1, kmeno_1, kpiu_2), [0 1e-8], [s0 e0 c0 p0]);

figure('Name', 'Rapid Kinetics (No Recombination)');
subplot(3, 2, 1); plot(T, X(:, 1)); title('[S] Concentration vs Time'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 2); plot(T, X(:, 2)); title('[E] Concentration vs Time'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 3); semilogy(T, X(:, 4)); title('[P] Concentration vs Time (Semilog Scale)'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 4); plot(T, X(:, 3)); title('[ES] Concentration vs Time'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 5); plot(T, X(:, 4)); title('[P] Concentration vs Time'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 6); plot(T, gradient(X(:, 4), T)); title('Product Formation Rate over Time (dp/dt)'); xlabel('Time [sec]'); ylabel('Velocity [mM*sec^-1]');

% Simulate concentration profiles over 10 milliseconds (Slow Kinetics)
[T, X] = ode45(@(t,x) fun_enzimi_co2(t, x, kpiu_1, kmeno_1, kpiu_2), [0 1e-3], [s0 e0 c0 p0]);

figure('Name', 'Slow Kinetics (No Recombination)');
subplot(3, 2, 1); plot(T, X(:, 1)); title('[S] Concentration vs Time'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 2); plot(T, X(:, 2)); title('[E] Concentration vs Time'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 3); semilogy(T, X(:, 4)); title('[P] Concentration vs Time (Semilog Scale)'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 4); plot(T, X(:, 3)); title('[ES] Concentration vs Time'); xlabel('Time [sec]'); ylabel('Concentration [mM]');
subplot(3, 2, 5); plot(T, X(:, 4)); title('[P] Concentration vs Time'); xlabel('Time [sec]'); ylabel('Concentration [mM]');

velocita_prodotto = kpiu_2 * X(:, 3);
subplot(3, 2, 6); plot(T, velocita_prodotto); title('Product Formation Rate over Time (dp/dt)'); xlabel('Time [sec]'); ylabel('Velocity [mM*sec^-1]');

% Free and bound enzyme at steady-state:
enzima_libero_regime = X(end, 2) / max(X(:, 2)) * 100; 
enzima_legato_regime = 100 - enzima_libero_regime;      

% Evaluations on tau0:
tau0 = (kpiu_1 * s0 + kmeno_1 + kpiu_2)^-1;
esaurimento_transitorio = 5 * tau0;
velocita_complesso = kpiu_1 .* X(:, 2) .* X(:, 1) - kpiu_2 .* X(:, 3) - kmeno_1 .* X(:, 3);

figure('Name', 'Transient Analysis (tau0)');
plot(T(1:100), velocita_complesso(1:100)); 
title('Complex Formation Velocity: Transient Analysis (dc/dt)');
xlabel('Time [sec]');
ylabel('Velocity [mM*sec^-1]');
xline(esaurimento_transitorio, '--r', '5*tau0');

% MICHAELIS-MENTEN LAW AS A FUNCTION OF s0:
s0_mm = linspace(0, 100, 200);
km = (kpiu_2 + kmeno_1) / kpiu_1;
Vmax = kpiu_2 * e0;

for i = 1:length(s0_mm)
    [TMM, XMM] = ode45(@(t,x) fun_enzimi_co2(t, x, kpiu_1, kmeno_1, kpiu_2), [0 5e-6], [s0_mm(i) e0 c0 p0]);
    v(i) = (Vmax * s0_mm(i)) / (km + s0_mm(i));
end

figure('Name', 'Michaelis-Menten Kinetics vs s0');
plot(s0_mm, v);
hold on;
plot(X(:, 1), velocita_prodotto);
title('Michaelis-Menten Velocity vs s0');
xlabel('s0 -- s(t)');
ylabel('Velocity [mM*sec^-1]');
yline(Vmax, '--r', 'Vmax');
xline(km, '--r', 'Km');
yline(Vmax / 2, '--r', 'Vmax/2');
legend('MM', 'Vp');


%% =========================================================================
% 4. CASE 4: CARBONIC ANHYDRASE WITH INHIBITION MECHANISMS
% =========================================================================

% --- COMPETITIVE INHIBITION ---
kpiu_1 = 75e3;  % mM^-1 sec^-1
kmeno_1 = 75;   % sec^-1
kpiu_2 = 600e3; % sec^-1
kpiu_3 = 75e2;  % mM^-1 sec^-1
kmeno_3 = 75;   % sec^-1
kpiu_4 = 600e3; % sec^-1
s0 = 100;       % mM
e0 = 1;         % mM
c10 = 0;        % mM
c20 = 0;        % mM
p10 = 0;        % mM
p20 = 0;        % mM
i0 = [0 10 20 50 100];

figure('Name', 'Competitive Inhibition vs i0');
for i = 1:length(i0)
    [T, X] = ode45(@(t,x) fun_enzimi_co2_competitiva(t, x, kpiu_1, kmeno_1, kpiu_2, kpiu_3, kmeno_3), [0 10e-6], [s0 e0 i0(i) c10 c20 p10]);
    
    velocita_prodotto = kpiu_2 * X(:, 4);
    Velocita_prodotto{i} = velocita_prodotto;
    enzima_libero_regime(i) = X(end, 2) / max(X(:, 2)) * 100;
    enzima_legato_regime(i) = 100 - enzima_libero_regime(i);
    
    subplot(3, 3, 1); plot(T, velocita_prodotto); hold on;
    subplot(3, 3, 2); plot(T, X(:, 5)); hold on;
    subplot(3, 3, 3); plot(T, X(:, 4)); hold on;
    subplot(3, 3, 4); plot(T, X(:, 1)); hold on;
    subplot(3, 3, 5); plot(T, X(:, 2)); hold on;
    subplot(3, 3, 6); plot(T, X(:, 3)); hold on;
    subplot(3, 3, 7); plot(T, X(:, 6)); hold on;
end

cell_legends = arrayfun(@(val) sprintf('i0=%d', val), i0, 'UniformOutput', false);
subtitles = {'Vp Product Formation vs i0', '[EI] vs i0', '[ES] vs i0', '[S] vs i0', '[E] vs i0', '[I] vs i0', '[P] vs i0'};
for k = 1:7
    subplot(3, 3, k);
    title(subtitles{k});
    xlabel('Time [sec]');
    ylabel('Concentration/Velocity');
    legend(cell_legends, 'Location', 'best');
end

% --- COMPETITIVE INHIBITION vs k3 ---
kpiu_3_vals = [2500 5000 7500 10000];
i0_fixed = 10; 
figure('Name', 'Competitive Inhibition vs k3');
for i = 1:length(kpiu_3_vals)
    [Tk, Xk] = ode45(@(t,x) fun_enzimi_co2_competitiva(t, x, kpiu_1, kmeno_1, kpiu_2, kpiu_3_vals(i), kmeno_3), [0:1e-9:1e-5], [s0 e0 i0_fixed c10 c20 p10]);
    
    velocita_prodottok = gradient(Xk(:, 6), Tk);
    
    subplot(3, 3, 1); plot(Tk(1:100), velocita_prodottok(1:100)); hold on;
    subplot(3, 3, 2); plot(Tk, Xk(:, 5)); hold on;
    subplot(3, 3, 3); plot(Tk, Xk(:, 4)); hold on;
    subplot(3, 3, 4); plot(Tk, Xk(:, 1)); hold on;
    subplot(3, 3, 5); plot(Tk, Xk(:, 2)); hold on;
    subplot(3, 3, 6); plot(Tk(1:50), Xk(1:50, 3)); hold on;
    subplot(3, 3, 7); plot(Tk, Xk(:, 6)); hold on;
end

cell_legends_k3 = arrayfun(@(val) sprintf('k3=%d', val), kpiu_3_vals, 'UniformOutput', false);
subtitles_k3 = {'Vp Product Formation vs k3', '[EI] vs k3', '[ES] vs k3', '[S] vs k3', '[E] vs k3', '[I] vs k3', '[P] vs k3'};
for k = 1:7
    subplot(3, 3, k);
    title(subtitles_k3{k});
    xlabel('Time [sec]');
    ylabel('Concentration/Velocity');
    legend(cell_legends_k3, 'Location', 'best');
end


% --- UNCOMPETITIVE INHIBITION ---
kpiu_3 = 75e2;  
figure('Name', 'Uncompetitive Inhibition vs i0');
for i = 1:length(i0)
    [Tant, Xant] = ode45(@(t,x) fun_enzimi_co2_anticompetitiva(t, x, kpiu_1, kmeno_1, kpiu_2, kpiu_3, kmeno_3), [0 1e-5], [s0 e0 i0(i) c10 c20 p10]);
    
    velocita_prodottoant = kpiu_2 * Xant(:, 4);
    
    subplot(3, 3, 1); plot(Tant, velocita_prodottoant); hold on;
    subplot(3, 3, 2); plot(Tant, Xant(:, 5)); hold on;
    subplot(3, 3, 3); plot(Tant, Xant(:, 4)); hold on;
    subplot(3, 3, 4); plot(Tant, Xant(:, 1)); hold on;
    subplot(3, 3, 5); plot(Tant, Xant(:, 2)); hold on;
    subplot(3, 3, 6); plot(Tant(1:50), Xant(1:50, 3)); hold on;
    subplot(3, 3, 7); plot(Tant, Xant(:, 6)); hold on;
end

for k = 1:7
    subplot(3, 3, k);
    title(subtitles{k});
    xlabel('Time [sec]');
    ylabel('Concentration/Velocity');
    legend(cell_legends, 'Location', 'best');
end

% --- UNCOMPETITIVE INHIBITION vs k3 ---
figure('Name', 'Uncompetitive Inhibition vs k3');
for i = 1:length(kpiu_3_vals)
    [Tkant, Xkant] = ode45(@(t,x) fun_enzimi_co2_anticompetitiva(t, x, kpiu_1, kmeno_1, kpiu_2, kpiu_3_vals(i), kmeno_3), [0 1e-5], [s0 e0 i0_fixed c10 c20 p10]);
    
    velocita_prodottokant = kpiu_2 * Xkant(:, 4);
    
    subplot(3, 3, 1); plot(Tkant, velocita_prodottokant); hold on;
    subplot(3, 3, 2); plot(Tkant, Xkant(:, 5)); hold on;
    subplot(3, 3, 3); plot(Tkant, Xkant(:, 4)); hold on;
    subplot(3, 3, 4); plot(Tkant, Xkant(:, 1)); hold on;
    subplot(3, 3, 5); plot(Tkant, Xkant(:, 2)); hold on;
    subplot(3, 3, 6); plot(Tkant(1:50), Xkant(1:50, 3)); hold on;
    subplot(3, 3, 7); plot(Tkant, Xkant(:, 6)); hold on;
end

for k = 1:7
    subplot(3, 3, k);
    title(subtitles_k3{k});
    xlabel('Time [sec]');
    ylabel('Concentration/Velocity');
    legend(cell_legends_k3, 'Location', 'best');
end


% --- NON-COMPETITIVE INHIBITION ---
c30 = 0;        
figure('Name', 'Non-Competitive Inhibition vs i0');
for i = 1:length(i0)
    [Tnon, Xnon] = ode45(@(t,x) fun_enzimi_co2_noncompetitiva(t, x, kpiu_1, kmeno_1, kpiu_2, kpiu_3, kmeno_3, kpiu_4), [0 1e-5], [s0 e0 i0(i) c10 c20 c30 p10]);
    
    velocita_prodottonon = kpiu_2 * Xnon(:, 4);
    
    subplot(3, 3, 1); plot(Tnon, velocita_prodottonon); hold on;
    subplot(3, 3, 2); plot(Tnon, Xnon(:, 5)); hold on;
    subplot(3, 3, 3); plot(Tnon, Xnon(:, 6)); hold on;
    subplot(3, 3, 4); plot(Tnon, Xnon(:, 4)); hold on;
    subplot(3, 3, 5); plot(Tnon, Xnon(:, 1)); hold on;
    subplot(3, 3, 6); plot(Tnon, Xnon(:, 2)); hold on;
    subplot(3, 3, 7); plot(Tnon(1:50), Xnon(1:50, 3)); hold on;
    subplot(3, 3, 8); plot(Tnon, Xnon(:, 7)); hold on;
end

subtitles_non = {'Vp Product Formation vs i0', '[EI] vs i0', '[EIS] vs i0', '[ES] vs i0', '[S] vs i0', '[E] vs i0', '[I] vs i0', '[P] vs i0'};
for k = 1:8
    subplot(3, 3, k);
    title(subtitles_non{k});
    xlabel('Time [sec]');
    ylabel('Concentration/Velocity');
    legend(cell_legends, 'Location', 'best');
end

% --- NON-COMPETITIVE INHIBITION vs k3 ---
figure('Name', 'Non-Competitive Inhibition vs k3');
for i = 1:length(kpiu_3_vals)
    [Tknon, Xknon] = ode45(@(t,x) fun_enzimi_co2_noncompetitiva(t, x, kpiu_1, kmeno_1, kpiu_2, kpiu_3_vals(i), kmeno_3, kpiu_4), [0 1e-5], [s0 e0 i0_fixed c10 c20 c30 p10]);
    
    velocita_prodottoknon = kpiu_2 * Xknon(:, 4);
    
    subplot(3, 3, 1); plot(Tknon, velocita_prodottoknon); hold on;
    subplot(3, 3, 2); plot(Tknon, Xknon(:, 5)); hold on;
    subplot(3, 3, 3); plot(Tknon, Xknon(:, 6)); hold on;
    subplot(3, 3, 4); plot(Tknon, Xknon(:, 4)); hold on;
    subplot(3, 3, 5); plot(Tknon, Xknon(:, 1)); hold on;
    subplot(3, 3, 6); plot(Tknon, Xknon(:, 2)); hold on;
    subplot(3, 3, 7); plot(Tknon(1:50), Xknon(1:50, 3)); hold on;
    subplot(3, 3, 8); plot(Tknon, Xknon(:, 7)); hold on;
end

subtitles_knon = {'Vp Product Formation vs k3', '[EI] vs k3', '[EIS] vs k3', '[ES] vs k3', '[S] vs k3', '[E] vs k3', '[I] vs k3', '[P] vs k3'};
for k = 1:8
    subplot(3, 3, k);
    title(subtitles_knon{k});
    xlabel('Time [sec]');
    ylabel('Concentration/Velocity');
    legend(cell_legends_k3, 'Location', 'best');
end