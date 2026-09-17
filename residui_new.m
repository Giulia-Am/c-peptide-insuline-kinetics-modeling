function [res, y] = residui_new(p, conc, t, w)
d = 49650; % Bolus of 49650 pmol defined in Exercise 1
A = [-p(1)-p(3), p(2); p(3), -p(2)];
B = [1; 0];
C = [1/p(4), 0];
D = 0;

% Linear Dynamical System
SDL = ss(A, B, C, D);

% Find the predicted output as the impulse response of this SDL system:
[y, T] = impulse(d * SDL, 0:0.1:t(length(t)));

% PROBLEM: We would like to compute the impulse response on the time vector 
% provided by the data (t), but it does not have a constant step size.
% Therefore, compute the impulse response for all times from 0:0.1 up to 
% the final time provided by the data:
for i = 1:length(t)
    indici(i) = find(T == t(i)); % Save only the times that coincide with those provided by the data (t)
end

y = y(indici); % Consequently, save only the predictions corresponding to the times provided by the dataset (t)

if length(w) > 1
    res = (conc - y) ./ w; % Residuals would be (conc-y)^2, but we don't square them since lsqnonlin does it internally
else
    res = conc - y; % Residuals would be (conc-y)^2, but we don't square them since lsqnonlin does it internally
end
end
