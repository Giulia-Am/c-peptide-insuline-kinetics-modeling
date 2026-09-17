function [GCV] = crossvalidation(G_regolarizzazione, SIGMA_V, P, y, gamma) 
% Given the chosen gamma, calculate the estimate and deconvolution
uhat = (G_regolarizzazione' * inv(SIGMA_V) * G_regolarizzazione + gamma * (P' * P))^-1 * (G_regolarizzazione' * inv(SIGMA_V) * y);
yhat = G_regolarizzazione * uhat;

% Calculate the residuals 
residuihat = y - yhat; % res = z - y
RSShat = sum(residuihat.^2);

% Influence matrix (Hat matrix) 
Hhat = G_regolarizzazione * inv((G_regolarizzazione' * inv(SIGMA_V) * G_regolarizzazione) + (gamma * (P' * P))) * G_regolarizzazione' * inv(SIGMA_V);
qhat = trace(Hhat);
n = length(y);
GCV = (n * residuihat' * inv(SIGMA_V) * residuihat) / (n - qhat).^2;
end