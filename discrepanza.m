function [discrep] = discrepanza(G_regolarizzazione, SIGMA_V, P, y, gamma)
% Calculate the estimate for the chosen gamma and deconvolution
uhat = (G_regolarizzazione' * inv(SIGMA_V) * G_regolarizzazione + gamma * (P' * P))^-1 * (G_regolarizzazione' * inv(SIGMA_V) * y);
yhat = G_regolarizzazione * uhat;

% Calculate the residual vector 
residuihat = y - yhat; % res = z - y

% Calculate the squared norm of the residuals 
RSShat = sum(residuihat.^2);

% Repeat until RSS is approximately the trace of the covariance matrix 
discrep = abs(RSShat - trace(SIGMA_V));
end