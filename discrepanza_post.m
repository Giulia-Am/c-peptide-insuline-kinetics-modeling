function [uhat, RSShat, residuihat] = discrepanza_post(G_regolarizzazione, SIGMA_V, P, y, gamma)
% Use this function to return uhat once the optimal gamma is found

uhat = (G_regolarizzazione' * inv(SIGMA_V) * G_regolarizzazione + gamma * (P' * P))^-1 * (G_regolarizzazione' * inv(SIGMA_V) * y);
yhat = G_regolarizzazione * uhat;
residuihat = y - yhat; % res = z - y
RSShat = sum(residuihat.^2);
end