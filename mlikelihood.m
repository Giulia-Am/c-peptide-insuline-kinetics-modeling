function [ML] = mlikelihood(G_regolarizzazione, SIGMA_V, P, y, gamma)
uhat = (G_regolarizzazione' * inv(SIGMA_V) * G_regolarizzazione + gamma * (P' * P))^-1 * (G_regolarizzazione' * inv(SIGMA_V) * y);
yhat = G_regolarizzazione * uhat;
residuihat = y - yhat; % res = z - y
RSShat = sum(residuihat.^2);
Hhat = G_regolarizzazione * inv((G_regolarizzazione' * inv(SIGMA_V) * G_regolarizzazione) + (gamma * (P' * P))) * G_regolarizzazione' * inv(SIGMA_V);
qhat = trace(Hhat);
n = length(y);
ML = abs((residuihat' * inv(SIGMA_V) * residuihat * qhat) / ((n - qhat) * uhat' * P' * P * uhat));
end