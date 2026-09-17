function [r2adj, Ftest, bhat, yhat] = regressione(n, k, x, y)
    % REGRESSION Analysis using Ordinary Least Squares (OLS)
    % Computes regression coefficients, fitted values, and model quality metrics.
    %
    % Inputs:
    %   n - Number of observations
    %   k - Number of regressors
    %   x - n x k Matrix of regressors
    %   y - n x 1 Response vector (dependent variable)
    %
    % Outputs:
    %   r2adj - Adjusted R-squared of the model
    %   Ftest - F-test statistic for overall model significance
    %   bhat  - Estimated regression coefficients (parameters)
    %   yhat  - Fitted response values

    % Construct the design matrix H by adding a column of ones for the intercept
    H = [ones(n, 1) x]; % Dimensions: n x (k+1)
    
    % Estimate coefficients using OLS formula: b = (H' * SIGMA_V * H)^-1 * H' * SIGMA_V * y
    % Here, SIGMA_V is the identity matrix (eye(n)) assuming homoscedasticity 
    % (constant variance of residuals).
    bhat = inv(H' * eye(n) * H) * H' * eye(n) * y; 
    
    % Calculate fitted values
    yhat = H * bhat; 
    
    % Evaluate the model performance using the external adjusted R-squared function
    [r2adj, Ftest] = r2adjusted(n, k, y, yhat);
    
end