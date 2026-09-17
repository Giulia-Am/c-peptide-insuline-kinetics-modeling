function [r2adj, Ftest] = r2adjusted(n, k, y, yhat)
    % R2ADJUSTED
    % Calculates the adjusted R-squared and the F-test statistic for a regression model.
    %
    % Inputs:
    %   n    - Number of observations
    %   k    - Number of regressors (parameters)
    %   y    - Actual response vector (dependent variable)
    %   yhat - Fitted response values from the model
    %
    % Outputs:
    %   r2adj - Adjusted R-squared value
    %   Ftest - F-test statistic for overall model significance

    SSR = sum((y - mean(y)).^2) - sum((y - yhat).^2); % Regression Sum of Squares (Explained variance)
    SSE = sum((y - yhat).^2);                         % Error Sum of Squares (Residual variance)
    SST = SSR + SSE;                                  % Total Sum of Squares (Total variance)
    
    r2 = SSR / SST;                                   % Standard R-squared
    
    % Compute adjusted R-squared based on observations (n) and regressors (k)
    r2adj = r2 - ((k / (n - k - 1)) * (1 - r2));
    
    % Evaluate overall model significance using the F-statistic (MSR / MSE)
    Ftest = (SSR / k) / (SSE / (n - k - 1)); 
    
end