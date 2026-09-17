function res = residui(p, conc, t, w)
% MONO-EXPONENTIAL MODEL
if length(p) == 2
    y = p(1) * exp(-p(2) .* t); % y = A*exp(-alpha*t)
end

% BI-EXPONENTIAL MODEL
if length(p) == 4
    y = (p(1) * exp(-p(2) .* t)) + (p(3) * exp(-p(4) .* t)); % y = A*exp(-alpha*t) + B*exp(-beta*t)
end

% TRI-EXPONENTIAL MODEL
if length(p) == 6
    y = (p(1) * exp(-p(2) .* t)) + (p(3) * exp(-p(4) .* t)) + (p(5) * exp(-p(6) .* t)); % y = A*exp(-alpha*t) + B*exp(-beta*t) + C*exp(-gamma*t)
end

if length(w) > 1
    res = (conc - y) ./ w; % WLS CASE: WEIGHTED_RESIDUALS = (z-y)/w = (z-h(t,phat))/w
    % Division by w is performed because w contains the standard deviation of V. 
    % By definition, we should multiply the residuals by the matrix. 
    % If W = inverse of variance -> res*W coincides with doing res./w 
    % (where w = std(V) and not std(V)^-1)
else
    res = conc - y; % LS CASE: RESIDUALS = (z-y) = (z-h(t,phat)) NOTE: WE DO NOT COMPUTE THE SQUARE BECAUSE lsqnonlin DOES IT INTERNALLY!
end
end