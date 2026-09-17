function dx = fun_enzimi_co2_anticompetitiva(t, x, kpiu_1, kmeno_1, kpiu_2, kpiu_3, kmeno_3)
    % State vector x = [s, e, i, c1, c2, p]:
    % Corresponds to the initial conditions [s0, e0, i0, c10, c20, p0] passed to ode45.
    % s  = substrate
    % e  = free enzyme
    % i  = inhibitor
    % c1 = enzyme-substrate complex (ES)
    % c2 = enzyme-substrate-inhibitor complex (ESI) - uncompetitive complex
    % p  = product
    
    dx = zeros(size(x));
    
    dx(1) = kmeno_1 * x(4) - kpiu_1 * x(1) * x(2);                               % ds/dt: Rate of change of substrate concentration
    dx(2) = kmeno_1 * x(4) - kpiu_1 * x(1) * x(2) + kpiu_2 * x(4);               % de/dt: Rate of change of free enzyme concentration
    dx(3) = kmeno_3 * x(5) - kpiu_3 * x(4) * x(3);                               % di/dt: Rate of change of inhibitor concentration
    dx(4) = kpiu_1 * x(1) * x(2) - kpiu_2 * x(4) - kmeno_1 * x(4) ...
            + kmeno_3 * x(5) - kpiu_3 * x(4) * x(3);                             % dc1/dt: Rate of change of the enzyme-substrate complex (ES)
    dx(5) = -kmeno_3 * x(5) + kpiu_3 * x(4) * x(3);                              % dc2/dt: Rate of change of the ESI complex
    dx(6) = kpiu_2 * x(4);                                                       % dp1/dt: Rate of change of product concentration
    
end