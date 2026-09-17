function dx = fun_enzimi_co2_noncompetitiva(t, x, kpiu_1, kmeno_1, kpiu_2, kpiu_3, kmeno_3, kpiu_4)
    % State vector x = [s, e, i, c1, c2, c3, p]:
    % Corresponds to the initial conditions [s0, e0, i0, c10, c20, c30, p0] passed to ode45.
    % s  = substrate
    % e  = free enzyme
    % i  = inhibitor
    % c1 = enzyme-substrate complex (ES)
    % c2 = enzyme-inhibitor complex (EI) or intermediate complex
    % c3 = ternary/other complex (e.g., EIS)
    % p  = product
    
    dx = zeros(size(x));
    
    dx(1) = kmeno_1 * x(4) - kpiu_1 * x(1) * x(2) - kpiu_4 * x(5) * x(1);       % ds/dt: Rate of change of substrate concentration
    dx(2) = kmeno_1 * x(4) - kpiu_1 * x(1) * x(2) + kpiu_2 * x(4) + kmeno_3 * x(1) * x(5) - kpiu_3 * x(2) * x(3); % de/dt: Rate of change of enzyme concentration
    dx(3) = kmeno_3 * x(5) * x(1) - kpiu_3 * x(2) * x(3);                       % di/dt: Rate of change of inhibitor concentration
    dx(4) = kpiu_1 * x(1) * x(2) - kpiu_2 * x(4) - kmeno_1 * x(4);               % dc1/dt: Rate of change of the first complex (ES)
    dx(5) = -kmeno_3 * x(5) * x(1) + kpiu_3 * x(2) * x(3) - kpiu_4 * x(1) * x(5); % dc2/dt: Rate of change of the second complex (EI)
    dx(6) = kpiu_4 * x(1) * x(5);                                               % dc3/dt: Rate of change of the third complex (EIS)
    dx(7) = kpiu_2 * x(4);                                                      % dp/dt: Rate of change of product concentration
    
end