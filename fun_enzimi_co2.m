function dx = fun_enzimi_co2(t, x, kpiu_1, kmeno_1, kpiu_2)
    % State vector x = [s, e, c, p]: 
    % Corresponds to the initial conditions [s0, e0, c0, p0] passed to ode45.
    % s = substrate, e = enzyme, c = complex, p = product.
    
    dx = zeros(size(x));
    
    dx(1) = -kpiu_1 * x(1) * x(2) + kmeno_1 * x(3);       % ds/dt: Rate of change of substrate concentration
    dx(2) = kmeno_1 * x(3) - kpiu_1 * x(1) * x(2) + kpiu_2 * x(3); % de/dt: Rate of change of enzyme concentration
    dx(3) = kpiu_1 * x(1) * x(2) - kpiu_2 * x(3) - kmeno_1 * x(3); % dc/dt: Rate of change of the enzyme-substrate complex
    dx(4) = kpiu_2 * x(3);                               % dp/dt: Rate of change of product concentration
    
end