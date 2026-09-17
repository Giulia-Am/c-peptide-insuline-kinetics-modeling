function dx = fun_enzimi(t, x, kpiu_1, kmeno_1, kpiu_2, kmeno_2)
    % State vector x = [s, e, c, p]:
    % Corresponds to the initial conditions [s0, e0, c0, p0] passed to ode45.
    % s = substrate
    % e = free enzyme
    % c = enzyme-substrate complex (ES)
    % p = product
    %
    % Reaction scheme: 
    % E + S <-> C <-> P + E (with reversible steps for both formation and dissociation)
    
    dx = zeros(size(x));
    
    dx(1) = -kpiu_1 * x(1) * x(2) + kmeno_1 * x(3);                         % ds/dt: Rate of change of substrate concentration
    dx(2) = kmeno_1 * x(3) - kpiu_1 * x(1) * x(2) + kpiu_2 * x(3) ...
            - kmeno_2 * x(2) * x(4);                                       % de/dt: Rate of change of free enzyme concentration
    dx(3) = kpiu_1 * x(1) * x(2) + kmeno_2 * x(2) * x(4) ...
            - kpiu_2 * x(3) - kmeno_1 * x(3);                              % dc/dt: Rate of change of the enzyme-substrate complex
    dx(4) = kpiu_2 * x(3) - kmeno_2 * x(2) * x(4);                         % dp/dt: Rate of change of product concentration
    
end