function [soln] = SHO_solver(m, k, t0, tf, x0, v0)
    % Parameters: m, k
    % Time interval: t0 to tf
    % Initial conditions: x0, v0

    w = sqrt(k/m); % Angular Frequency

    % 2nd Order ODE: mx'' = -kx  ->  x'' = -(w^2)x
    syms t x(t)
    ode = diff(x,2) == -(w^2)*x;

    % Rewrite as a 1st Order ODE: X2' = -(w^2)X1
    V = odeToVectorField(ode); % V = [X2 ; -(w^2)X1]

    %Generate function: f = @(t,X)[X(2) ;-(w.^2).*X(1)] 
    f = matlabFunction(V, 'vars', {'t', 'Y'}); 

    % Use ode45 to solve the 1st Order ODE
    soln = ode45(f, [t0 tf], [x0 v0]);

    % Evaluate first component of solution at each point and plot
    fplot(@(x)deval(soln, x, 1), [t0, tf]);

end



