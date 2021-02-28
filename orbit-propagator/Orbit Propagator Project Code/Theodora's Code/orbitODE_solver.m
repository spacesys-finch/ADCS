function orbitODE_solver(a,e,i,O,w,t0,tf)
    % Constants
    G = 6.67430e-20;    % Gravitational Constant (in km)
    ME = 5.9725e24;     % Mass of Earth
    mu = G*ME;          % Gravitational Parameter for Earth
    
    % Rotation Matrices
    C1_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
    C3_w = [cos(w), sin(w), 0; -sin(w), cos(w), 0; 0, 0, 1];
    C3_O = [cos(O), sin(O), 0; -sin(O), cos(O), 0; 0, 0, 1];

    % Move from perifocal to ECI frame
    C_gp = (C3_w*C1_i*C3_O)'; % Transpose of C_pg (equal to inverse)

    % Eccentric Anomaly (E)
    t = t0;
    M = sqrt(mu/a^3)*(t-t0); % Mean anomaly
    E = M;                  % First guess for Newton's method
    n = 0;                  % Number of iterations
    n_max = 20;             % Max number of iterations
    err = 1;                % Initial error value
    
    while (err>=1e-6 && n<=n_max)   % Continue if error is too large or reached max iterations
        f = E - e*sin(E) - M;       % Function
        ff = 1 - e*cos(E);          % Function derivative
        E_new = E - f/ff;           % Newton's Method
        err = abs(E-E_new);         % Check error value
        E = E_new;                  % Update E
        n = n+1;                    % Count iteration    
    end
    
    % Calculate True Anomaly (theta) and Radial Distance (R)
    theta = 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); 
    l = a*(1-e^2);                  % Semilatus rectum
    R = l/(1+e*cos(theta)); 
   
    % Find r0 and v0
    r0 = C_gp*[R*cos(theta); R*sin(theta); 0]; 
    v0 = C_gp*[-sin(theta); cos(theta) + e; 0]*sqrt(mu/l);
    
    y0 = [r0 v0]; % Initial position and velocity vector
    [t,y] = ode45(@ode, [t0 tf], y0);
    
    output 
    
    return
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function dydt = ode(t,f)
        x = f(1);
        y = f(2);
        z = f(3);
        vx = f(4);
        vy = f(5);
        vz = f(6);
        r = norm([x y z]); % Magnitude of position vector
        ax = -mu*x/r^3;
        ay = -mu*y/r^3;
        az = -mu*z/r^3;
        dydt = [vx vy vz ax ay az]';
        end 
    
    function output
        figure
        plot3(y(:,1), y(:,2), y(:,3), "r--", 'LineWidth', 2)
        hold on
        [X,Y,Z] = sphere;
        surf(6371*X,6371*Y,6371*Z) % Radius of the Earth = 6371km
        axis equal
        title("Satellite Position (km)")
    end 
end

