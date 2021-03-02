function [r, v, s] = plot_orbit(a,e,i,O,w,t0,tf)
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

    % Initialize vectors 
    r = zeros(3,tf); % Position vector
    v = zeros(3,tf); % Velocity vector
    s = zeros(1,tf); % Speed vector
    
    % Calculations
    for t = 1:tf              
        M = sqrt(mu/a^3)*(t-t0); % Mean anomaly
        E = M;                  % First guess for Newton's method
        n = 0;                  % Number of iterations
        n_max = 20;             % Max number of iterations
        err = 1;                % Initial error value
  
        % Find Eccentric Anomaly (E)
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
    l = a*(1-e^2);
    R = l/(1+e*cos(theta));        
    
    % Calculate Position and Velocity in ECI from Perifocal Frame
    r(:,t) = C_gp*[R*cos(theta); R*sin(theta); 0]; 
    v(:,t) = C_gp*[-sin(theta); cos(theta) + e; 0]*sqrt(mu/l); 
    
    % Calculate Speed (magnitude of velocity vector)
    s(t) = sqrt(v(1,t)^2 + v(2,t)^2 + v(3,t)^2);
    end
 
    % Plot Satellite Position
    figure
    plot3(r(1,:), r(2,:), r(3,:), '--r', 'LineWidth', 2)
    hold on
    [X,Y,Z] = sphere;
    surf(6371*X,6371*Y,6371*Z) % Radius of the Earth = 6371km
    axis equal
    title("Satellite Position (km)")
    %hold off
    
    % Plot Satellite Speed
    figure
    plot(1:tf, s, "r", 'LineWidth', 2)
    title ("Satellite Speed (km/s)")
    xlabel("Time (s)")
    ylabel("Speed (km/s)")   
end

% Resources:
% [1] Spacecraft Dynamics and Control (AER506) Course Notes
