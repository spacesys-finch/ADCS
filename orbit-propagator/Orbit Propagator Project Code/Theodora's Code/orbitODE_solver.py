import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

def orbitODE_solver(a,e,i,O,w,t0,tf,ts):
    '''
    Integrates ODE given by the two-body problem.
    Returns position and velocity vectors of a satellite.
    
        Parameters:
            a (float):  Semi-major axis
                        Largest radius measured from centre of ellipse
            e (float):  Eccentricity
                        Indicates elliptical shape of orbit (0 <= e < 1, where e = 0 is a circle)
            i (float):  Inclination
                        Angle between elliptical plane and reference plane
            O (float):  Right Ascension of the Ascending Node (RAAN) / Longitude of Ascending Node
                        Angle between reference direction and ascending node in travelling direction of satellite from ascending node
                        Generally represented by Omega
            w (float):  Argument of periapsis / Argument of perigee
                        Angle between ascending node and semi-major axis
            t0 (int):   Time of periapsis passage
            tf (int):   Length of orbit path
            ts (float): Desired time step between points
    
        Returns:
            position (list):   Satellite position vectors in the form [[x], [y], [z]]
            velocity (list):   Satellite velocity vectors in the form [[u], [v], [w]]
            t_list (list):     Time vector associated with position and velocity vectors
    '''

    # Define Constants
    G = 6.67430e-20     # Gravitational Constant (km)
    ME = 5.9725e24      # Mass of Earth (kg)
    mu = G*ME           # Gravitational Parameter for Earth (km*kg)

    # Define Rotation Matrices
    C1_i = np.array([[1,0,0],[0, np.cos(i), np.sin(i)], [0,-np.sin(i), np.cos(i)]])
    C3_w = np.array([[np.cos(w),np.sin(w),0],[-np.sin(w), np.cos(w), 0], [0,0, 1]])
    C3_O = np.array([[np.cos(O),np.sin(O),0],[-np.sin(O), np.cos(O), 0], [0,0, 1]])
    
    # Move from perifocal to ECI frame
    C_gp = np.linalg.inv(np.matmul(np.matmul(C3_w, C1_i), C3_O))
    
    # Find Eccentric Anomaly (E) using Newton's method
    t = t0
    M = np.sqrt(mu/a**3)*(t-t0) # Mean Anomaly
    E = M                       # First guess for Newton's method
    n = 0                       # Number of iterations
    n_max = 20                  # Max number of iterations (arbitrary)
    err = 1                     # Initial error value

    while(err>=1e-6 and n<=n_max):
        f = E - e*np.sin(E) - M
        ff = 1 - e*np.cos(E)
        E_new = E - f/ff
        err = abs(E-E_new)
        E = E_new 
        n = n+1

    # Calculate True Anomaly (theta) and Radial Distance (R)
    theta = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
    l = a*(1-e**2)
    R = l/(1+e*np.cos(theta))

    # Find r0 and v0
    r0 = np.matmul(C_gp, [[R*np.cos(theta)],[R*np.sin(theta)], [0]])
    v0 = np.matmul(C_gp, [[-np.sin(theta)],[np.cos(theta) + e], [0]]) * np.sqrt(mu/l)

    # Initial position and velocity vector
    y0 = np.concatenate((r0, v0), axis=None)        
    
    # Define ODE to be solved
    def ode(t,x):  
        r = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
        q = -mu/r**3
        A = np.array([[0,0,0,1,0,0], [0,0,0,0,1,0], [0,0,0,0,0,1], [q, 0,0,0,0,0], [0,q,0,0,0,0], [0,0,q,0,0,0]])
        dx = np.matmul(A,x)
        return dx  
    
    # Time list
    t_list = np.arange(t0,tf,ts)

    # Solve ODE using 8th order integration method
    data = solve_ivp(ode,[t0,tf], y0,'DOP853', t_list)  

    position = np.concatenate(([data.y[0]],[data.y[1]],[data.y[2]]), axis=0)
    velocity = np.concatenate(([data.y[3]],[data.y[4]],[data.y[5]]), axis=0)
    
    return position, velocity, t_list

'''
References:
[1] Spacecraft Dynamics and Control (AER506) Course Notes
[2] H. D. Curtis, Orbital Mechanics for Engineering Students (Appendix D.6)
'''