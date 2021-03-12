import numpy as np
from scipy.integrate import solve_ivp

def orbit_from_equations(a,e,i,O,w,t0,tf,ts):
    '''
    Propagates orbit using direct equations.
    Returns position and velocity vectors of a satellite.
    
        Parameters:
            a (float):  Semi-major axis
                        Largest radius measured from centre of ellipse
            e (float):  Eccentricity
                        Indicates elliptical shape of orbit 
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
            position (array):   Satellite position vectors in the form [[x0, x1...], [y0, y1...], [z0, z1...]]
            velocity (array):   Satellite velocity vectors in the form [[u0, u0...], [v0, v0...], [w0, w0...]]
            t_list (array):     Time vector associated with position and velocity vectors
            y0 (array):         Initial condition vector in the form [x0, y0, z0, u0, v0, w0]
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
    
    # Initializa position and velocity vectors
    r = []
    v = []
    t_list = np.arange(t0,tf,ts)

    for t in t_list:
        M = np.sqrt(mu/a**3)*(t-t0) # Mean Anomaly
        E = M                       # First guess for Newton's method
        n = 0                       # Number of iterations
        n_max = 20                  # Max number of iterations (arbitrary)
        err = 1                     # Initial error value
    
        # Find Eccentric Anomaly (E) using Newton's method
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
        r.append(np.matmul(C_gp, [[R*np.cos(theta)],[R*np.sin(theta)], [0]]))
        v.append(np.matmul(C_gp, [[-np.sin(theta)],[np.cos(theta) + e], [0]]) * np.sqrt(mu/l))

    # Initial position and velocity vector
    y0 = np.concatenate((r[0], v[0]), axis=None)

    position = [[], [], []]
    velocity = [[], [], []]
    
    for i in r:
        position[0].append(float(i[0]))
        position[1].append(float(i[1]))
        position[2].append(float(i[2]))
    
    for j in v:
        velocity[0].append(float(j[0]))
        velocity[1].append(float(j[1]))
        velocity[2].append(float(j[2]))

    position = np.array(position)
    velocity = np.array(velocity)
    
    return position, velocity, t_list, y0

'''
References:
[1] Spacecraft Dynamics and Control (AER506) Course Notes
'''
