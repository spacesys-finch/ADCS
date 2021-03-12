import matplotlib.pyplot as plt
import numpy as np

def plot_orbit2(pos1, vel1, t1, pos2, vel2, t2):
    '''
    Plots two satellites' orbit paths around the Earth.
    Plots two satellites' speeds versus time.
    
        Parameters:
            pos1 (array):   Satellite position vectors in the form [[x0, x1...], [y0, y1...], [z0, z1...]]
            vel1 (array):   Satellite velocity vectors in the form [[u0, u0...], [v0, v0...], [w0, w0...]]
            t1 (list):      Time vector associated with position and velocity vectors of satellite 1
            pos2 (array):   Satellite position vectors in the form [[x0, x0...], [y0, y0...], [z0, z0...]]
            vel2 (array):   Satellite velocity vectors in the form [[u0, u0...], [v0, v0...], [w0, w0...]]
            t2 (array):     Time vector associated with position and velocity vectors of satellite 2
    
        Returns:
            Figure 1 (3D):   Two satellites' paths around the Earth
            Figure 2 (2D):   Two satellites' speed versus time in km/s      
    '''
    
    x1, y1, z1 = pos1[0], pos1[1], pos1[2]
    u1, v1, w1 = vel1[0], vel1[1], vel1[2]

    x2, y2, z2 = pos2[0], pos2[1], pos2[2]
    u2, v2, w2 = vel2[0], vel2[1], vel2[2]

    fig1 = plt.figure()
    fig2 = plt.figure()
    
    # Spherical Coordinates
    theta = np.linspace(0, 2*np.pi)
    phi = np.linspace(0, np.pi)
    
    # Plot Earth (Radius = 6371 km)
    a = 6371 * np.outer(np.cos(theta), np.sin(phi))
    b = 6371 * np.outer(np.sin(theta), np.sin(phi))
    c = 6371 * np.outer(np.ones(np.size(theta)), np.cos(phi))

    earth = fig1.add_subplot(111, projection = '3d')
    earth.plot_surface(a, b, c, cmap='gist_earth_r', linewidth = 0.5)

    # Plot Position Data
    earth.plot3D(x1, y1, z1, c = 'red', label = 'Orbit Propagator')
    earth.plot3D(x2, y2, z2, c = 'blue', label = 'Poliastro')
    earth.set_title('Satellite Position (km)')
    earth.legend()
    #zoom = 15000
    #earth.set(xlim=(-zoom, zoom), ylim=(-zoom, zoom), zlim=(-zoom, zoom))
    
    # Plot Speed Data
    velocity = fig2.add_subplot(111)    
    speed1 = np.sqrt(u1**2 + v1**2 + w1**2)
    speed2 = np.sqrt(u2**2 + v2**2 + w2**2)
    velocity.plot(t1, speed1, c = 'red', label = 'Orbit Propagator')
    velocity.plot(t2, speed2, c = 'blue', label = 'Poliastro')
    velocity.set_title('Satellite Speed (km/s)')
    velocity.set_xlabel('Time (s)')
    velocity.set_ylabel('Speed (km/s)')
    velocity.legend()
    
    plt.show()

    return None
