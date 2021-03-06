import matplotlib.pyplot as plt
import numpy as np

def plot_orbit(position, velocity, t_list):
    '''
    Plots a satellite's orbit path around the Earth.
    Plots satellite speed versus time.
    
        Parameters:
            position (list):   Satellite position vectors in the form [[x], [y], [z]]
            velocity (list):   Satellite velocity vectors in the form [[u], [v], [w]]
            t_list (list):     Time vector associated with position and velocity vectors
    
        Returns:
            Figure 1 (3D):   Satellite path around the Earth
            Figure 2 (2D):   Satellite speed versus time in km/s      
    '''
    x, y, z = position[0], position[1], position[2]
    u, v, w = velocity[0], velocity[1], velocity[2]

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
    earth.scatter(x,y,z, c = 'red')
    earth.set_title('Satellite Position (km)')
    #zoom = 15000
    #earth.set(xlim=(-zoom, zoom), ylim=(-zoom, zoom), zlim=(-zoom, zoom))
    
    # Plot Speed Data
    velocity = fig2.add_subplot(111)    
    speed = np.sqrt(u**2 + v**2 + w**2)
    velocity.plot(t_list, speed)
    velocity.set_title('Satellite Speed (km/s)')
    velocity.set_xlabel('Time (s)')
    velocity.set_ylabel('Speed (km/s)')
    
    plt.show()

    return None
