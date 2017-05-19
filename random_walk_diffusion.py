import numpy as np

def random_spin3d (size=1,seed=None) :
    """
    Generates a random three-dimensional vector, normalized on the
    unit sphere.
    """
    if seed is not None :
        np.random.seed(seed)
    phi = 2*np.pi*np.random.random(size=size)
    z = -1.0 + 2.0*np.random.random(size=size)
    s = np.sqrt (1.-z*z)
    return np.array ([s*np.cos(phi),s*np.sin(phi),z]).T

def random_walk_3D (nsteps,delta,origin=np.zeros(3)) :
    """
    Generates an array containing the 3D coordinates of a random walk starting
    at 'origin', lasting 'nsteps', with each step that is of length 'delta'.
    """
    trajectory = np.zeros((nsteps,3))
    trajectory[0,:] = origin
    trajectory[1:,:] = origin + np.cumsum(delta*random_spin3d(nsteps-1),axis=0)
    return trajectory
