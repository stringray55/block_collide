import numpy as np

def get_delt_v(m,I_inv,r_c,vi,Li):
    """
    returns delta v. plane of collision must have normal in k_hat direction
    inputs:
    m) mass,
    I_inv) inverted inertia tensor,
    r_c) matrix rep of cross prod of point of collision relative to com. ,
    vi) initial com velocity
    Li) initial ang mom about com
    """
    a=0.5*(m+(m**2)*np.dot(I_inv,r_c).transpose()[2].dot(r_c.transpose()[2]))
    b=m*vi[2]+0.5*m*(np.dot(I_inv,Li).dot(r_c.transpose()[2])+np.dot(I_inv,r_c).transpose()[2].dot(Li))
    c=0
    return np.roots(np.array([a,b,c]))
