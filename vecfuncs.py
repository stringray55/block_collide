import numpy as np

def mag(v):
    """returns magnitude of vector v"""
    return np.sqrt(v.dot(v))

def unit(v):
    """returns the unit vector as np array in direction of vector v"""
    return v*mag(v)**(-1)

def R(r1, r2):
    """returns np array from r1 to r2"""
    return r2-r1

def Rhat(r1,r2):
    """returns unit np array from r1 to r2"""
    return unit(R(r1,r2))

def cross_mat(v):
    return np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
