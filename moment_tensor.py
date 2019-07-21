import numpy as np

def get(lopm):
    """returns moment of inertia tensor for a lopm [[p, m],...]"""
    Ixx=0
    Iyy=0
    Izz=0
    Ixy=0
    Ixz=0
    Iyz=0
    for pm in lopm:
        xx=pm[0][1:3]
        yy=pm[0][0:3:2]
        zz=pm[0][0:2]
        Ixx+=xx.dot(xx)*pm[1]
        Iyy+=yy.dot(yy)*pm[1]
        Izz+=zz.dot(zz)*pm[1]
        Ixy+=pm[0][0]*pm[0][1]*pm[1]
        Ixz+=pm[0][0]*pm[0][2]*pm[1]
        Iyz+=pm[0][1]*pm[0][2]*pm[1]
    I=np.array([[Ixx,-Ixy,-Ixz],[-Ixy,Iyy,-Iyz],[-Ixz,-Iyz,Izz]])
    return I

def get_nodes(lon):
    """returns moment of inertia tensor for a lon"""
    Ixx=0
    Iyy=0
    Izz=0
    Ixy=0
    Ixz=0
    Iyz=0
    for n in lon:
        xx=n.p[1:3]
        yy=n.p[0:3:2]
        zz=n.p[0:2]
        Ixx+=xx.dot(xx)*n.m
        Iyy+=yy.dot(yy)*n.m
        Izz+=zz.dot(zz)*n.m
        Ixy+=n.p[0]*n.p[1]*n.m
        Ixz+=n.p[0]*n.p[2]*n.m
        Iyz+=n.p[1]*n.p[2]*n.m
    I=np.array([[Ixx,-Ixy,-Ixz],[-Ixy,Iyy,-Iyz],[-Ixz,-Iyz,Izz]])
    return I
