import numpy as np
from shapes import Node

class Quat:
    def __init__(self,real,vect):
        self.real=real
        self.vect=vect

    def mult(self,q2):
        real=self.real*q2.real-self.vect.dot(q2.vect)
        vect=self.real*q2.vect+q2.real*self.vect+np.cross(q2.vect,self.vect)
        return Quat(real,vect)

    def inverse(self):
        real=self.real
        vect=-self.vect
        return Quat(real,vect)


def rot_1_point(q,p):
    """inputs"""
    rot=q.inverse().mult(p).mult(q)
    return rot.vect

def point(point,w_hat,rad):
    q=Quat(np.cos(rad/2),np.sin(rad/2)*w_hat)
    p=Quat(0,point)
    return rot_1_point(q,p)

def points(lop,w_hat,rad):
    q=Quat(np.cos(rad/2),np.sin(rad/2)*w_hat)
    p_list=[Quat(0,point) for point in lop]
    rot_points=[rot_1_point(q,p) for p in p_list]
    return rot_points

def nodes(lon,w_hat,rad):
    """returns list of rotated nodes without changing original lon"""
    q=Quat(np.cos(rad/2),np.sin(rad/2)*w_hat)
    rot_nodes=[Node(rot_1_point(q,Quat(0,n.p)),n.m) for n in lon]
    return rot_nodes
