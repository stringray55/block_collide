import numpy as np
from shapes import Node

def get(lon):
    com=np.array([0,0,0])
    tot_mass=0
    for n in lon:
        com=com+n.p*n.m
        tot_mass+=n.m
    return com/tot_mass

def get_node(lon):
    com=np.array([0,0,0])
    tot_mass=0
    for n in lon:
        com=com+n.p*n.m
        tot_mass+=n.m
    return Node(com/tot_mass,tot_mass)
