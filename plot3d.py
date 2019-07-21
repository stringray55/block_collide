import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def scatter(lop,size=8,xlim=(-10,10),ylim=(-10,10),zlim=(-10,10),lw=1,elev=35,azim=45):
    """plots lop (list of positions) in 3D"""
    x=[n.p[0] for n in lon]
    y=[n.p[1] for n in lon]
    z=[n.p[2] for n in lon]
    fig=plt.figure(figsize=(size,size))
    ax=plt.axes(projection='3d',xlim=xlim,ylim=ylim,zlim=zlim)
    ax.elev=elev
    ax.azim=azim
    ax.scatter(x,y,z,lw=lw)
    plt.xlabel('x')
    plt.ylabel('y')

def scatter_frame(lop,xlim,ylim,zlim,lw=1,elev=35,azim=45):
    """plots lop (list of positions) in 3D"""
    x=[n.p[0] for n in lon]
    y=[n.p[1] for n in lon]
    z=[n.p[2] for n in lon]
    ax=plt.axes(projection='3d',xlim=xlim,ylim=ylim,zlim=zlim)
    ax.elev=elev
    ax.azim=azim
    ax.scatter(x,y,z,lw=lw)
    plt.xlabel('x')
    plt.ylabel('y')
    return ax

def scatter_nodes(lon,size=8,xlim=(-10,10),ylim=(-10,10),zlim=(-10,10),lw=1,elev=35,azim=45):
    """plots lop (list of positions) in 3D"""
    x=[n.p[0] for n in lon]
    y=[n.p[1] for n in lon]
    z=[n.p[2] for n in lon]
    fig=plt.figure(figsize=(size,size))
    ax=plt.axes(projection='3d',xlim=xlim,ylim=ylim,zlim=zlim)
    ax.elev=elev
    ax.azim=azim
    ax.scatter(x,y,z,lw=lw)
    plt.xlabel('x')
    plt.ylabel('y')

def scatter_nodes_frame(lon,xlim=(-10,10),ylim=(-10,10),zlim=(-10,10),lw=1,elev=35,azim=45):
    """plots lop (list of positions) in 3D"""
    x=[n.p[0] for n in lon]
    y=[n.p[1] for n in lon]
    z=[n.p[2] for n in lon]
    ax=plt.axes(projection='3d',xlim=xlim,ylim=ylim,zlim=zlim)
    ax.elev=elev
    ax.azim=azim
    ax.scatter(x,y,z,lw=lw)
    plt.xlabel('x')
    plt.ylabel('y')
    return ax

def scatter_nodes_frame2(lon,ax,xlim=(-10,10),ylim=(-10,10),zlim=(-10,10),lw=1,elev=35,azim=45):
    """plots lop (list of positions) in 3D"""
    x=[n.p[0] for n in lon]
    y=[n.p[1] for n in lon]
    z=[n.p[2] for n in lon]
    ax.elev=elev
    ax.azim=azim
    ax.scatter(x,y,z,lw=lw)
    plt.xlabel('x')
    plt.ylabel('y')
    return ax

def line_nodes_frame(lon,xlim=(-10,10),ylim=(-10,10),zlim=(-10,10),lw=1,elev=35,azim=45):
    """plots lop (list of positions) in 3D"""
    x=[n.p[0] for n in lon]
    y=[n.p[1] for n in lon]
    z=[n.p[2] for n in lon]
    ax=plt.axes(projection='3d',xlim=xlim,ylim=ylim,zlim=zlim)
    ax.elev=elev
    ax.azim=azim
    ax.plot(x,y,z,lw=lw)
    plt.xlabel('x')
    plt.ylabel('y')
    return ax
