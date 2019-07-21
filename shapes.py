import numpy as np
import rotate

class Node:
    def __init__(self,p,m):
        self.p=p
        self.m=m

class Tee:
    def __init__(self,long_nodes,short_nodes,spacing,init_w):
        """short_nodes is num nodes on each side of long axis"""
        self.long_nodes=long_nodes
        self.short_nodes=short_nodes
        self.spacing=spacing
        self.make_tee()

    def make_long_ax(self):
        self.tee=[Node(p=np.array([i*self.spacing,0,0]),m=1) for i in range(0,self.long_nodes)]

    def make_short_ax(self):
        x_pos=self.tee[-1].p[0]
        half1=[Node(p=np.array([x_pos,i*self.spacing,0]),m=1.5) for i in range(1,self.short_nodes+1)]
        half2=[Node(p=np.array([x_pos,-i*self.spacing,0]),m=1.5) for i in range(1,self.short_nodes+1)]
        self.tee=self.tee+half1+half2

    def make_tee(self):
        self.make_long_ax()
        self.make_short_ax()

class Cube:
    def __init__(self,nodes_per_side,spacing,corner,node_mass):
        self.nps=nodes_per_side
        self.nm=node_mass
        self.spacing=spacing
        self.corner=corner
        self.cube=self.make_cube()

    def make_line(self,corner):
        line=[Node(corner+np.array([i*self.spacing,0,0]),self.nm) for i in range(0,self.nps)]
        return line

    def make_sheet(self,corner):
        sheet=[]
        for i in range(0,self.nps):
            sheet=sheet+self.make_line(corner+np.array([0,i*self.spacing,0]))
        return sheet

    def make_cube(self):
        cube=[]
        for i in range(0,self.nps):
            cube=cube+self.make_sheet(self.corner+np.array([0,0,i*self.spacing]))
        return cube

class RectPrism:
    def __init__(self,xnodes,ynodes,znodes,spacing,corner,node_mass):
        self.xnodes=xnodes
        self.ynodes=ynodes
        self.znodes=znodes
        self.spacing=spacing
        self.corner=corner
        self.nm=node_mass
        self.make_rp

    def make_line(self,corner):
        line=[Node(corner+np.array([i*self.spacing,0,0]),self.nm) for i in range(0,self.xnodes)]
        return line

    def make_sheet(self,corner):
        sheet=[]
        for i in range(0,self.ynodes):
            sheet=sheet+self.make_line(corner+np.array([0,i*self.spacing,0]))
        return sheet

    def make_rp(self):
        rp=[]
        for i in range(0,self.znodes):
            rp=rp+self.make_sheet(self.corner+np.array([0,0,i*self.spacing]))
        return rp

class Block:
    def __init__(self,xnodes,ynodes,znodes,spacing,corner,node_mass=1):
        self.xnodes=xnodes
        self.ynodes=ynodes
        self.znodes=znodes
        self.spacing=spacing
        self.corner=corner
        self.nm=node_mass
        self.make_block()

    def make_line(self,corner):
        line=[Node(corner+np.array([i*self.spacing,0,0]),self.nm) for i in range(0,self.xnodes)]
        return line

    def make_corner_line(self,corner):
        line=[Node(corner+np.array([i*self.spacing,0,0]),self.nm) for i in range(1,self.xnodes-1)]
        corners=[Node(corner+np.array([i*self.spacing,0,0]),self.nm) for i in [0,self.xnodes-1]]
        return line,corners

    def make_sheet(self,corner):
        sheet=[]
        for i in range(0,self.ynodes):
            sheet=sheet+self.make_line(corner+np.array([0,i*self.spacing,0]))
        return sheet

    def make_corner_sheet(self,corner):
        sheet,corners=self.make_corner_line(corner)
        for i in range(1,self.ynodes-1):
            sheet=sheet+self.make_line(corner+np.array([0,i*self.spacing,0]))
        final_corner=corner+np.array([0,(self.ynodes-1)*self.spacing,0])
        sheet2,corners2=self.make_corner_line(final_corner)
        sheet=sheet+sheet2
        corners=corners+corners2
        return sheet,corners

    def make_block(self):
        block,corners=self.make_corner_sheet(self.corner)
        for i in range(1,self.znodes-1):
            block=block+self.make_sheet(self.corner+np.array([0,0,i*self.spacing]))
        final_corner=self.corner+np.array([0,0,(self.znodes-1)*self.spacing])
        top_sheet,top_corners=self.make_corner_sheet(final_corner)
        self.block=block+top_sheet
        self.corners=corners+top_corners

class Wheel:
    def __init__(self, len, rad, node_mass):
        self.len=len
        self.rad=rad
        self.nm=node_mass
        self.wheel=self.make_wheel()

    def make_wheel(self):
        wheel = [0 for i in range(0,26)]
        wheel[0]=Node(np.array([0,0,0]),self.nm)
        wheel[1]=Node(np.array([self.len,0,0]),self.nm)
        point=np.array([self.len,0,self.rad])
        for i in range(2,26):
            wheel[i]=Node(point,self.nm)
            point=rotate.point(point,w_hat=np.array([1,0,0]),deg=15)
        return wheel
