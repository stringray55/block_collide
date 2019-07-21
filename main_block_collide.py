#%%
import numpy as np
import matplotlib.pyplot as plt
import moment_tensor as mt
import animate
import plot3d as p3d
import rotate
import vecfuncs as vf
from shapes import Node
import center_of_mass as com
from collision import get_delt_v as gdv
import time
ar=np.array

class ComNode(Node):
    def __init__(self,p,m,v,g):
        super().__init__(p,m)
        self.v=v
        self.a=ar([0,0,g])

class Block:
    def __init__(self,xnodes,ynodes,znodes,spacing,corner,com_v,ang_mom,Elc,col_N,floor_tol,g,node_mass=1):
        self.xnodes=xnodes
        self.ynodes=ynodes
        self.znodes=znodes
        self.spacing=spacing
        self.corner=corner
        self.nm=node_mass
        self.col_N=col_N
        self.floor_tol=floor_tol
        self.make_block()
        self.com_node=self.get_com_node(com_v,g)
        self.princ_block=self.get_princ_block() #never changes, used with o_tens to get centered block
        self.centered_corners=self.get_centered_corners(self.corners,self.com_node)
        self.mom_tens=mt.get_nodes(self.princ_block+self.centered_corners)
        self.o_tens=np.array([[1,0,0],[0,1,0],[0,0,1]]) #actually the transpose of the matrix w/ columns representing principal axes
        self.ang_mom=ang_mom
        self.check_height=self.get_check_height()
        self.Elc=Elc
        self.kin_energies=[]
        self.pot_energies=[]
        self.times=[]

    def get_check_height(self):
        v=self.spacing*ar([self.xnodes,self.ynodes,self.znodes])
        return np.sqrt(v.dot(v))+1

    def make_line(self,corner):
        line=[Node(corner+ar([i*self.spacing,0,0]),self.nm) for i in range(0,self.xnodes)]
        return line

    def make_corner_line(self,corner):
        line=[Node(corner+ar([i*self.spacing,0,0]),self.nm) for i in range(1,self.xnodes-1)]
        corners=[Node(corner+ar([i*self.spacing,0,0]),self.nm) for i in [0,self.xnodes-1]]
        return line,corners

    def make_sheet(self,corner):
        sheet=[]
        for i in range(0,self.ynodes):
            sheet=sheet+self.make_line(corner+ar([0,i*self.spacing,0]))
        return sheet

    def make_corner_sheet(self,corner):
        sheet,corners=self.make_corner_line(corner)
        for i in range(1,self.ynodes-1):
            sheet=sheet+self.make_line(corner+ar([0,i*self.spacing,0]))
        final_corner=corner+ar([0,(self.ynodes-1)*self.spacing,0])
        sheet2,corners2=self.make_corner_line(final_corner)
        sheet=sheet+sheet2
        corners=corners+corners2
        return sheet,corners

    def make_block(self):
        block,corners=self.make_corner_sheet(self.corner)
        for i in range(1,self.znodes-1):
            block=block+self.make_sheet(self.corner+ar([0,0,i*self.spacing]))
        final_corner=self.corner+ar([0,0,(self.znodes-1)*self.spacing])
        top_sheet,top_corners=self.make_corner_sheet(final_corner)
        self.block=block+top_sheet
        self.corners=corners+top_corners

    def get_com_node(self,v,g):
        node=com.get_node(self.block+self.corners)
        return ComNode(node.p,node.m,v,g)

    def get_princ_block(self):
        return [Node(n.p-self.com_node.p,n.m) for n in self.block]

    def get_returned_block(self,centered_block):
        return [Node(n.p+self.com_node.p,n.m) for n in centered_block]

    def get_returned_corners(self,centered_corners):
        return [Node(n.p+self.com_node.p,n.m) for n in centered_corners]

    def get_centered_corners(self,corners,com_node):
        return [Node(n.p-com_node.p,n.m) for n in corners]

    def get_centered_block(self,o_tens):
        return [Node(np.dot(o_tens.transpose(),n.p),n.m) for n in self.princ_block]

    def get_mt(self,o_tens):
        return np.dot(np.dot(o_tens.transpose(),self.mom_tens),o_tens)

    def get_w(self,o_tens,ang_mom):
        mom_tens=self.get_mt(o_tens)
        return np.linalg.solve(mom_tens,ang_mom)

    def get_kin_energy(self,o_tens,ang_mom,com_node):
        mt=self.get_mt(o_tens)
        w=np.linalg.solve(mt,ang_mom)
        return 0.5*(com_node.m*com_node.v.dot(com_node.v)+0.5*w.dot(ang_mom))

    def get_pot_energy(self,o_tens,ang_mom,com_node):
        mt=self.get_mt(o_tens)
        w=np.linalg.solve(mt,ang_mom)
        return (-com_node.m*com_node.p[2]*com_node.a[2]+0.5*vf.mag(np.cross(w,ang_mom)))

    def get_corr_w(self,o_tens,ang_mom):
        k1=self.get_w(o_tens,ang_mom)
        k2=self.get_w(self.get_rot_o_tens(o_tens,k1/2),ang_mom)
        k3=self.get_w(self.get_rot_o_tens(o_tens,k2/2),ang_mom)
        k4=self.get_w(self.get_rot_o_tens(o_tens,k3),ang_mom)
        return (1/6)*(k1+2*k2+2*k3+k4)

    def get_translate(self,com_node,dt):
        com_p=com_node.p+(com_node.v+com_node.a*dt/2)*dt
        com_v=com_node.v+com_node.a*dt
        return ComNode(p=com_p,m=self.com_node.m,v=com_v,g=self.com_node.a[2])

    def get_rot_o_tens(self,o_tens,w):
        return ar([vf.unit(rotate.point(point=ax,w_hat=vf.unit(w),rad=vf.mag(w*self.dt))) for ax in o_tens])

    def get_rot_o_tens2(self,o_tens,w_hat,rad):
        return ar([vf.unit(rotate.point(point=ax,w_hat=w_hat,rad=rad)) for ax in o_tens])

    def get_rotate_corners(self,corners,w_hat,rad):
        return rotate.nodes(lon=corners,w_hat=w_hat,rad=rad)

    def nodes_below(self,com_node,corners):
        for n in corners:
            if n.p[2]+com_node.p[2]<=self.floor_tol:
                return True
        return False

    def get_nodes_below(self,step_com_node,step_corners,corners):
        nodes_below=[]
        for i in range(0,len(step_corners)):
            if step_corners[i].p[2]+step_com_node.p[2]<=self.floor_tol:
                nodes_below.append(corners[i])
        return nodes_below

    def get_nodes_above(self,step_com_node,step_corners,corners):
        nodes_above=[]
        for i in range(0,len(step_corners)):
            if step_corners[i].p[2]+step_com_node.p[2]>self.floor_tol:
                nodes_above.append(corners[i])
        return nodes_above

    def get_delt_v_L(self,r,o_tens,ang_mom,com_node):
        mt=self.get_mt(o_tens)
        I_inv=np.linalg.inv(mt)
        r_c=vf.cross_mat(r)
        Ei=self.get_kin_energy(o_tens,ang_mom,com_node)
        delt_v=gdv(m=com_node.m,I_inv=I_inv,r_c=r_c,vi=com_node.v,Li=ang_mom,Ei=Ei,Elc=self.Elc)
        if delt_v[1] > delt_v[0]:
            delt_v=ar([0,0,delt_v[1]])
        else:
            delt_v=ar([0,0,delt_v[0]])
        delt_L=com_node.m*np.dot(r_c,delt_v)
        return delt_v,delt_L

    def collision_check(self,step_com_node,step_corners,w_hat,rad):
        if step_com_node.p[2]<self.check_height:
            if not self.nodes_below(step_com_node,step_corners):
                self.com_node=step_com_node
                self.centered_corners=step_corners
                self.o_tens=self.get_rot_o_tens2(self.o_tens,w_hat,rad)
            else:
                (self.com_node,self.o_tens,self.centered_corners,self.ang_mom)=self.collide()
        else:
            self.com_node=step_com_node
            self.centered_corners=step_corners
            self.o_tens=self.get_rot_o_tens2(self.o_tens,w_hat,rad)

    def collide(self):
        self.dt_col=self.dt/self.col_N
        com_node=self.com_node
        corners=self.centered_corners
        o_tens=self.o_tens
        ang_mom=self.ang_mom
        print('collide initiated')
        for i in range(0,self.col_N):
            w=self.get_corr_w(o_tens,ang_mom)
            w_hat=vf.unit(w)
            rad=vf.mag(w)*self.dt_col
            step_com_node=self.get_translate(com_node,self.dt_col)
            step_corners=self.get_rotate_corners(corners,w_hat,rad)
            if not self.nodes_below(step_com_node,step_corners):
                com_node=step_com_node
                corners=step_corners
                o_tens=self.get_rot_o_tens2(o_tens,w_hat,rad)
            else:
                nodes_below=self.get_nodes_below(step_com_node,step_corners,corners)
                #print('nodes below: %s' %len(nodes_below))
                r=ar([0,0,0])
                count=0
                for n in nodes_below:
                    r=r+n.p
                    count+=1
                r=r/count
                #energy_before=self.get_kin_energy(o_tens,ang_mom,com_node)
                #print('energy before: %s' %energy_before)
                (delt_v,delt_L)=self.get_delt_v_L(r=r,o_tens=o_tens,ang_mom=ang_mom,com_node=com_node)
                #print('delt v: %s' %delt_v)
                #print('delt L: %s' %delt_L)
                com_node.v=com_node.v+delt_v
                ang_mom=ang_mom+delt_L
                #energy_after=self.get_kin_energy(o_tens,ang_mom,com_node)
                #print('energy after: %s' %energy_after)
                w=self.get_corr_w(o_tens,ang_mom)
                w_hat=vf.unit(w)
                rad=vf.mag(w)*self.dt_col
                step_com_node=self.get_translate(com_node,self.dt_col)
                step_corners=self.get_rotate_corners(corners,w_hat,rad)
                if not self.nodes_below(step_com_node,step_corners):
                    print('all good')
                    com_node=step_com_node
                    corners=step_corners
                    o_tens=self.get_rot_o_tens2(o_tens,w_hat,rad)
                else:
                    print('it happened')
                    nodes_below=self.get_nodes_below(step_com_node,step_corners,corners)
                    nodes_above=self.get_nodes_above(step_com_node,step_corners,corners)
                    r=ar([0,0,0])
                    count=0
                    for n in nodes_below:
                        r=r+n.p
                        count+=1
                    r=r/count
                    #energy_before=self.get_kin_energy(o_tens,ang_mom,com_node)
                    #print('energy before: %s' %energy_before)
                    (delt_v,delt_L)=self.get_delt_v_L(r=r,o_tens=o_tens,ang_mom=ang_mom,com_node=com_node)
                    com_node.v=com_node.v+delt_v
                    ang_mom=ang_mom+delt_L
                    #energy_after=self.get_kin_energy(o_tens,ang_mom,com_node)
                    #print('energy after: %s' %energy_after)
                    w=self.get_corr_w(o_tens,ang_mom)
                    w_hat=vf.unit(w)
                    rad=vf.mag(w)*self.dt_col
                    step_com_node=self.get_translate(com_node,self.dt_col)
                    step_corners=self.get_rotate_corners(corners,w_hat,rad)
                    if not self.nodes_below(step_com_node,step_corners):
                        print('all good 2')
                        com_node=step_com_node
                        corners=step_corners
                        o_tens=self.get_rot_o_tens2(o_tens,w_hat,rad)
                    else:
                        print('it happened again')
                        print('com node v is 10')
                        com_node.v[2]=10
                        com_node=self.get_translate(com_node,self.dt_col)
                        corners=self.get_rotate_corners(corners,w_hat,rad)
                        o_tens=self.get_rot_o_tens2(o_tens,w_hat,rad)
        print('collide complete')
        return com_node,o_tens,corners,ang_mom

    def step(self,times):
        for i in range(0,times):
            w=self.get_corr_w(self.o_tens,self.ang_mom)
            w_hat=vf.unit(w)
            rad=vf.mag(w)*self.dt
            step_com_node=self.get_translate(self.com_node,self.dt)
            step_corners=self.get_rotate_corners(self.centered_corners,w_hat,rad)
            self.collision_check(step_com_node,step_corners,w_hat,rad)
            self.t+=self.dt

    def print_block(self,limits):
        p3d.scatter_nodes(self.block+self.corners,xlim=limits,ylim=limits,zlim=(0,50))

    def print_frame(self,limits):
        return p3d.scatter_nodes_frame(self.block+self.corners,xlim=limits,ylim=limits,zlim=(0,100),elev=10,azim=(45+self.camera_rot*self.t))

    def print_frame2(self,block_ax,kin_energy_ax,pot_energy_ax,block_lim):
        block_ax=p3d.scatter_nodes_frame2(lon=self.block+self.corners,ax=block_ax,xlim=block_lim,ylim=block_lim,zlim=(0,100),elev=10,azim=(45+self.camera_rot*self.t))
        kin_energy_ax.plot(self.times,self.kin_energies,lw=2)
        kin_energy_ax.scatter(self.times[-1],self.kin_energies[-1],color='red',lw=5)
        plt.xlabel('Time')
        plt.ylabel('Energy')
        pot_energy_ax.plot(self.times,self.pot_energies,lw=2,color='black')
        kin_energy_ax.scatter(self.times[-1],self.pot_energies[-1],color='red',lw=5)
        plt.ylabel('Energy')
        return block_ax,kin_energy_ax,pot_energy_ax

    def anim_init(self):
        self.grid=plt.GridSpec(5,3,wspace=0.1,hspace=0.1)
        self.block_ax=self.fig.add_subplot(self.grid[:3,:3],projection='3d',xlim=(-50,50),ylim=(-50,50),zlim=(0,100))
        self.kin_energy_ax=self.fig.add_subplot(self.grid[3,:3],xlim=(0,self.tf))
        self.pot_energy_ax=self.fig.add_subplot(self.grid[4,:3],xlim=(0,self.tf))
        return self.block_ax,self.kin_energy_ax,self.pot_energy_ax

    def frame(self,f):
        if f%30==0:
            print("sim time: %s" %(self.t))
        self.kin_energies.append(self.get_kin_energy(self.o_tens,self.ang_mom,self.com_node))
        self.times.append(self.t)
        self.fig.clf()
        centered_block=self.get_centered_block(self.o_tens)
        self.block=self.get_returned_block(centered_block)
        self.corners=self.get_returned_corners(self.centered_corners)
        ax=self.print_frame(limits=(-50,50))
        self.step(self.steps_per_frame)
        return ax

    def frame2(self,f):
        if f%30==0:
            print("sim time: %s" %(self.t))
        self.kin_energies.append(self.get_kin_energy(self.o_tens,self.ang_mom,self.com_node))
        self.pot_energies.append(self.get_pot_energy(self.o_tens,self.ang_mom,self.com_node))
        self.times.append(self.t)
        self.block_ax.remove()
        self.kin_energy_ax.remove()
        self.pot_energy_ax.remove()
        centered_block=self.get_centered_block(self.o_tens)
        self.block=self.get_returned_block(centered_block)
        self.corners=self.get_returned_corners(self.centered_corners)
        self.block_ax=self.fig.add_subplot(self.grid[:3,:3],projection='3d',xlim=(-50,50),ylim=(-50,50),zlim=(0,100))
        self.pot_energy_ax=self.fig.add_subplot(self.grid[4,:3],xlim=(0,self.tf))
        self.kin_energy_ax=self.fig.add_subplot(self.grid[3,:3],sharex=self.pot_energy_ax)
        (self.block_ax,self.kin_energy_ax,self.pot_energy_ax)=self.print_frame2(self.block_ax,self.kin_energy_ax,self.pot_energy_ax,block_lim=(-50,50))
        self.step(self.steps_per_frame)
        return self.block_ax,self.kin_energy_ax

    def animate_block(self,save_name,size,dt,tf,steps_per_frame,camera_rot):
        self.t=0
        self.dt=dt
        self.tf=tf
        self.steps_per_frame=steps_per_frame
        self.camera_rot=camera_rot
        self.fig=plt.figure(figsize=(10,10))
        animate.do_60fps(save_name=save_name,fig=self.fig,anim_init=animate.anim_init_3d,frame=self.frame,tot_frames=int(tf/(dt*steps_per_frame)))

    def animate_block2(self,save_name,size,dt,tf,steps_per_frame,camera_rot):
        self.t=0
        self.dt=dt
        self.tf=tf
        self.steps_per_frame=steps_per_frame
        self.camera_rot=camera_rot
        self.fig=plt.figure(figsize=(9,15))
        animate.do_60fps(save_name=save_name,fig=self.fig,anim_init=self.anim_init,frame=self.frame2,tot_frames=int(tf/(dt*steps_per_frame)))
