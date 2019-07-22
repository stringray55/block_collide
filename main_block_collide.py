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
    def __init__(self,xnodes,ynodes,znodes,spacing,corner,com_v,ang_mom,vlc,col_N,floor_tol,g,ground_E_tol,fric,node_mass=1):
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
        self.vlc=vlc
        self.trans_KE=[]
        self.rot_KE=[]
        self.PE=[]
        self.TE=[]
        self.times=[]
        self.on_ground=False
        self.ground_E_tol=ground_E_tol
        self.fric=fric #high values for quick slowdown. must be less than mag(ang_mom)/dt or else ang_mom will oscillate

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

    def get_trans_KE(self,com_node):
        return 0.5*(com_node.m*com_node.v.dot(com_node.v))

    def get_rot_KE(self,o_tens,ang_mom):
        mt=self.get_mt(o_tens)
        w=np.linalg.solve(mt,ang_mom)
        return 0.5*w.dot(ang_mom)

    def get_alt_rot_KE(self,o_tens,ang_mom):
        #for use in computing energy minus z rotational energy for on ground detection
        mt=self.get_mt(o_tens)
        w=np.linalg.solve(mt,ang_mom)
        return 0.5*w[0:2].dot(ang_mom[0:2])

    def get_KE(self,o_tens,ang_mom,com_node):
        return self.get_trans_KE(com_node)+self.get_rot_KE(o_tens,ang_mom)

    def get_PE(self,com_node):
        return -com_node.m*(com_node.p[2])*com_node.a[2]

    def get_TE(self,o_tens,ang_mom,com_node):
        return self.get_trans_KE(com_node)+self.get_rot_KE(o_tens,ang_mom)+self.get_PE(com_node)

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

    def get_nodes_below2(self,step_com_node,step_corners):
        nodes_below=[]
        for i in range(0,len(step_corners)):
            if step_corners[i].p[2]+step_com_node.p[2]<=self.floor_tol:
                nodes_below.append(step_corners[i])
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
        delt_v=gdv(m=com_node.m,I_inv=I_inv,r_c=r_c,vi=com_node.v,Li=ang_mom)
        if delt_v[1] > delt_v[0]:
            delt_v=ar([0,0,self.vlc*delt_v[1]])
        else:
            delt_v=ar([0,0,self.vlc*delt_v[0]])
        delt_L=com_node.m*np.dot(r_c,delt_v)
        return delt_v,delt_L

    def do_collide_before(self,step_com_node,step_corners,corners,o_tens,ang_mom,com_node):
        #does collision as if happening at step before corner goes below floor
        nodes_below=self.get_nodes_below(step_com_node,step_corners,corners)
        r=ar([0,0,0])
        count=0
        for n in nodes_below:
            r=r+n.p
            count+=1
        r=r/count
        (delt_v,delt_L)=self.get_delt_v_L(r=r,o_tens=o_tens,ang_mom=ang_mom,com_node=com_node)
        com_node.v=com_node.v+delt_v
        ang_mom=ang_mom+delt_L
        w=self.get_corr_w(o_tens,ang_mom)
        w_hat=vf.unit(w)
        rad=vf.mag(w)*self.dt_col
        step_com_node=self.get_translate(com_node,self.dt_col)
        step_corners=self.get_rotate_corners(corners,w_hat,rad)
        return step_com_node,step_corners,ang_mom,com_node,w_hat,rad

    def do_collide_after(self,step_com_node,step_corners,corners,step_o_tens,ang_mom,):
        #does collision at step in which corner is below floor
        nodes_below=self.get_nodes_below2(step_com_node,step_corners)
        r=ar([0,0,0])
        count=0
        for n in nodes_below:
            r=r+n.p
            count+=1
        r=r/count
        (delt_v,delt_L)=self.get_delt_v_L(r=r,o_tens=step_o_tens,ang_mom=ang_mom,com_node=step_com_node)
        step_com_node.v=step_com_node.v+delt_v
        ang_mom=ang_mom+delt_L
        w=self.get_corr_w(step_o_tens,ang_mom)
        w_hat=vf.unit(w)
        rad=vf.mag(w)*self.dt_col
        step_com_node=self.get_translate(step_com_node,self.dt_col)
        step_corners=self.get_rotate_corners(corners,w_hat,rad)
        return step_com_node,step_corners,ang_mom,w_hat,rad

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
            if not self.on_ground:
                if not self.nodes_below(step_com_node,step_corners):
                    com_node=step_com_node
                    corners=step_corners
                    o_tens=self.get_rot_o_tens2(o_tens,w_hat,rad)
                else:
                    (step_com_node,step_corners,ang_mom,com_node,w_hat,rad)=self.do_collide_before(step_com_node,step_corners,corners,o_tens,ang_mom,com_node)
                    if not self.nodes_below(step_com_node,step_corners):
                        print('all good')
                        com_node=step_com_node
                        corners=step_corners
                        o_tens=self.get_rot_o_tens2(o_tens,w_hat,rad)
                    else:
                        if self.get_trans_KE(com_node)+self.get_alt_rot_KE(o_tens,ang_mom)<=self.ground_E_tol:
                            print('on ground, energy less than ground tol')
                            com_node.v=ar([0,0,0])
                            com_node.a=ar([0,0,0])
                            ang_mom=ar([0,0,ang_mom[2]])
                            self.on_ground=True
                        else:
                            step_o_tens=self.get_rot_o_tens2(o_tens,w_hat,rad)
                            (step_com_node,step_corners,ang_mom,w_hat,rad)=self.do_collide_after(step_com_node,step_corners,corners,step_o_tens,ang_mom)
                            if not self.nodes_below(step_com_node,step_corners):
                                print('all good 2')
                                com_node=step_com_node
                                corners=step_corners
                                o_tens=self.get_rot_o_tens2(o_tens,w_hat,rad)
                            else:
                                print('fail, on ground')
                                com_node.v=ar([0,0,0])
                                com_node.a=ar([0,0,0])
                                ang_mom=ar([0,0,ang_mom[2]])
                                self.on_ground=True
            else:
                w=self.get_corr_w(self.o_tens,self.ang_mom)
                rad=vf.mag(w)*self.dt
                if rad!=0:
                    w_hat=vf.unit(w)
                    corners=self.get_rotate_corners(self.centered_corners,w_hat,rad)
                    o_tens=self.get_rot_o_tens2(self.o_tens,w_hat,rad)
                    ang_mom=self.ang_mom-self.fric*self.ang_mom*self.dt
                else:
                    print('done spinning')
                    break
        print('collide complete')
        return com_node,o_tens,corners,ang_mom

    def step(self,times):
        for i in range(0,times):
            if not self.on_ground:
                w=self.get_corr_w(self.o_tens,self.ang_mom)
                w_hat=vf.unit(w)
                rad=vf.mag(w)*self.dt
                step_com_node=self.get_translate(self.com_node,self.dt)
                step_corners=self.get_rotate_corners(self.centered_corners,w_hat,rad)
                self.collision_check(step_com_node,step_corners,w_hat,rad)
            else:
                w=self.get_corr_w(self.o_tens,self.ang_mom)
                rad=vf.mag(w)*self.dt
                if rad!=0:
                    w_hat=vf.unit(w)
                    self.centered_corners=self.get_rotate_corners(self.centered_corners,w_hat,rad)
                    self.o_tens=self.get_rot_o_tens2(self.o_tens,w_hat,rad)
                    self.ang_mom=self.ang_mom-self.fric*self.ang_mom*self.dt
            self.t+=self.dt


    def print_block(self,limits):
        p3d.scatter_nodes(self.block+self.corners,xlim=limits,ylim=limits,zlim=(0,50))

    def print_frame(self,limits):
        return p3d.scatter_nodes_frame(self.block+self.corners,xlim=limits,ylim=limits,zlim=(0,100),elev=10,azim=(45+self.camera_rot*self.t))

    def print_frame_energy(self,block_ax,energy_ax,block_lim):
        block_ax=p3d.scatter_nodes_frame2(lon=self.block+self.corners,ax=block_ax,xlim=block_lim,ylim=block_lim,zlim=(0,100),elev=10,azim=(45+self.camera_rot*self.t))
        energy_ax=self.get_energy_ax(energy_ax)
        return block_ax,energy_ax

    def get_energy_ax(self,energy_ax):
        energy_ax.plot(self.times,self.trans_KE,lw=2,color='blue',label='translational kinetic')
        energy_ax.scatter(self.times[-1],self.trans_KE[-1],color='red',lw=5)
        energy_ax.plot(self.times,self.rot_KE,lw=2,color='black',label='rotational kinetic')
        energy_ax.scatter(self.times[-1],self.rot_KE[-1],color='green',lw=5)
        energy_ax.plot(self.times,self.PE,lw=2,color='orange',label='potential')
        energy_ax.scatter(self.times[-1],self.PE[-1],color='brown',lw=5)
        energy_ax.plot(self.times,self.TE,lw=2,color='purple',label='total')
        energy_ax.scatter(self.times[-1],self.TE[-1],color='grey',lw=5)
        energy_ax.legend(loc='upper right')
        plt.xlabel('Time')
        plt.ylabel('Energy')
        return energy_ax

    def append_energies(self):
        trans_KE=self.get_trans_KE(self.com_node)
        self.trans_KE.append(trans_KE)
        rot_KE=self.get_rot_KE(self.o_tens,self.ang_mom)
        self.rot_KE.append(rot_KE)
        PE=self.get_PE(self.com_node)
        self.PE.append(PE)
        TE=trans_KE+rot_KE+PE
        self.TE.append(TE)
        self.times.append(self.t)

    def anim_init(self):
        self.grid=plt.GridSpec(6,4,wspace=0.1,hspace=0.1)
        self.block_ax=self.fig.add_subplot(self.grid[:4,:4],projection='3d',xlim=(-50,50),ylim=(-50,50),zlim=(0,100))
        self.energy_ax=self.fig.add_subplot(self.grid[4:6,:4],xlim=(0,self.tf))
        return self.block_ax,self.energy_ax

    def frame(self,f):
        if f%30==0:
            print("sim time: %s" %(self.t))
        self.fig.clf()
        centered_block=self.get_centered_block(self.o_tens)
        self.block=self.get_returned_block(centered_block)
        self.corners=self.get_returned_corners(self.centered_corners)
        ax=self.print_frame(limits=(-50,50))
        self.step(self.steps_per_frame)
        return ax

    def frame_energy(self,f):
        if f%30==0:
            print("sim time: %s" %(self.t))
        self.append_energies()
        self.block_ax.remove()
        self.energy_ax.remove()
        centered_block=self.get_centered_block(self.o_tens)
        self.block=self.get_returned_block(centered_block)
        self.corners=self.get_returned_corners(self.centered_corners)
        self.block_ax=self.fig.add_subplot(self.grid[:4,:4],projection='3d',xlim=(-50,50),ylim=(-50,50),zlim=(0,100))
        self.energy_ax=self.fig.add_subplot(self.grid[4:6,:4],xlim=(0,self.tf))
        (self.block_ax,self.energy_ax)=self.print_frame_energy(self.block_ax,self.energy_ax,block_lim=(-50,50))
        self.step(self.steps_per_frame)
        return self.block_ax,self.energy_ax

    def animate_block(self,save_name,size,dt,tf,steps_per_frame,camera_rot):
        self.t=0
        self.dt=dt
        self.tf=tf
        self.steps_per_frame=steps_per_frame
        self.camera_rot=camera_rot
        self.fig=plt.figure(figsize=(10,10))
        animate.do_60fps(save_name=save_name,fig=self.fig,anim_init=animate.anim_init_3d,frame=self.frame,tot_frames=int(tf/(dt*steps_per_frame)))

    def animate_block_energy(self,save_name,size,dt,tf,steps_per_frame,camera_rot):
        self.t=0
        self.dt=dt
        self.tf=tf
        self.steps_per_frame=steps_per_frame
        self.camera_rot=camera_rot
        self.fig=plt.figure(figsize=(8,12))
        animate.do_60fps(save_name=save_name,fig=self.fig,anim_init=self.anim_init,frame=self.frame_energy,tot_frames=int(tf/(dt*steps_per_frame)))
