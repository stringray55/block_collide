import main_block_collide as mbc
import matplotlib.pyplot as plt
import numpy as np
ar=np.array

block=mbc.Block(xnodes=3,ynodes=9,znodes=18,spacing=2,corner=ar([0,0,80]),com_v=ar([0,0,-10]),ang_mom=ar([8000,-4500,2000]),vlc=1,col_N=10,floor_tol=0,g=-6,ground_E_tol=40,fric=1,node_mass=.1)
block.animate_block_energy(save_name="block_collide_e_cons2.mp4",size=8,dt=.01,tf=200,t_plot_range=50,steps_per_frame=4,camera_rot=2)
#%%
