import main_block_collide as mbc
import matplotlib.pyplot as plt
import numpy as np
ar=np.array

block=mbc.Block(xnodes=3,ynodes=9,znodes=18,spacing=2,corner=ar([0,0,80]),com_v=ar([0,0,-10]),ang_mom=ar([6000,-4500,2000]),vlc=.8,col_N=10,floor_tol=0,g=-6,ground_E_tol=10,fric=0.5,node_mass=.1)
block.animate_block_energy(save_name="block_collide2246.mp4",size=8,dt=.01,tf=30,steps_per_frame=4,camera_rot=2)
#%%
