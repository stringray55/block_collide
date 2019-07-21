import main_block_collide as mbc
import matplotlib.pyplot as plt
import numpy as np
ar=np.array

block=mbc.Block(xnodes=3,ynodes=9,znodes=18,spacing=2,corner=ar([0,0,60]),com_v=ar([0,0,-10]),ang_mom=ar([5000,-4500,2000]),Elc=1,col_N=10,floor_tol=0,g=-6,node_mass=.1)
block.animate_block2(save_name="block_collide2239.mp4",size=8,dt=.01,tf=100,steps_per_frame=4,camera_rot=3)
#%%
