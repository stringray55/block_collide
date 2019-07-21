import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import time

def anim_init_3d():
    ax = plt.axes(projection='3d')
    object = ax.scatter([],[],[])
    return object

def do_60fps(save_name, fig, anim_init, frame, tot_frames):
    print('animating')
    start=time.time()
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1800)
    anim = animation.FuncAnimation(fig, frame, init_func=anim_init, frames=tot_frames, interval=20, blit=False)
    anim.save(save_name, writer=writer)
    end=time.time()
    print('time elapsed: '+str(end-start)+' seconds')

def do_30fps(save_name, fig, anim_init, frame, tot_frames):
    print('animating')
    start=time.time()
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    anim = animation.FuncAnimation(fig, frame, init_func=anim_init, frames=tot_frames, interval=20, blit=False)
    anim.save(save_name, writer=writer)
    end=time.time()
    print('time elapsed: '+str(end-start)+' seconds')
