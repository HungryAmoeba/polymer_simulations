import matplotlib.pyplot as plt
import numpy as np
import re
import random
from tqdm import tqdm
import os
from mpl_toolkits.mplot3d import Axes3D
from animation import rotanimate

N = 60
rows = 30
cols = 30
locs = np.zeros((N,2))
ind = 0
arr_dict = {}
mcstep = 0
step_list = list()

with open('simulation_results/results.txt') as f:
    for line in f.readlines():
        line = line.strip()
        if 'Step is' in line:
            mcstep = int(re.findall(r'\d+', line)[0])
            arr_dict[mcstep] = np.zeros((N,2))
            ind = 0
            step_list.append(mcstep)
        else:
            arr_dict[mcstep][ind] = [int(s) for s in re.findall(r'\b\d+\b', line)]
            ind = ind + 1
print(step_list)
for s in step_list:
    fig = plt.figure(figsize=(15,15))
    r = cols/(2 * np.pi)
    ax = fig.add_subplot(111, projection= '3d')
    z = np.linspace(0, rows, 600)
    x = np.linspace(- round( r), round(r), 600)
    Xc, Zc = np.meshgrid(x,z)
    Yc = np.sqrt((r**2 - Xc**2))

    rstride = 20
    cstride = 10
    ax.plot_surface(Xc,Yc,Zc,alpha =0.2, rstride=rstride, cstride=cstride)
    ax.plot_surface(Xc,-Yc,Zc, alpha=0.2, rstride=rstride, cstride=cstride)

    # map the monomer_states to points on the cylinder
    thetas = 2 * np.pi/ cols * arr_dict[s][:,1]
    x_data = np.cos(thetas) * r
    y_data = np.sin(thetas) * r
    z_data = arr_dict[s][:,0]

    ax.plot(x_data, y_data,z_data)
    s_ax = ax
    plt.axis('off')
    angles = np.linspace(0,360,41)

    rotanimate(ax, angles, f'plotted_states/polymer_step_{s}.gif', delay=10)

