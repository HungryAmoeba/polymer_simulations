import matplotlib.pyplot as plt
import numpy as np
import re

N = 50
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
    fig, ax = plt.subplots()
    plt.xlim([0,N]); plt.ylim([0,N])
    plt.plot(arr_dict[s][:,0], arr_dict[s][:,1], marker='*')
    plt.savefig(f'plotted_states/sim{s}.jpg')