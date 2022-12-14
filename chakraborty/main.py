import numpy as np
import matplotlib.pyplot as plt

import sys

from polymer import polymer_chain
import cProfile
import re

argc = len(sys.argv)
if argc != 4:
    print(f"Usage: python main.py --N --persistence_length --box_length")
    print("N is equivalent to L in the paper")
    print(quit)
    quit()

N = int(sys.argv[1]); l_p = int(sys.argv[2]); W = int(sys.argv[3]);
print(f"N is {N}, l_p is {l_p}, W is {W}")
# initialize the grid

#kuhn length
b = 2.8

box_length = int(np.ceil(b * W))
N_steps = int(10e3)
polymer = polymer_chain(N, box_length, l_p)
polymer.populate_grid()
polymer.calculate_bending()
polymer.metropolis(N_steps = N_steps + 1, savefig = int(10e2))
print("success")
