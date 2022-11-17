import numpy as np
import matplotlib.pyplot as plt

import sys

from cylindrical_polymer import polymer_chain

argc = len(sys.argv)
if argc != 5:
    print(f"Usage: python cylindrical_main.py --N --persistence_length --C --h")
    print("N is equivalent to L in the paper")
    print(quit)
    quit()

N = int(sys.argv[1]); l_p = int(sys.argv[2]); C = int(sys.argv[3]); h = int(sys.argv[4])

# initialize the grid

#kuhn length
b = 2.8

N_steps = 10000
polymer = polymer_chain(N, (C,h), l_p)
polymer.populate_grid()
polymer.visualize_grid()

#polymer.calculate_bending()
#polymer.metropolis(N_steps = N_steps + 1, savefig = 1000)
print("success")
