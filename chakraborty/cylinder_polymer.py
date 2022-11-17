import numpy as np
import random
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

from polymer import permute_moves
from polymer import normalize

class cylindrical_polymer:
    def __init__(self, N, grid_dimension, l_p):
        self.N = N
        self.grid_dimension = grid_dimension
        self.pad_width = 4
        self.grid = np.concatenate(np.concatenate(np.ones((self.pad_width, grid_dimension[0])) *2, \
            np.zeros(grid_dimension)), np.ones((self.pad_width, grid_dimension[0])) * 2)
        self.l_p = l_p
        self.b = 2.8
        self.derivatives = np.zeros((N,2))
        self.bending_energy = 0
        self.allowed_bonds = permute_moves(([2,0],  

