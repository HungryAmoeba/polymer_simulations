import numpy as np
import random
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

def permute_moves(tuples):
    moves = list()
    for bond in tuples:
        for sign_1 in (1, -1):
            # assume that if there is a zero, it comes in the second spot:
            if bond[1] == 0:
                sign_2_arr = [1]
            else:
                sign_2_arr = [1,-1]
            for sign_2 in sign_2_arr:
                moves.append(bond * np.array([sign_1, sign_2]))
                if bond[0] != bond[1]:
                    moves.append(np.flip(bond * np.array([sign_1, sign_2])))
    return moves

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v/ norm

class polymer_chain:
    def __init__(self, N, grid_dimension, l_p):
        self.N = N
        self.grid_dimension = grid_dimension
        self.pad_width = 4
        # add padding to the grid, with out of bound values being assigned 2
        self.grid = np.pad(np.zeros((grid_dimension, grid_dimension), dtype = np.intc), pad_width = self.pad_width, mode = 'constant', constant_values = 2) 
        self.l_p = l_p
        self.b = 2.8 # average bond length
        self.derivatives = np.zeros((N,2)) # u in the paper
        self.bending_energy = 0
        self.allowed_bonds = permute_moves(([2,0], [2,1], [2,2], [3,0], [3,1], [3,2]))
        self.monomer_states = np.zeros((N,2), dtype = np.intc) # true positions. doesn't depend on boundary
        self.beta = 1
        self.energy_history = np.zeros(1)
    
    def populate_grid(self, population_attempts = 0):
        if population_attempts >= 50:
            print(f"could not populate board in {population_attempts} attempts")
            raise Exception("Could not populate board. Try again, or reduce N")
        # randomly choose a starting point, and mark it in the grid
        grid_offset = np.array([self.pad_width, self.pad_width])
        self.monomer_states[0] = np.random.randint(0, high = self.grid_dimension, size = (1,2))
        self.grid[tuple(self.monomer_states[0] + grid_offset)] = 1
        
        # loop through the rest of the monomer states
        for index in range(1,self.N):
            point_unassigned = 1; tries = 0
            move_order = np.random.permutation(self.allowed_bonds) 
            while point_unassigned and tries < 8:
                # propose a random move
                move = move_order[tries]
                #ensure that this move is not out of bounds
                proposed_grid_location = self.monomer_states[index - 1] + move + grid_offset 
                legal_move = 1
                if self.grid[tuple(proposed_grid_location)] != 0:
                    legal_move = 0
                for x_offset in range(-1, 2):
                    for y_offset in range(-1,2):
                        if self.grid[tuple(proposed_grid_location + [x_offset, y_offset])] == 1:
                            legal_move = 0
            # proposed move is ok if we get to here
                if legal_move:
                    self.grid[tuple(proposed_grid_location)] = 1
                    self.monomer_states[index] = self.monomer_states[index - 1] + move
                    point_unassigned = 0
                tries = tries + 1
            if point_unassigned:
                # populating the board failed
                population_attempts = population_attempts + 1
                #clear everything and try again TODO
                for monomer in self.monomer_states:
                    self.grid[tuple(monomer + grid_offset)] = 0
                self.monomer_states = np.zeros((self.N,2), dtype = np.intc) 
                self.populate_grid(population_attempts = population_attempts)

                 

    def visualize_grid(self, n_steps, directory = 'None'):
        plt.figure()
        plt.plot(self.monomer_states[:,0], self.monomer_states[:,1], linewidth = 1, marker = "*")
        plt.xticks([],[])
        plt.yticks([],[])
        plt.title(f'N = {self.N}, l_p = {self.l_p}, W = {self.grid_dimension},E = {round(self.bending_energy)}, step = {n_steps} ')
        if directory == 'None':
            plt.savefig(f'{self.N}_{self.l_p}_{self.grid_dimension}_{n_steps}.png')
        else:
            plt.savefig(directory + f'/states/{self.N}_{self.l_p}_{self.grid_dimension}_{n_steps}.png')

    def visualize_derivatives(self, n_steps, directory = 'None'):
        plt.figure()
        # add the derivatives together end to end
        connected_derivatives = np.zeros(self.derivatives.shape)
        connected_derivatives[0] = self.derivatives[0]
        for i in range(1,self.N):
            connected_derivatives[i] = connected_derivatives[i-1] + self.derivatives[i]
        plt.plot(connected_derivatives[:,0], connected_derivatives[:,1], linewidth = 1, marker = "*")
        plt.xticks([],[]); plt.yticks([],[])
        plt.title(f'Derivatives, step = {n_steps}')
        plt.savefig(directory + f'/derivatives/derivatives_{n_steps}.png')

    def print_grid(self):
        print(self.grid)

    def calculate_bending_energy(self, derivative_array):
        # assumes that the derivative array is of length N
        #end_zero_derivative = np.append(derivative_array, np.zeros((1,2)), axis = 0)
        #start_zero_derivative = np.insert(derivative_array, 0, np.zeros((1,2)), axis=0)
        #norms = np.linalg.norm(end_zero_derivative - start_zero_derivative, axis = 1)
        #relevant_diffs = norms[0:self.N - 1]
        #bending_energy = self.l_p / (2 * self.b) * sum( np.square(relevant_diffs))
        sum = 0
        for ind in range(0, self.N-1):
            sum += np.linalg.norm( derivative_array[ind + 1] - derivative_array[ind] ) **2
        bending_energy = sum * self.l_p / (2 * self.b)
        return bending_energy

    def save_states(self, iteration, directory = 'None'):
        # save the monomer_states
        np.save(directory + f'/states/_{iteration}_monomer_states.npy', self.monomer_states)
        # save the derivatives 
        np.save(directory + f'/derivatives/{iteration}_derivatives.npy', self.derivatives)

    def calculate_bending(self):
        for ind in range(1, self.N - 1):
            self.derivatives[ind] = normalize( normalize(self.monomer_states[ind] - self.monomer_states[ind - 1]) + normalize(self.monomer_states[ind + 1] - self.monomer_states[ind]))
        self.derivatives[0] = normalize( self.monomer_states[1] - self.monomer_states[0] )
        self.derivatives[self.N-1] = normalize( self.monomer_states[self.N - 1] - self.monomer_states[self.N - 2] )
        self.bending_energy = self.calculate_bending_energy(self.derivatives)

    def metropolis(self, N_steps = 10000, savefig = 0):
        self.energy_history = np.zeros(N_steps // 10 + 1)
        if savefig != 0:
            dir_name = f'simulation_results/{self.N}_{self.l_p}_{self.grid_dimension}_simulation'
            os.mkdir(dir_name)
            os.mkdir(dir_name + '/states')
            os.mkdir(dir_name + '/derivatives')
        allowed_moves = permute_moves(([1,1], [1,0]))
        grid_offset = np.array([self.pad_width, self.pad_width])
        for iteration in tqdm(range(N_steps)):
            if iteration % savefig == 0:
                self.visualize_grid(iteration, directory = dir_name)
                self.save_states(iteration, directory = dir_name)
                self.visualize_derivatives(iteration, directory = dir_name)
            # at each iteration, visit every point
            if iteration % 10 == 0:
                self.energy_history[iteration//10] = self.bending_energy
            point_order_list = np.random.randint(0, high = self.N, size = self.N)
            #__import__('pdb').set_trace()
            for point in point_order_list:
                point_unassigned = 1
                tries = 0
                move_order_choice = np.random.permutation(8)
                while point_unassigned and tries < 8:
                    # propose a random move
                    move = allowed_moves[move_order_choice[tries]]
                    tries = tries + 1
                    #ensure that this move is not out of bounds
                    proposed_grid_location = self.monomer_states[point] + move + grid_offset
                    legal_move = 1
                    if self.grid[tuple(proposed_grid_location)] != 0:
                        legal_move = 0
                    # check that there are no points valued 1 in a 2x2 
                    
                    for x_offset in range(-1, 2):
                        for y_offset in range(-1,2):
                            if not np.array_equal([x_offset, y_offset], -move):
                                if self.grid[tuple(proposed_grid_location + [x_offset, y_offset])] == 1:
                                    legal_move = 0
                    # check that the proposed grid location is not too far away:
                    # check new distance to prior equation
                    if point != 0:
                        dist = self.monomer_states[point] + move - self.monomer_states[point - 1]
                        if np.linalg.norm(dist) >= 4:
                            legal_move = 0
                    if point != self.N - 1:
                        dist = self.monomer_states[point] + move - self.monomer_states[point + 1]
                        if np.linalg.norm(dist) >= 4:
                            legal_move = 0
                    
                    if legal_move: 
                        # proposed move is legal if we get to here
                        point_unassigned = 0
                        # calculate energy of proposed move
                        new_monomer_states = np.copy(self.monomer_states)
                        new_monomer_states[point] = new_monomer_states[point] + move
                        new_monomer_derivatives = np.copy(self.derivatives)
                        # only recalculate the relevant derivatives
                        for ind in range(np.max((0,point - 2)), np.min((self.N-1, point + 2))):
                            new_monomer_derivatives[ind] = normalize( normalize(new_monomer_states[ind] - new_monomer_states[ind - 1]) + normalize(new_monomer_states[ind + 1] - new_monomer_states[ind])) 
                        new_energy = self.calculate_bending_energy(new_monomer_derivatives)

                        energy_diff = new_energy - self.bending_energy
                        if np.exp(-1 * energy_diff) > np.random.uniform():
                            # reset the previous board state
                            self.grid[tuple(self.monomer_states[point] + grid_offset)] = 0
                            # accept the move 
                            self.grid[tuple(proposed_grid_location)] = 1
                            self.monomer_states[point] = self.monomer_states[point] + move
                            self.derivatives = new_monomer_derivatives
                            self.bending_energy = new_energy


        fig, ax = plt.subplots()
        ax.plot(self.energy_history)
        ax.set_title('Energy history')
        fig.show()
    
