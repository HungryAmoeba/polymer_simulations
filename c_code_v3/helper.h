#ifndef __HELPER_H__
#define __HELPER_H__

#include <stdio.h>

typedef struct System
{
    int geometry;               // 0 for planar or 1 for cylindrical
    int n_monomers;             // N
    int board_rows;             // The number of rows
    int board_cols;             // number of columns
    double l_p;                 // persistence length
    double b;                   // kuhn length
    int* monomer_locations;     // Array storing the coordinates of the monomers
    int* lattice;               // Lattice showing open and occupied squares
    double* derivatives;        // First derivatives of the monomer array
    double* diff_derivatives;   // Second derivatives of monomer array
    double curr_energy;         // Current bending energy
} System;

System* create_system(int N, int board_rows, int board_cols, double l_p, double b, int geometry);

void print_board(System *sys);

void destroy_system(System* sys);

System* deep_copy_system(System *sys);

void initialize_starting_state(System* sys);

void calc_d_ind_cylinder(System *sys, int index);

void calc_d_ind(System *sys, int index);

/**
 * Random shuffling of an array
 * following Fisher-Yates
 *
 * @param array an integer array
 * @param n the size of that array
 * @return nothing. Shuffles in places
 */
void shuffle(int *array, size_t n);

/**
 * Gets one move that's nearest neighbors and diagonals
 *
 * @param move_num an integer 0 to 7
 * @return a 2D array with i,j coordinates
 */
int *get_move(int move_num);

int check_legal_cylinder(int ind, System* sys, int* mv_arr);

/**
 * Check that a proposed move is legal, meaning that it is
 * not out of bounds, or does not violate excluded volume or
 * exceed maximum allowed bond length
 *
 * @param ind the index of the monomer which is proposed to be moved
 * @param sys The system containing the lattice and monomer locations
 * @param mv_arr An array with the index of the proposed move.
 * @return 1 if the move is legal, 0 otherwise
 */
int check_legal(int ind, System* sys, int* mv_arr);

double update_system(System *sys, int ind);

void sys_copy_over(System* sys_1, System* sys_2, int ind);

int min(int a, int b);

#endif
