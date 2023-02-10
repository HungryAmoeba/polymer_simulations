#ifndef __HELPER_H__
#define __HELPER_H__

#include <stdio.h>
#include "vectors.h"

void print_array(int **arr, int max_i, int max_j);

void initialize_starting_state(int_array* board, ivec* monomer_array, int N);

int check_legal(int index, int_array* monomer_locations, int_array *board_states, int *move_arr);

void calc_derivatives(ivec* monomer_locations, vec_array* derivative_array);

void calc_dd(vec_array* derivative_array, vec_array* dd_array, double* dd_energy_array);

#endif
