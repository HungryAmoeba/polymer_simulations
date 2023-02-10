#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "helper.h"

void print_array(int **arr, int max_i, int max_j) {
    for (int i = 0; i < max_i; i++) {
        for (int j = 0; j < max_j; j++) {
            printf("%d, ", arr[i][j]);
        }
        printf("\n");
    }
}

// takes in the board and the number of monomers
void initialize_starting_state(int_array* board, ivec* monomer_locations, int N) {
  // the starting state cannot be too low energy
  // populate starting monomer locations
  monomer_locations[0] = create_ivec(board->rows / 2, board->columns / 2);
  int monomer_ind = 0;
  int segment_length = 1;
  int segment_counter = 0;
  int segment_leg = 0;
  // draw a spiral originating from the center of the board
  for (int i = 1; i < N; i++) {
      ivec direction = create_ivec((int) - cos(M_PI/2 * (double) segment_counter), (int) sin(M_PI/2 * (double) segment_counter));
      monomer_locations[i] = add_ivec(monomer_locations[i - 1], scalar_multiplication(2, direction));
      monomer_ind++;
      if (monomer_ind == segment_length) {
          segment_counter++;
          monomer_ind = 0;
          segment_leg++; 
      }
      if (segment_leg == 2) {
          segment_leg = 0;
          segment_length++;
      }
  }

  for (int i = 0; i < N; i++) {
    int i_ind = monomer_locations[i].i_ind;
    int j_ind = monomer_locations[i].j_ind;
    board->array[i_ind][j_ind] = 1;
  }
}

int check_legal(int index, int_array* monomer_locations, int_array *board_states, int *move_arr) {
  // board states carries the lattice, monomer_locations is N x 2 array
  int N = monomer_locations->rows;
  // check that the proposed move is not out of bounds
  int pi = board_states->array[index][0] + move_arr[0];
  int pj = board_states->array[index][1] + move_arr[1];
  if (pi < 0 || pi >= board_states->rows) {
    return 0;
  }
  if (pj < 0 || pj >= board_states->columns) {
    return 0;
  }
  if (board_states->array[pi][pj] == 1) {
    return 0;
  }
  // move is not out of bounds, and proposed square is not occupied

  // check for excluded volume
  // i am so sorry for this horrible code style
  for (int i = -1; i < 2; i ++) {
    for (int j = -1; j < 2; j ++) {
      if (i != -move_arr[0] || j != -move_arr[1]) {
        if (i + pi >= 0 && i + pi < board_states->rows) {
          if (j + pj >= 0 && j + pj < board_states->columns) {
            if (board_states->array[pi +i][pj + j] == 1) {
              return 0;
            }
          }
        }
        
      }
    }
  }

  // check length does not exceed 3
  if (index != 0) {
    // compare to previous monomer
    int dist_x = pi - monomer_locations->array[index-1][0];
    int dist_y = pj - monomer_locations->array[index-1][1];
    //printf("%d is pi, %d is prev, %d is sq\n", pi, states[index-1][0], dist_x*dist_x);
    //printf("%d is the distance on left \n", dist_x*dist_x + dist_y*dist_y);
    if (dist_x*dist_x + dist_y*dist_y > 13) {
      //printf("rejected \n \n");
      return 0;
    }
  }
  if (index != N-1) {
    // compare to next monomer
    int dist_x = pi - monomer_locations->array[index+1][0];
    int dist_y = pj - monomer_locations->array[index+1][1];
    if (dist_x*dist_x + dist_y*dist_y > 13) {
      return 0;
    } 
  }
  return 1;
}

void calc_derivatives(ivec* monomer_locations, vec_array* derivative_array) {
  int N = derivative_array->n_elements;
  ivec derivative_0 = subtract_ivec(monomer_locations[1], monomer_locations[0]);
  derivative_array->array[0] = normalize_vec(ivec_to_vec(derivative_0));
  ivec derivative_last = subtract_ivec(monomer_locations[N-1], monomer_locations[N-2]);
  derivative_array->array[N-1] = normalize_vec(ivec_to_vec(derivative_last));
  for (int i = 1; i < N-1; i ++) {
    vec derivative_left = normalize_vec(ivec_to_vec(subtract_ivec(monomer_locations[i], monomer_locations[i-1])));
    vec derivative_right = normalize_vec(ivec_to_vec(subtract_ivec(monomer_locations[i+1], monomer_locations[i])));
    derivative_array->array[i] = normalize_vec(add_vec(derivative_left, derivative_right));
  }
}

void calc_dd(vec_array* derivative_array, vec_array* dd_array, double* dd_energy_array) {
  int N = derivative_array->n_elements;
  for (int i = 0; i < N-1; i++) {
    printf("i =%d in calc_dd\n", i);
    vec dd_vec = subtract_vec(derivative_array->array[i+1], derivative_array->array[i]);
    printf("Finished calculating dd = %f, %f \n", dd_vec.i_ind, dd_vec.j_ind);
    // dd_array->array[i].i_ind = dd_vec.i_ind;
    // dd_array->array[i].j_ind = dd_vec.j_ind;
    dd_array->array[i] = dd_vec;
    printf("dd_array is %f, %f \n", dd_array->array[i].i_ind, dd_array->array[i].j_ind);
    printf("Finished assigning dd_array\n");
    dd_energy_array[i] = square_vec(dd_array->array[i]);
  }
  printf("exiting loop \n");
}
