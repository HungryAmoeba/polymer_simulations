#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helper.h"

System* create_system(int N, int board_rows, int board_cols, double l_p, double b, int geometry) {
    System* sys = malloc(sizeof(System));
    sys->n_monomers = N;
    sys->board_cols = board_cols;
    sys->board_rows = board_rows;
    sys->b = b;
    sys->l_p = l_p;
    sys->curr_energy = 0;
    sys->geometry = geometry;

    sys->monomer_locations = calloc(2 * N , sizeof(int));
    sys->lattice = calloc(board_cols * board_rows, sizeof(int));
    sys->derivatives = calloc(2 * N , sizeof(double));
    sys->diff_derivatives = calloc((N-1) , sizeof(double)); // stores the square length of adjacent derivative terms

    return sys;
}

void destroy_system(System* sys) {
    free(sys->monomer_locations);
    free(sys->lattice);
    free(sys->derivatives);
    free(sys->diff_derivatives);
    free(sys);
}

System* deep_copy_system(System *sys) {
    System* sys_copy = create_system(sys->n_monomers, sys->board_rows, sys->board_cols, sys->l_p, sys->b, sys->geometry);
    for (int i = 0; i < 2 * sys->n_monomers; i++) {
        sys_copy->monomer_locations[i] = sys->monomer_locations[i];
        sys_copy->derivatives[i] = sys->derivatives[i];
    }
    for (int i = 0; i < sys->board_cols*sys->board_rows; i++) {
        sys_copy->lattice[i] = sys->lattice[i];
    }
    for (int i = 0; i < sys->n_monomers-1; i++) {
        sys_copy->diff_derivatives[i] = sys->diff_derivatives[i];
    }
    sys_copy->curr_energy = sys->curr_energy;
    return sys_copy;
}

void initialize_starting_state(System* sys) {
    // the starting state cannot be too low energy
    // populate starting monomer locations
    sys->monomer_locations[0] = sys->board_rows / 2;
    sys->monomer_locations[sys->n_monomers] = sys->board_cols / 2;
    int monomer_ind = 0;
    int segment_length = 1;
    int segment_counter = 0;
    int segment_leg = 0;
    // draw a spiral originating from the center of the board
    for (int i = 1; i < sys->n_monomers; i++) {
        int direction_x = (int) -cos(M_PI/2 * (double) segment_counter);
        int direction_y = (int) sin(M_PI/2 * (double) segment_counter); 
        sys->monomer_locations[i] = sys->monomer_locations[i-1] + 2 * direction_x;
        sys->monomer_locations[i + sys->n_monomers] = sys->monomer_locations[ i -1 + sys->n_monomers] + 2 * direction_y;
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

    // update the lattice to reflect the monomer locations
    for (int i = 0; i < sys->n_monomers; i++) {
        sys->lattice[sys->monomer_locations[i]* sys->board_rows + sys->monomer_locations[i + sys->n_monomers]] = 1;
    }

    // print out the lattice
    print_board(sys);

    // calculate the initial derivatives
    for (int i = 0; i < sys->n_monomers; i++) {
        calc_d_ind(sys, i);
    }

    // calculate the difference in derivatives
    for (int i = 0; i < sys->n_monomers - 1; i++) {
        double i_ind = sys->derivatives[i + 1] - sys->derivatives[i];
        double j_ind = sys->derivatives[i + 1 + sys->n_monomers] - sys->derivatives[i + sys->n_monomers];
        sys->diff_derivatives[i] = i_ind*i_ind + j_ind*j_ind;
    }

    // calculate the initial starting energy
    for (int i = 0; i < sys->n_monomers-1; i++) {
        sys->curr_energy += - sys->l_p/(2 * sys->b) * sys->diff_derivatives[i];
    }
}

void print_board(System *sys) {
    for (int i = 0; i < sys->board_cols * sys->board_rows; i++) {
        printf("%d, ", sys->lattice[i]);
        if ((i+1) % sys->board_cols == 0) {
            printf("\n");
        }
    }   
}

void calc_d_ind_cylinder(System *sys, int index) {
    // TODO
}
void calc_d_ind(System *sys, int index) {
    if (index < 0 || index > sys->n_monomers) {
        printf("error, out of bounds\n");
    }
    double i_ind, j_ind, norm;

    if (index == 0) {
        // normalize R2 - R1
        i_ind = (double) (sys->monomer_locations[1] - sys->monomer_locations[0]);
        j_ind = (double) (sys->monomer_locations[sys->n_monomers + 1] - sys->monomer_locations[sys->n_monomers]);
    } else if (index == (sys->n_monomers - 1)) {
        // normalize RN - R(N-1)
        i_ind = (double) (sys->monomer_locations[index] - sys->monomer_locations[index - 1]);
        j_ind = (double) (sys->monomer_locations[index + sys->n_monomers] - sys->monomer_locations[index -1 + sys->n_monomers]);
    } else {
        // calculate left hand side terms
        double i_ind_left = (double) (sys->monomer_locations[index] - sys->monomer_locations[index -1]);
        double j_ind_left = (double) (sys->monomer_locations[index + sys->n_monomers] - sys->monomer_locations[index + sys->n_monomers - 1] );
        double left_norm = sqrt(i_ind_left*i_ind_left + j_ind_left*j_ind_left); 
        double i_ind_right = (double) (sys->monomer_locations[index + 1] - sys->monomer_locations[index]);
        double j_ind_right = (double) (sys->monomer_locations[index + sys->n_monomers + 1] - sys->monomer_locations[index + sys->n_monomers] );
        double right_norm = sqrt(i_ind_right*i_ind_right + j_ind_right*j_ind_right);
        i_ind = (double) i_ind_left/left_norm + i_ind_right/right_norm;
        j_ind = (double) j_ind_left/left_norm + j_ind_right/right_norm;
    }
    // complete the calculation and make the assignment
    norm = sqrt(i_ind*i_ind + j_ind*j_ind);
    sys->derivatives[index] = i_ind/norm;
    sys->derivatives[index + sys->n_monomers] = j_ind/norm;
}

void shuffle(int *array, size_t n)
{
    if (n > 1) 
    {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

int *get_move(int move_num) {
  int *move_arr = malloc(2 * sizeof(int));
  switch(move_num) {
    case 0:
      move_arr[0] = -1; move_arr[1] = -1;
      break;
    case 1:
      move_arr[0] = -1; move_arr[1] = 0;
      break;
    case 2:
      move_arr[0] = -1; move_arr[1] = 1;
      break;
    case 3:
      move_arr[0] = 0; move_arr[1] = -1;
      break;
    case 4:
      move_arr[0] = 0; move_arr[1] = 1;
      break;
    case 5:
      move_arr[0] = 1; move_arr[1] = -1;
      break;
    case 6:
      move_arr[0] = 1; move_arr[1] = 0;
      break;
    case 7:
      move_arr[0] = 1; move_arr[1] = 1;
      break;
  }
  return move_arr;
}

int check_legal_cylinder(int ind, System* sys, int* mv_arr) {
    int pi = sys->monomer_locations[ind] + mv_arr[0];
    int pj = sys->monomer_locations[ind + sys->n_monomers] + mv_arr[1];

    pi = (pi + sys->n_monomers)%sys->n_monomers;

    if (pj < 0 || pj >= sys->board_cols) {
        return 0;
    }
    if (sys->lattice[pi*sys->board_cols + pj] == 1) {
        return 0;
    }
    // check for excluded volume in neighbors and diagonals
    // don't check move array index
    // only check if in bounds

    // i index search = 
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            if (i != -mv_arr[0] || j != -mv_arr[1]) {
                // check that the j-index is not out of bounds
                if (j + pj >= 0 && j + pj < sys->board_cols) {
                    int row_num = (pi + i + sys->n_monomers)%sys->n_monomers;
                    int col_num = pj + j;
                    if (sys->lattice[row_num* sys->board_cols + col_num] == 1) {
                        return 0;
                    }
                }
            }
        }
    }

    // check the distance works
    if (ind != 0) {
        int dist_x = pi - sys->monomer_locations[ind - 1]; 

        int dist_y_no_wrap = pj - sys->monomer_locations[ind - 1 + sys->n_monomers];
        int dist_y_wrap = sys->n_monomers - dist_y_no_wrap;
        int dist_y = min(abs(dist_y_no_wrap), abs(dist_y_wrap));
    
        if (dist_x*dist_x + dist_y*dist_y > 13) {
            return 0;
        }
    }
    if (ind != sys->n_monomers -1) {
        // compare to next monomer
        int dist_x = pi - sys->monomer_locations[ind + 1];
        int dist_y_no_wrap = pj - sys->monomer_locations[ind + 1 + sys->n_monomers];
        int dist_y_wrap = sys->n_monomers - dist_y_no_wrap;
        int dist_y = min(abs(dist_y_no_wrap), abs(dist_y_wrap));
        if (dist_x*dist_x + dist_y*dist_y > 13) {
            return 0;
        } 
    }
    return 1;
}

int check_legal(int ind, System* sys, int* mv_arr) {
    // call other function if it's a cylinder
    if (sys->geometry == 2) {
        int result = check_legal_cylinder(ind, sys, mv_arr);
        return result;
    }
    int pi = sys->monomer_locations[ind] + mv_arr[0];
    int pj = sys->monomer_locations[ind + sys->n_monomers] + mv_arr[1];

    if (pi < 0 || pi >= sys->board_rows) {
        return 0;
    }
    if (pj < 0 || pj >= sys->board_cols) {
        return 0;
    }
    if (sys->lattice[pi*sys->board_cols + pj] == 1) {
        return 0;
    }
    // check for excluded volume
    // i am so sorry for this horrible code style
    for (int i = -1; i < 2; i ++) {
        for (int j = -1; j < 2; j ++) {
            if (i != -mv_arr[0] || j != -mv_arr[1]) {
                if (i + pi >= 0 && i + pi < sys->board_rows) {
                    if (j + pj >= 0 && j + pj < sys->board_cols) {
                        if (sys->lattice[(pi + i)* sys->board_cols + pj + j] == 1) {
                            return 0;
                        }
                    }
                }
            }
        }
    }

    if (ind != 0) {
        // compare to previous monomer
        int dist_x = pi - sys->monomer_locations[ind - 1]; 
        int dist_y = pj - sys->monomer_locations[ind - 1 + sys->n_monomers];
        //printf("%d is pi, %d is prev, %d is sq\n", pi, states[index-1][0], dist_x*dist_x);
        //printf("%d is the distance on left \n", dist_x*dist_x + dist_y*dist_y);
        if (dist_x*dist_x + dist_y*dist_y > 13) {
            //printf("rejected \n \n");
            return 0;
        }
    }
    if (ind != sys->n_monomers -1) {
        // compare to next monomer
        int dist_x = pi - sys->monomer_locations[ind + 1];
        int dist_y = pj - sys->monomer_locations[ind + 1 + sys->n_monomers];
        if (dist_x*dist_x + dist_y*dist_y > 13) {
            return 0;
        } 
    }
    return 1;
}

double update_system(System *sys, int ind) {
    // assume that system has been updated at a certain location. 
    // update the derivatives, dd, and energy
    int left_ind = ind - 1;
    int right_ind = ind + 1;
    if (left_ind < 0) {
        left_ind = 0;
    }
    if (right_ind >= sys->n_monomers) {
        right_ind = ind;
        // right_ind = sys->n_monomers - 1;
    }
    for (int i = left_ind; i < right_ind + 1; i++) {
        // update the derivatives
        calc_d_ind(sys, i);
    }

    int dd_left_ind = left_ind- 1;
    int dd_right_ind = right_ind + 1;
    if (dd_left_ind < 0) {
        dd_left_ind = 0;
    }
    if (dd_right_ind >= sys->n_monomers - 1) {
        // the index should go up to N - 2
        dd_right_ind = sys->n_monomers - 2;
    }
    double difference_ce = 0; // new energy - old energy
    for (int i = dd_left_ind; i < dd_right_ind + 1; i ++) {
        double i_ind = sys->derivatives[i + 1] - sys->derivatives[i];
        double j_ind = sys->derivatives[i + 1 + sys->n_monomers] - sys->derivatives[i + sys->n_monomers];
        difference_ce += - sys->l_p/(2 * sys->b) * (i_ind*i_ind + j_ind*j_ind - sys->diff_derivatives[i]);
        sys->diff_derivatives[i] = (i_ind*i_ind + j_ind*j_ind);
    }
    
    // calculate the difference in current energy and return it 
    sys->curr_energy = sys->curr_energy + difference_ce;
    return sys->curr_energy;
}

// Copies from sys_1 to sys_2
void sys_copy_over(System* sys_1, System* sys_2, int ind) {
    // update the monomer location
    sys_2->monomer_locations[ind] = sys_1->monomer_locations[ind];
    sys_2->monomer_locations[ind + sys_2->n_monomers] = sys_1->monomer_locations[ind + sys_1->n_monomers];
    
    // update the derivatives
    int left_ind = ind - 1;
    int right_ind = ind + 1;
    if (left_ind < 0) {
        left_ind = 0;
    }
    if (right_ind >= sys_1->n_monomers) {
        right_ind = sys_1->n_monomers - 1;
    }
    for (int i = left_ind; i < right_ind + 1; i++) {
        // update the derivatives
        sys_2->derivatives[i] = sys_1->derivatives[i];
        sys_2->derivatives[i + sys_2->n_monomers] = sys_1->derivatives[i + sys_1->n_monomers];
    }

    // update dd_array
    int dd_left_ind = left_ind- 1;
    int dd_right_ind = right_ind + 1;
    if (dd_left_ind < 0) {
        dd_left_ind = 0;
    }
    if (dd_right_ind >= sys_1->n_monomers - 1) {
        // the index should go up to N - 2
        dd_right_ind = sys_1->n_monomers - 2;
    }

    for (int i = dd_left_ind; i < dd_right_ind + 1; i ++) {
        sys_2->diff_derivatives[i] = sys_1->diff_derivatives[i];
    }

    // update current energy
    sys_2->curr_energy = sys_1->curr_energy;
}

int min(int a, int b) {
    if (a < b) {
        return a;
    }
    return b;
}
