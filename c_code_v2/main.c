#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helper.h"
#include "vectors.h"

int main(int argc, char **argv) {

    // expect to get Main rows, columns, N
    int rows, cols, N;
    if (argc != 4) {
        printf("Usage: ./Main {lattice rows} {lattice columns} {number of monomers}\n");
        exit(0);
    }
    else {
        rows = atoi(argv[1]);     // lattice board rows
        cols = atoi(argv[2]);     // lattice board columns
        N = atoi(argv[3]);        // number of monomers
    }
    double l_p = 2.0;                   // persistence length
    double b = 2.8;                      // average bond length

    int_array* board = create_int_array(rows,cols);
    // int_array* monomer_locations = create_int_array(N, 2);
    ivec* monomer_locations = malloc(N * sizeof(ivec));

    initialize_starting_state(board, monomer_locations, N);

    vec_array* derivative_array = create_vec_array(N);
    double* dd_energy_array = malloc(N-1 * sizeof(double)); // this is NOT multiplied by l_p/2b yet
    vec_array* dd_array = create_vec_array(N-1);

    calc_derivatives(monomer_locations, derivative_array);
    // for (int i = 0; i < N; i++) {
    //     printf("Derivative %d is %f, %f \n", i, derivative_array->array[i].i_ind, derivative_array->array[i].j_ind);
    // }
    calc_dd(derivative_array, dd_array, dd_energy_array);

    print_array(board->array, board->rows, board->columns);
    

    destroy_int_array(board);
    free(monomer_locations);

    free(dd_energy_array);
    destroy_vec_array(derivative_array);
    destroy_vec_array(dd_array);

    //populate_array(board_states, rows, columns);
    //print_array(arr, rows, columns);


    // int_array* arr2 = create_int_array(5,3);
    // print_array(arr2->array, arr2->rows, arr2->columns);

    // destroy_int_array(arr2);


}
