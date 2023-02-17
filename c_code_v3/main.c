#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helper.h"

int main(int argc, char **argv) {
    // expect to get Main rows, columns, N
    int n_steps = 1 + 1E8; // Should use 1  + E9 for full experiments
    int save_freq = 1E7; // Frequency at which the board state gets written
    
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
    
    // open file to write to
    FILE *fptr = fopen("simulation_results/results.txt", "w");
    if (fptr == NULL) {
        printf("Error in opening file pointer \n");
        exit(1);
    }
    printf("Found file to write to\n");

    System* sys = create_system(N, rows, cols, l_p, b);
    initialize_starting_state(sys);
    System* sys_copy = deep_copy_system(sys);

    srand(12); // initialization for random number. set to 12 for reproducibility
    clock_t start, curr;
    double cpu_time_used;
    start = clock();

    // Perform n_steps monte carlo sweeps
    for (int step = 0; step < n_steps; step++) {
        // every save_freq steps, record the current state of the system
        if (step%save_freq == 0) {
            curr = clock();
            cpu_time_used = ((double) (curr - start)) / CLOCKS_PER_SEC;
            printf("Step is %d after time %f, ce is %f\n", step, cpu_time_used, sys->curr_energy);
            // TODO: check for floating point drift and fix if necessary

            // double recalculated_energy = calc_energy(states, N, l_p, b);
            // printf("Current recalculated energy is %f, difference is %lf\n", recalculated_energy, ce-recalculated_energy);
            // save states to the file
            fprintf(fptr, "Step is %d \n", step);
            for (size_t i = 0; i < N; i++) {
                fprintf(fptr, "%d,%d\n", sys->monomer_locations[i], sys->monomer_locations[i + sys->n_monomers]);
            }
            print_board(sys);
        }
        
        // for each step, perform n_monomers proposals
        for (size_t pt = 0; pt < N; pt ++) {
            int ind = rand() % N;
            // generate potential indicies 
            int *so = malloc(8*sizeof(int));
            for (int s = 0; s < 8; s ++) {
                so[s] = s;
            }
            shuffle(so,8);
            int made_mv = 0;
            int final_move_ind = -1;

            int s = 0;

            while ((!made_mv) && (s < 8)) {
                int *mv_arr = get_move(so[s]);
                int legal_move = check_legal(ind, sys, mv_arr);
                if (legal_move) {
                    made_mv = 1;
                    final_move_ind = s;
                }
                free(mv_arr);
                s++;
            }

            if (made_mv) {
                int *pm = get_move(so[final_move_ind]);
                // calculate the energy difference of the new move
                int original_i = sys->monomer_locations[ind];
                int original_j = sys->monomer_locations[ind + sys->n_monomers];
                int pi = original_i + pm[0];
                int pj = original_j + pm[1];

                sys_copy->monomer_locations[ind] = pi;
                sys_copy->monomer_locations[ind + sys->n_monomers] = pj;

                // update the derivatives of the new array
                double proposed_e = update_system(sys_copy, ind);

                double delta_e = proposed_e - sys->curr_energy;
                double runif = (double)rand()/(double)RAND_MAX;

                if (exp(-delta_e) > runif) {
                    // accept the state

                    // update the lattice
                    sys->lattice[original_i * sys->board_cols + original_j] = 0;
                    sys->lattice[pi * sys->board_cols + pj] = 1;

                    // update sys from the copy
                    sys_copy_over(sys_copy, sys, ind);
                } else {
                    sys_copy_over(sys, sys_copy, ind);
                }

            }

            free(so);
        }
    }


    destroy_system(sys_copy);
    destroy_system(sys);

}

