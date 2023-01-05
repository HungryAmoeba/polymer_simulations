#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helper.h"

int main(int argc, char **argv) {

    int lw = 50;                       // lattice board width
    int N = 50;                         // number of monomers
    double l_p = 15.0;                   // persistence length
    double b = 2.8;                      // average bond length

    int states[N][2];                   // create 2D array for the monomers
    // int lat[lw][lw];                 // create lattice
    double d[N][2];                      // derivatives
    double dd[N-1][2];                   // store differences in adjacent derivatives

    int n_steps =  1 + 1E7;//1 + 1E7;              // number of monte carlo sweeps
    int save_freq = 1E6;//1E6;
    double ce = 0;                       // current energy
    int total_accept = 0;
    int total_reject = 0;

    FILE *fptr = fopen("simulation_results/results.txt", "w");
    if (fptr == NULL) {
        printf("Error in opening file pointer \n");
        exit(1);
    }
    printf("Found file to write to\n");

    // create the lattice
    int** lat = (int**)malloc(lw *sizeof(int*));
    for (size_t i = 0; i < lw; i++) {
        lat[i] = (int*)malloc(lw * sizeof(int));
    } 

    // populate the lattice
    int fs = 0; // number of filled states
    for (int i = 0; i < lw; i++) {
        if (i%2 == 0) {
            for (int j = 0; j < lw; j++) {
                if (fs < N && i % 3 == 0 && j % 3 == 0) {
                    lat[i][j] = 1;
                    states[fs][0] = i;
                    states[fs][1] = j;
                    fs++;
                } else {
                    lat[i][j] = 0;
                }
            }
        }

        if (i %2 == 1) {
            for (int j = (lw - 1); j >= 0; j--) {
                if (fs < N && i % 3 == 0 && j % 3 == 0) {
                    lat[i][j] = 1;
                    states[fs][0] = i;
                    states[fs][1] = j;
                    fs++;
                } else {
                    lat[i][j] = 0;
                }
            }
        }
    }

    // calculate initial derivatives and differences in derivatives
    for (int i = 0; i < N; i++) {
        // first handle edge cases (start or end of array)
        double *derivative = get_derivative(i, N, states);
        d[i][0] = derivative[0];
        d[i][1] = derivative[1];
        free(derivative);
    }

    // calculate differences in derivatives
    // and update current energy
    
    for (size_t i = 0; i < N-1; i++) {
        dd[i][0] = d[i+1][0] - d[i][0];
        dd[i][1] = d[i+1][1] - d[i][1];
        //printf("For index %zu, additional term is %f \n",i,(dd[i][0] * dd[i][0] + dd[i][1] * dd[i][1]) * l_p/(2*b) );
        //printf("dd here is %f, %f\n", dd[i][0], dd[i][1]);
        ce += (dd[i][0] * dd[i][0] + dd[i][1] * dd[i][1]) * l_p/(2*b);
    }

    // do metropolis
    srand(12); // initialization for random number. set to 12 for reproducibility
    clock_t start, curr;
    double cpu_time_used;
    start = clock();
    for (size_t step = 0; step < n_steps; step++) {
        if (step % save_freq == 0) {
            curr = clock();
            cpu_time_used = ((double) (curr - start)) / CLOCKS_PER_SEC;
            printf("Step is %zu after time %f, ce is %f\n", step, cpu_time_used, ce);
            double recalculated_energy = calc_energy(states, N, l_p, b);
            printf("Current recalculated energy is %f, difference is %lf\n", recalculated_energy, ce-recalculated_energy);
            // save states to the file
            fprintf(fptr, "Step is %zu \n", step);
            for (size_t i = 0; i < N; i++) {
                fprintf(fptr, "%d,%d\n", states[i][0],states[i][1]);
            }
        }
        // choose N points to shift
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
                int legal_move = check_legal(ind, states, lat, mv_arr, lw, N);
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
                int original_i = states[ind][0];
                int original_j = states[ind][1];
                int pi = states[ind][0] + pm[0];
                int pj = states[ind][1] + pm[1];
                states[ind][0] = pi;
                states[ind][1] = pj;

                double prop_energy = calc_energy(states, N, l_p, b);
                double delta_e = prop_energy - ce;
                double runif = (double)rand()/(double)RAND_MAX;
                if (exp(-delta_e) > runif) {
                    ce = prop_energy;
                    // recalculate the necessary derivatives
                    for (int i = fmax(0, ind-1); i < fmin(N-1, ind + 1); i++) {
                        double *derivative = get_derivative(i, N, states);
                        d[ind][0] = derivative[0];
                        d[ind][1] = derivative[1];
                        free(derivative);
                    }
                    // recalculate the necessary differences in derivative
                    for (int i = fmax(0, ind-2); i < fmin(N-2, ind + 1); i++) {
                        dd[i][0] = d[i + 1][0] - d[i][0];
                        dd[i][1] = d[i + 1][1] - d[i][1];
                    }
                    total_accept++;
                } else {
                    states[ind][0] = original_i;
                    states[ind][1] = original_j;
                    total_reject++;
                }
                free(pm);

            }
            free(so);
        }
    }

    for (size_t i = 0; i < lw; i++) {
        free(lat[i]);
    }
    free(lat);
    fclose(fptr);
    printf("total running has %d acceptances, %d rejections\n", total_accept, total_reject);
    printf("the acceptance rate is %f\n", (float) total_accept/(total_accept + total_reject));
    return 0;
}

