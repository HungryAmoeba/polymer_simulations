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

    int n_steps =  1 + 1E7;              // number of monte carlo sweeps
    int save_freq = 1E6;
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

    for (size_t i = 0; i < N; i++) {
        if (i == 0) {
            double i_diff = states[i + 1][0] - states[i][0];
            double j_diff = states[i + 1][1] - states[i][1];
            double norm = sqrt(i_diff*i_diff + j_diff * j_diff);
            d[i][0] = i_diff/norm;
            d[i][1] = j_diff/norm;
        } else if (i == N-1) {
            double i_diff = states[i][0] - states[i-1][0];
            double j_diff = states[i][1] - states[i-1][1];
            double norm = sqrt(i_diff*i_diff + j_diff * j_diff);
            d[i][0] = i_diff/norm;
            d[i][1] = j_diff/norm;            
        } else {
            // calculate derivatives using the normal formula
            double fti = states[i][0] - states[i-1][0];
            double ftj = states[i][1] - states[i-1][1];
            double norm = sqrt(fti*fti + ftj*ftj);
            double ni = fti/norm; double nj = ftj/norm;

            double sti = states[i+1][0] - states[i][0];
            double stj = states[i+1][1] - states[i][1];
            norm = sqrt(sti*sti + stj*stj);
            ni = ni + sti/norm; nj = nj + stj/norm;

            norm = sqrt(ni*ni + nj*nj);
            d[i][0] = ni/norm; d[i][1] = nj/norm;
        }
        //printf("Derivative at ind %zu is %f, %f\n", i, d[i][0], d[i][1]);
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

            for (int s = 0; s < 8; s++) {
                if (!made_mv) {
                    int *mv_arr = get_move(so[s]);
                    int legal_move = check_legal(ind, states, lat, mv_arr, lw, N);
                    if (legal_move) {
                        made_mv = 1;
                        final_move_ind = s;
                        //printf("in here, %d, %d\n", mv_arr[0],mv_arr[1]);
                    } 
                    free(mv_arr);
                }
            }

            if (made_mv) { 
                //printf("%d",final_move_ind);
                int *pm = get_move(so[final_move_ind]);
                //printf("proposed move is %d,%d\n", pm[0], pm[1]);
                // calculate the energy difference of the new move
                int pi = states[ind][0] + pm[0];
                int pj = states[ind][1] + pm[1];
                //printf("at ind %d, previous is %d, %d, proposed is %d, %d\n", ind, states[ind][0], states[ind][1], pi, pj);
                if (ind == 0) {
                    double pdi = states[ind + 1][0] - pi; // proposed derivative i term
                    double pdj = states[ind + 1][1] - pj; // proposed derivative j term
                    double norm = sqrt(pdi*pdi + pdj*pdj);
                    pdi = pdi/norm; pdj = pdj/norm;
                    double pddi = d[ind + 1][0] - pdi; // difference in derivative i term
                    double pddj = d[ind + 1][1] - pdj; // difference in derivative j term
                    // compare to previous energy
                    double prev_e = (dd[ind][0]*dd[ind][0] + dd[ind][1]*dd[ind][1]) * l_p/(2 * b); 
                    double prop_e = (pddi*pddi + pddj*pddj) * l_p/(2 * b);
                    double delta_e = prop_e - prev_e;
                    double runif =((double) rand()/(double)RAND_MAX);
                    total_reject++;
                    if (exp(-delta_e) > runif) {
                        //printf("accepting the state at ind 0\n");
                        // accept the state
                        lat[states[ind][0]][states[ind][1]] = 0;
                        states[ind][0] = pi;
                        states[ind][1] = pj;
                        lat[states[ind][0]][states[ind][1]] = 1;
                        d[ind][0] = pdi; d[ind][1] = pdj;
                        dd[ind][0] = pddi; dd[ind][1] = pddj;
                        ce = ce + delta_e;
                        total_accept++;
                        total_reject--;   
                    }
                } else if (ind == N-1) {
                    double pdi = pi - states[ind-1][0];
                    double pdj = pj - states[ind-1][1];
                    double norm = (pdi*pdi + pdj*pdj);
                    pdi = pdi/norm; pdj = pdj/norm;
                    double pddi = pdi - d[ind-1][0];
                    double pddj = pdi - d[ind-1][1];
                    double prev_e = (dd[ind-1][0]*dd[ind-1][0] + dd[ind-1][1]*dd[ind-1][1]) * l_p/(2 * b); 
                    double prop_e = (pddi*pddi + pddj*pddj) * l_p/(2 * b);
                    double delta_e = prop_e -prev_e;
                    double runif = ((double) rand()/(double)RAND_MAX);
                    total_reject++;
                    if (exp(-delta_e) > runif) {
                        // accept the state
                        lat[states[ind][0]][states[ind][1]] = 0;
                        states[ind][0] = pi;
                        states[ind][1] = pj;
                        lat[states[ind][0]][states[ind][1]] = 1;
                        d[ind][0] = pdi; d[ind][1] = pdj;
                        dd[ind-1][0] = pddi; dd[ind-1][1] = pddj;
                        ce = ce + delta_e;
                        total_accept++;
                        total_reject--;                     
                    }
                } else {
                    // calculate derivatives using the normal formula
                    double fti = pi - states[ind-1][0]; // first term i
                    double ftj = pj - states[ind-1][1]; // first term j
                    double norm = sqrt(fti*fti + ftj*ftj);
                    double ni = fti/norm; double nj = ftj/norm; //normalized i, j

                    double sti = states[ind+1][0] - pi; // second term i
                    double stj = states[ind+1][1] - pj; // second term j
                    norm = sqrt(sti*sti + stj*stj);
                    ni = ni + sti/norm; nj = nj + stj/norm;

                    norm = sqrt(ni*ni + nj*nj);
                    double pdi = ni/norm; double pdj = nj/norm;

                    // right of monomer changes
                    double R_pddi = d[ind + 1][0] - pdi; // difference in derivative i term
                    double R_pddj = d[ind + 1][1] - pdj; // difference in derivative j term
                    // compare to previous energy
                    double R_prev_e = (dd[ind][0]*dd[ind][0] + dd[ind][1]*dd[ind][1]) * l_p/(2 * b); // rightside previous energy
                    double R_prop_e = (R_pddi*R_pddi + R_pddj*R_pddj) * l_p/(2 * b); // rightside proposed energy
                    double R_delta_e = R_prop_e - R_prev_e; // rightside proposed energy difference

                    // left of monomer changes
                    double L_pddi = pdi - d[ind-1][0];
                    double L_pddj = pdi - d[ind-1][1];
                    double L_prev_e = (dd[ind-1][0]*dd[ind-1][0] + dd[ind-1][1]*dd[ind-1][1]) * l_p/(2 * b); 
                    double L_prop_e = (L_pddi*L_pddi + L_pddj*L_pddj) * l_p/(2 * b);
                    double L_delta_e = L_prop_e - L_prev_e;

                    double delta_e = R_delta_e + L_delta_e;
                    double runif = ((double) rand()/(double)RAND_MAX);
                    total_reject++;
                    if (exp(-delta_e) > runif) {
                        // accept the state
                        lat[states[ind][0]][states[ind][1]] = 0;
                        states[ind][0] = pi;
                        states[ind][1] = pj;
                        lat[pi][pj] = 1;
                        d[ind][0] = pdi; d[ind][1] = pdj;
                        dd[ind-1][0] = L_pddi; dd[ind-1][1] = L_pddj;
                        dd[ind][0] = R_pddi; dd[ind][1] = R_pddj;
                        ce = ce + delta_e;   
                        total_accept++;
                        total_reject--;                     
                    }
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
    printf("success! \n");
    return 0;
}

