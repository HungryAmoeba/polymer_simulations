#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "helper.h"

int main(int argc, char **argv) {

    int lw = 120;                       // lattice board width
    int N = 50;                         // number of monomers
    float l_p = 15.0;                   // persistence length
    float b = 2.8;                      // average bond length

    int states[N][2];                   // create 2D array for the monomers
    // int lat[lw][lw];                    // create lattice
    float d[N][2];                      // derivatives
    float dd[N-1][2];                   // store differences in adjacent derivatives

    int n_steps = 10000000;                    // number of steps
    float ce = 0;                       // current energy

    // create the lattice
    int** lat = (int**)malloc(lw *sizeof(int*));
    for (size_t i = 0; i < lw; i++) {
        lat[i] = (int*)malloc(lw * sizeof(int));
    } 

    // populate the lattice
    int fs = 0; // number of filled states
    for (size_t i = 0; i < lw; i++) {
        for (size_t j = 0; j < lw; j++) {
            if (fs < N && i % 3 == 0 && j % 3 == 0) {
                lat[i][j] = 1;
                states[fs][0] = i;
                states[fs][1] = j;
            }
            else {
                lat[i][j] = 0;
            }
        }
    }

    for (size_t i = 0; i < N; i++) {
        if (i == 0) {
            float i_diff = states[i + 1][0] - states[i][0];
            float j_diff = states[i + 1][1] - states[i][1];
            float norm = sqrt(i_diff*i_diff + j_diff * j_diff);
            d[i][0] = i_diff/norm;
            d[i][1] = j_diff/norm;
        } else if (i == N-1) {
            float i_diff = states[i][0] - states[i-1][0];
            float j_diff = states[i][1] - states[i-1][1];
            float norm = sqrt(i_diff*i_diff + j_diff * j_diff);
            d[i][0] = i_diff/norm;
            d[i][1] = j_diff/norm;            
        } else {
            // calculate derivatives using the normal formula
            float fti = states[i][0] - states[i-1][0];
            float ftj = states[i][1] - states[i-1][1];
            float norm = sqrt(fti*fti + ftj*ftj);
            float ni = fti/norm; float nj = ftj/norm;

            float sti = states[i+1][0] - states[i][0];
            float stj = states[i+1][1] - states[i][0];
            norm = sqrt(sti*sti + stj*stj);
            ni = ni + sti/norm; nj = nj + stj/norm;

            norm = sqrt(ni*ni + nj*nj);
            d[i][0] = ni/norm; d[i][1] = nj/norm;
        }
    }

    // calculate differences in derivatives
    // and update current energy
    
    for (size_t i = 0; i < N-1; i++) {
        dd[i][0] = d[i+1][0] - d[i][0];
        dd[i][1] = d[i+1][1] - d[i][1];
        ce += (dd[i][0] * dd[i][0] + dd[i][1] * dd[i][1]) * l_p/(2*b);
    }

    // do metropolis
    srand(12); // initialization for random number
    clock_t start, curr;
    double cpu_time_used;
    start = clock();
    for (size_t step = 0; step < n_steps; step++) {
        if (step % 1000000 == 0) {
            curr = clock();
            cpu_time_used = ((double) (curr - start)) / CLOCKS_PER_SEC;
            printf("Step is %zu after time %f\n", step, cpu_time_used);
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
                if (made_mv == 0) {
                    int *mv_arr = get_move(so[s]);
                    int legal_move = check_legal(ind, states, lat, mv_arr, lw);
                    if (legal_move) {
                        made_mv = 1;
                        final_move_ind = s;
                    } 
                    free(mv_arr);
                }
            }

            if (made_mv) {
                int *pm = get_move(so[final_move_ind]);
                // calculate the energy difference of the new move
                int pi = states[ind][0] + pm[0];
                int pj = states[ind][1] + pm[1];
                if (ind == 0) {
                    float pdi = states[ind + 1][0] - pi; // proposed derivative i term
                    float pdj = states[ind + 1][1] - pj; // proposed derivative j term
                    float norm = sqrt(pdi*pdi + pdj*pdj);
                    pdi = pdi/norm; pdj = pdj/norm;
                    float pddi = d[ind + 1][0] - pdi; // difference in derivative i term
                    float pddj = d[ind + 1][1] - pdj; // difference in derivative j term
                    // compare to previous energy
                    float prev_e = (dd[ind][0]*dd[ind][0] + dd[ind][1]*dd[ind][1]) * l_p/(2 * b); 
                    float prop_e = (pddi*pddi + pddj*pddj) * l_p/(2 * b);
                    float delta_e = prop_e - prev_e;
                    float runif =((float) rand()/(float)RAND_MAX);
                    if (exp(-delta_e) > runif) {
                        // accept the state
                        lat[states[ind][0]][states[ind][1]] = 0;
                        states[ind][0] = pi;
                        states[ind][0] = pj;
                        lat[states[ind][0]][states[ind][1]] = 1;
                        d[ind][0] = pdi; d[ind][1] = pdj;
                        dd[ind][0] = pddi; dd[ind][1] = pddj;
                        ce = ce + delta_e;
                    }
                } else if (ind == N-1) {
                    float pdi = pi - states[ind-1][0];
                    float pdj = pj - states[ind-1][1];
                    float norm = (pdi*pdi + pdj*pdj);
                    pdi = pdi/norm; pdj = pdj/norm;
                    float pddi = pdi - d[ind-1][0];
                    float pddj = pdi - d[ind-1][1];
                    float prev_e = (dd[ind-1][0]*dd[ind-1][0] + dd[ind-1][1]*dd[ind-1][1]) * l_p/(2 * b); 
                    float prop_e = (pddi*pddi + pddj*pddj) * l_p/(2 * b);
                    float delta_e = prop_e -prev_e;
                    float runif = ((float) rand()/(float)RAND_MAX);
                    if (exp(-delta_e) > runif) {
                        // accept the state
                        lat[states[ind][0]][states[ind][1]] = 0;
                        states[ind][0] = pi;
                        states[ind][0] = pj;
                        lat[states[ind][0]][states[ind][1]] = 1;
                        d[ind][0] = pdi; d[ind][1] = pdj;
                        dd[ind-1][0] = pddi; dd[ind-1][1] = pddj;
                        ce = ce + delta_e;                     
                    }
                } else {
                    // calculate derivatives using the normal formula
                    float fti = pi - states[ind-1][0]; // first term i
                    float ftj = pj - states[ind-1][1]; // first term j
                    float norm = sqrt(fti*fti + ftj*ftj);
                    float ni = fti/norm; float nj = ftj/norm; //normalized i, j

                    float sti = states[ind+1][0] - pi; // second term i
                    float stj = states[ind+1][1] - pj; // second term j
                    norm = sqrt(sti*sti + stj*stj);
                    ni = ni + sti/norm; nj = nj + stj/norm;

                    norm = sqrt(ni*ni + nj*nj);
                    float pdi = ni/norm; float pdj = nj/norm;

                    // right of monomer changes
                    float R_pddi = d[ind + 1][0] - pdi; // difference in derivative i term
                    float R_pddj = d[ind + 1][1] - pdj; // difference in derivative j term
                    // compare to previous energy
                    float R_prev_e = (dd[ind][0]*dd[ind][0] + dd[ind][1]*dd[ind][1]) * l_p/(2 * b); // rightside previous energy
                    float R_prop_e = (R_pddi*R_pddi + R_pddj*R_pddj) * l_p/(2 * b); // rightside proposed energy
                    float R_delta_e = R_prop_e - R_prev_e; // rightside proposed energy difference

                    // left of monomer changes
                    float L_pddi = pdi - d[ind-1][0];
                    float L_pddj = pdi - d[ind-1][1];
                    float L_prev_e = (dd[ind-1][0]*dd[ind-1][0] + dd[ind-1][1]*dd[ind-1][1]) * l_p/(2 * b); 
                    float L_prop_e = (L_pddi*L_pddi + L_pddj*L_pddj) * l_p/(2 * b);
                    float L_delta_e = L_prop_e - L_prev_e;

                    float delta_e = R_delta_e + L_delta_e;
                    float runif = ((float) rand()/(float)RAND_MAX);
                    if (exp(-delta_e) > runif) {
                        // accept the state
                        lat[states[ind][0]][states[ind][1]] = 0;
                        states[ind][0] = pi;
                        states[ind][0] = pj;
                        lat[states[ind][0]][states[ind][1]] = 1;
                        d[ind][0] = pdi; d[ind][1] = pdj;
                        dd[ind-1][0] = L_pddi; dd[ind-1][1] = L_pddj;
                        dd[ind][0] = R_pddi; dd[ind][1] = R_pddj;
                        ce = ce + delta_e;                     
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

    printf("success! \n");
}

