#ifndef __HELPER_H__
#define __HELPER_H__

#include <stdio.h>

struct vec;
typedef struct vec vec;

vec add_vec(vec x, vec y);
vec normalize_vec(vec x);


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

int check_legal(int index, int states[][2], int **lat, int *move_arr, int lw, int N);

double *get_derivative(int i, int N, int states[][2]);

double calc_energy(int states[][2], int N, double l_p, double b);

#endif
