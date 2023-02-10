#ifndef __VECTORS_H__
#define __VECTORS_H__

#include <stdio.h>

struct vec;
struct vec {
  double i_ind;
  double j_ind;
};

typedef struct vec vec;

struct int_array;
struct int_array {
    int **array;
    int rows;
    int columns;
};
typedef struct int_array int_array;

struct vec_array;
struct vec_array {
  vec* array;
  int n_elements;
};
typedef struct vec_array vec_array;

struct ivec; // integer vec
struct ivec {
  int i_ind;
  int j_ind;
};
typedef struct ivec ivec;

vec_array* create_vec_array(int num_vecs);

void destroy_vec_array(vec_array *arr);

int_array* create_int_array(int rows, int columns);

void destroy_int_array(int_array *arr);

vec create_vec(double i_ind, double j_ind);
vec add_vec(vec x, vec y);
vec subtract_vec(vec x, vec y);
vec normalize_vec(vec x);
double square_vec(vec x);

ivec create_ivec(int i_ind, int j_ind);
ivec add_ivec(ivec a, ivec b);
ivec subtract_ivec(ivec a, ivec b);
ivec scalar_multiplication(int lambda, ivec x);

vec ivec_to_vec(ivec a);

#endif
