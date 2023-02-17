#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"

int_array* create_int_array(int rows, int columns) {
    int_array* arr = malloc(sizeof(int_array));
    arr->array = (int**) malloc( rows * sizeof(int*));
    for (int i = 0; i < rows; i++){
        arr->array[i] = (int*) calloc(columns, sizeof(int));
    }
    arr->rows = rows;
    arr->columns = columns;
    return arr;
}

void destroy_int_array(int_array* arr) {
    for (int i = 0; i < arr->rows; i++) {
        free(arr->array[i]);
    }
    free(arr->array);
    free(arr);
}

vec_array* create_vec_array(int num_vecs) {
    vec_array* arr = malloc(sizeof(vec_array));
    arr->array = (vec*) malloc(num_vecs * sizeof(vec));
    arr->n_elements = num_vecs;
    return arr; 
}

void destroy_vec_array(vec_array *arr) {
    free(arr->array);
    free(arr);
}

vec create_vec(double i_ind, double j_ind) {
  vec z;
  z.i_ind = i_ind;
  z.j_ind = j_ind;
  return z;
}

vec add_vec(vec x, vec y) {
  vec z;
  z.i_ind = x.i_ind + y.i_ind;
  z.j_ind = x.j_ind + y.j_ind;
  return z;
}

vec subtract_vec(vec x, vec y) {
    vec z;
    z.i_ind = x.i_ind - y.i_ind;
    z.j_ind = x.j_ind - y.j_ind;
    return z;
}

vec normalize_vec(vec x) {
  vec n;
  double norm = sqrt(x.i_ind * x.i_ind + x.j_ind * x.j_ind);
  n.i_ind = x.i_ind/norm;
  n.j_ind = x.j_ind/norm;
  return n;
}

double square_vec(vec x) {
    return x.i_ind*x.i_ind + x.j_ind*x.j_ind;
}

ivec create_ivec(int i_ind, int j_ind) {
  ivec n;
  n.i_ind = i_ind;
  n.j_ind = j_ind;
  return n;
}

ivec add_ivec(ivec a, ivec b) {
  ivec summed;
  summed.i_ind = a.i_ind + b.i_ind;
  summed.j_ind = a.j_ind + b.j_ind;
  return summed;
}

ivec subtract_ivec(ivec a, ivec b) {
  ivec subtract;
  subtract.i_ind = a.i_ind - b.i_ind;
  subtract.j_ind = a.j_ind - b.j_ind;
  return subtract;
}

ivec scalar_multiplication(int lambda, ivec x) {
  ivec product; 
  product.i_ind = x.i_ind * lambda;
  product.j_ind = x.j_ind * lambda;
  return product;
}

vec ivec_to_vec(ivec a) {
    vec vec_a;
    vec_a.i_ind = (double) a.i_ind;
    vec_a.j_ind = (double) a.j_ind;
    return vec_a;
}
