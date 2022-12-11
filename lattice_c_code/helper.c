#include <stdio.h>
#include <stdlib.h>

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

int check_legal(int index, int states[][2], int **lat, int *move_arr, int lw) {
  // check that the proposed move is not out of bounds
  int pi = states[index][0] + move_arr[0];
  int pj = states[index][1] + move_arr[1];
  if (pi < 0 || pi >= lw) {
    return 0;
  }
  if (pj < 0 || pj >= lw) {
    return 0;
  }
  if (lat[pi][pj] == 1) {
    return 0;
  }
  // move is not out of bounds, and proposed square is not occupied
  for (size_t i = -1; i < 2; i ++) {
    for (size_t j = -1; j < 2; j ++) {
      if (i != -move_arr[0] || j != -move_arr[1]) {
        if (lat[pi +i][pj + j] == 1) {
          return 0;
        }
      }
    }
  }
  return 1;
}
