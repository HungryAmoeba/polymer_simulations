#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

int check_legal(int index, int states[][2], int **lat, int *move_arr, int lw, int N) {
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

  // check for excluded volume
  // i am so sorry for this horrible code style
  for (int i = -1; i < 2; i ++) {
    for (int j = -1; j < 2; j ++) {
      if (i != -move_arr[0] || j != -move_arr[1]) {
        if (i + pi >= 0 && i + pi < lw) {
          if (j + pj >= 0 && j + pj < lw) {
            if (lat[pi +i][pj + j] == 1) {
              return 0;
            }
          }
        }
        
      }
    }
  }

  // check length does not exceed 3
  if (index != 0) {
    // compare to previous monomer
    int dist_x = pi - states[index-1][0];
    int dist_y = pj - states[index-1][1];
    //printf("%d is pi, %d is prev, %d is sq\n", pi, states[index-1][0], dist_x*dist_x);
    //printf("%d is the distance on left \n", dist_x*dist_x + dist_y*dist_y);
    if (dist_x*dist_x + dist_y*dist_y > 13) {
      //printf("rejected \n \n");
      return 0;
    }
  }
  if (index != N-1) {
    // compare to next monomer
    int dist_x = pi - states[index+1][0];
    int dist_y = pj - states[index+1][1];
    if (dist_x*dist_x + dist_y*dist_y > 13) {
      return 0;
    } 
  }
  return 1;
}

double calc_energy(int states[][2], int N, double l_p, double b) {
  double ce = 0;
  double d[N][2];
  double dd[N-1][2];

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

    return ce;
}
