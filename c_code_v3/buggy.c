#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct System
{
    int n_monomers;             // N
    int board_rows;             // The number of rows
    int board_cols;             // number of columns
    int* monomer_locations;     // Array storing the coordinates of the monomers
} System;

System* create_system(int N, int board_rows, int board_cols) {
    System* sys = malloc(sizeof(System));
    sys->n_monomers = N;
  
    sys->monomer_locations = calloc(2 * N , sizeof(int));
    return sys;
}

void initialize_starting_state(System* sys) {
    sys->monomer_locations[0] = 2; // sys->board_rows / 2;
    printf("made first assignment\n");
    sys->monomer_locations[sys->n_monomers] = 2; //sys->board_cols / 2;
    printf("made starting\n");
    int monomer_ind = 0;
    int segment_length = 1;
    int segment_counter = 0;
    int segment_leg = 0;
    // draw a spiral originating from the center of the board
    for (int i = 1; i < sys->n_monomers; i++) {
        int direction_x = (int) -cos(M_PI/2 * (double) segment_counter);
        int direction_y = (int) sin(M_PI/2 * (double) segment_counter); 
        sys->monomer_locations[i] = sys->monomer_locations[i-1] + 2 * direction_x;
        printf("Assigned value %d\n", sys->monomer_locations[i-1 + sys->n_monomers]);
        printf("Will assign value %d\n", 2 + sys->monomer_locations[i -1 + sys->n_monomers]);
        sys->monomer_locations[i + sys->n_monomers] = sys->monomer_locations[ i -1 + sys->n_monomers] + 2 * direction_y;
        monomer_ind++;
        if (monomer_ind == segment_length) {
            segment_counter++; 
            monomer_ind = 0;
            segment_leg++;
        }
        if (segment_leg == 2) {
            segment_leg = 0;
            segment_length++;
        }
    }
}

int main(void) {
    System* sys = create_system(20, 25, 25);
    printf("made system\n");
    initialize_starting_state(sys);
    free(sys->monomer_locations);
    free(sys);    
}
