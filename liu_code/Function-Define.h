#ifndef Fucntion-Define-H
#define Fucntion-Define-H

#define NumPolymer 2   // the number of polymers
#define NumMonomer 101 // the number of monomers for each chain; monomers are labeled from 0 to 100 and the number of bonds thus is 100

#define NumMonomerFloat 101.0 // the number of monomers for accurate calculations


#define TotalMCS   50000 //the total Monte Carlo steps
#define RouseRelaxationTime 50000
// the time for relaxing chains; data are not collected until this time


#define RelaxTime 50000  // the time interval between two sampled data
// to get rid of the correlation; it is usually chosen as the same
// magnitude of Rouse relaxation time
//

#define ConfigWaitTime 1000 // the frequency to record the coordinatesof monomers
#define InteractionCoefficient2 0.1  // a pseudo harmonic interaction coefficient; it is used for overlapping two chains
#define Stiffness  1     // the stiffness of chains; it is used for semiflexible chains
#define Box_size_x 60    //the length of the confinement
			 //
#define Box_size_y 10

// for generating random numbers
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
extern long Seed;

typedef struct PolymerDef // Information of polymers
{
int MonomerX;
int MonomerY;
int MonomerFormerX;
// current x coordinate
// current y coordinate
// x coordinate of a monomer before MCS
int MonomerFormerY;

double BendingEnergy;
double BendingFormerEnergy;
double BondLengthSquare;
double BondLengthFormerSquare;
} Polymer;

Polymer Monomer[NumPolymer][NumMonomer];

double RandomFunc2(long *); // random number generator
long   SeedGenerator(void); // seed generator for RandomFunc2()
void   BondInitial(int,int); // use to generate the configurations of
polymers; (polymer, site)

int    BondInitialChecking(int,int); //check self avoidance during
// initializing chains; (polymer,site)
void   MonomerMoving(int,int,int,int);  // move monomers in MCS
//(polymer, site, choice,direction)
int    JudgeNewBond(int,int,double);   //check self avoidance for
//the new bonds in Monte Carlo steps; (polymer, site, random number)
int    JudgeBoundray(int,int,int,int);  // check the boundary conditions;
//(order,polymer,x,y), order=0 means checking the conditions during
//initializing chain. Order = 1 mean checking the conditions after that.

#endif

