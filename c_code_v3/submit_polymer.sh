#!/bin/bash

#SBATCH --job-name=my_sim_job
#SBATCH --output=results_%A_%a.out
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --reservation=sumry2023
#SBATCH --cpus-per-task=1
#SBATCH --array=1-36

# Load any modules or software if necessary

LROWS_ARR=(30 60)
LCOLS_ARR=(30 60)
MONOMERS_ARR=(50 100 150)
LP_ARR=(1.0 2.0 4.0)
GEOMETRY=2

index=$SLURM_ARRAY_TASK_ID

lrows_idx=$(( (index-1) / 18 ))
lcols_idx=$(( (index-1 - 18*lrows_idx) / 9 ))
monomers_idx=$(( (index-1 - 18*lrows_idx - 9*lcols_idx) / 3 ))
lp_idx=$(( (index-1 - 18*lrows_idx - 9*lcols_idx - 3*monomers_idx) ))
geometry_idx=0  # since GEOMETRY is constant at 2

./Main ${LROWS_ARR[$lrows_idx]} ${LCOLS_ARR[$lcols_idx]} ${MONOMERS_ARR[$monomers_idx]} ${LP_ARR[$lp_idx]} $GEOMETRY
