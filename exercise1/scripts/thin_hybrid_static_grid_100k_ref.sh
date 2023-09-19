#!/bin/bash
#SBATCH --job-name=thin_hybrid_static_grid_100k_ref
#SBATCH --output=thin_hybrid_static_grid_100k_ref.csv
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=THIN
#SBATCH --nodelist=thin010

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthreads_per_process,total_nthreads,time"

export OMP_NUM_THREADS=1
mpirun -np 1 thin_bin/conway_hybrid_pgm.v2.out -r -f grid_100k -n 100 -s 0 -e 1


