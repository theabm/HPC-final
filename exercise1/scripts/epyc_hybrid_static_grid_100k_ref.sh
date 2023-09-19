#!/bin/bash
#SBATCH --job-name=epyc_serial_static_big_ref
#SBATCH --output=epyc_serial_static_big_grid_ref.csv
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=EPYC
#SBATCH --exclude=epyc002

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthreads_per_process,total_nthreads,time"

export OMP_NUM_THREADS=1
mpirun -np 1 epyc_bin/conway_hybrid_pgm.v2.out -r -f big_grid -n 100 -s 0 -e 1


