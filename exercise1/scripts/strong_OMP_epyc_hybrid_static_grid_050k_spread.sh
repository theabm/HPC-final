#!/bin/bash
#SBATCH --job-name=strong_OMP_epyc_hybrid_static_spread_050k
#SBATCH --output=strong_OMP_epyc_hybrid_static_grid_050k_spread.csv
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=2
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=64
#SBATCH --time=2:00:00
#SBATCH --partition=EPYC
#SBATCH --exclude=epyc002,epyc003

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthreads_per_process,total_nthreads,time"

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

for j in 1 8 16 24 32 40 48 56 64 ; do
export OMP_NUM_THREADS=$j
for i in {0..2} ; do
 	mpirun --map-by socket epyc_bin/conway_hybrid_pgm.v2.out -r -f grid_050k -n 200 -s 0 -e 1
done 
done
