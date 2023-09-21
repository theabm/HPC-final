#!/bin/bash
#SBATCH --job-name=strong_MPI_epyc_hybrid_static_grid_070k
#SBATCH --output=strong_MPI_epyc_hybrid_static_grid_070k.csv
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=EPYC
#SBATCH --exclude=epyc002,epyc003

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthreads_per_process,total_nthreads,time"

export OMP_NUM_THREADS=1

for j in 64 128 192 256 320 384 448 512 ; do
for i in {0..4} ; do
 	mpirun -np $j epyc_bin/conway_hybrid_pgm.v2.out -r -f grid_070k -n 100 -s 0 -e 1
done 
done

mpirun -np 1 epyc_bin/conway_hybrid_pgm.v2.out -r -f grid_070k -n 100 -s 0 -e 1
mpirun -np 1 epyc_bin/conway_hybrid_pgm.v2.out -r -f grid_070k -n 100 -s 0 -e 1