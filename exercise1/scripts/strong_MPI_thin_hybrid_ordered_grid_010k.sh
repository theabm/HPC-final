#!/bin/bash
#SBATCH --job-name=strong_MPI_thin_hybrid_ordered_grid_010k
#SBATCH --output=strong_MPI_thin_hybrid_ordered_grid_010k.csv
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=THIN
#SBATCH --nodelist=thin007

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthreads_per_process,total_nthreads,time"

export OMP_NUM_THREADS=1

for j in 1 4 8 12 16 20 24 ; do
for i in {0..4} ; do
 	mpirun -np $j thin_bin/conway_hybrid_pgm.v2.out -r -f grid_010k -n 200 -s 0 -e 0
done 
done
