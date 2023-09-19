#!/bin/bash
#SBATCH --job-name=strong_MPI_thin_hybrid_static_grid_070k
#SBATCH --output=strong_MPI_thin_hybrid_static_grid_070k.csv
#SBATCH --nodes=3
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=THIN
#SBATCH --nodelist=thin007,thin008,thin010

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthreads_per_process,total_nthreads,time"

export OMP_NUM_THREADS=1

for j in 72 60 48 36 24 12 ; do
for i in {0..4} ; do
 	mpirun -np $j thin_bin/conway_hybrid_pgm.v2.out -r -f grid_070k -n 100 -s 0 -e 1
done 
done
