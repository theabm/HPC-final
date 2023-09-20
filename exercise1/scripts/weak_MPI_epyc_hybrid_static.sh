#!/bin/bash
#SBATCH --job-name=weak_MPI_epyc_hybrid_static
#SBATCH --output=weak_MPI_epyc_hybrid_static.csv
#SBATCH --nodes=4
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=EPYC
#SBATCH --exclude=epyc002,epyc003

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthreads_per_process,total_nthreads,time"

export OMP_NUM_THREADS=64
export OMP_PLACES=cores
export OMP_PROC_BIND=close

for j in 1 2 3 4 5 6 7 8 ; do
for i in {0..9} ; do
 	mpirun --map-by socket -np $j epyc_bin/conway_hybrid_pgm.v2.out -r -f grid_${j} -n 200 -s 0 -e 1
done 
done
