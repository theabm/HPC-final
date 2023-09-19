#!/bin/bash
#SBATCH --job-name=strong_OMP_thin_hybrid_static_extra_3
#SBATCH --output=strong_OMP_thin_hybrid_static_extra_3.csv
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=2
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2:00:00
#SBATCH --partition=THIN
#SBATCH --nodelist=thin007

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthreads_per_process,total_nthreads,time"

export OMP_PLACES=cores
export OMP_PROC_BIND=close

for j in 1 2 4 6 8 10 12 ; do
export OMP_NUM_THREADS=$j
for i in {0..4} ; do
 	mpirun --map-by socket thin_bin/conway_hybrid_pgm.v2.out -r -f grid_3 -n 200 -s 0 -e 1
done 
done
