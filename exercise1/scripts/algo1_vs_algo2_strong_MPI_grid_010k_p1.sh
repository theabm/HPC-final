#!/bin/bash
#SBATCH --job-name=a1_v_a2_strong_MPI_grid_010k_p1
#SBATCH --output=a1_v_a2_strong_MPI_grid_010k_p1.csv
#SBATCH --nodes=3
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=THIN
#SBATCH --nodelist=thin007,thin008,thin010

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "size,nthread_per_process,total_nthreads,time"

export OMP_NUM_THREADS=1

for i in {0..2} ; do 
	mpirun -np 1 thin_bin/conway_hybrid_pgm.v2.out -r -f grid_010k -n 200 -s 0 -e 1
done

for i in {0..2} ; do 
	mpirun -np 1 thin_bin/conway_hybrid_pgm.v2.out -r -f glider_010k.pgm -n 200 -s 0 -e 1
done

for j in 12 24 36 48 60 72; do
for i in {0..9} ; do
 	mpirun -np $j thin_bin/conway_hybrid_pgm.v2.out -r -f grid_010k -n 200 -s 0 -e 1
done 
done

for j in 12 24 36 48 60 72; do
for i in {0..9} ; do
 	mpirun -np $j thin_bin/conway_hybrid_pgm.v2.out -r -f glider_010k.pgm -n 200 -s 0 -e 1
done 
done
