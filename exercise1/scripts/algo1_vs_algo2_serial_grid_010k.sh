#!/bin/bash
#SBATCH --job-name=a1_v_a2_serial_grid_010k
#SBATCH --output=a1_v_a2_serial_grid_010k.csv
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=THIN
#SBATCH --nodelist=thin007

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "time"

for i in {0..9} ; do
 	thin_bin/conway_serial_pgm.out -r -f grid_010k -n 200 -s 0 -e 1
done 

for i in {0..9} ; do
 	thin_bin/algo2_conway_serial_pgm.out -r -f grid_010k -n 200 -s 0 -e 1
done 

for i in {0..9} ; do
 	thin_bin/conway_serial_pgm.out -r -f glider_010k.pgm -n 200 -s 0 -e 1
done 

for i in {0..9} ; do
 	thin_bin/algo2_conway_serial_pgm.out -r -f glider_010k.pgm -n 200 -s 0 -e 1
done 
