#!/bin/bash
#SBATCH --job-name=thin_serial_static_big
#SBATCH --output=thin_serial_static_big_grid.csv
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=THIN
#SBATCH --nodelist=thin010

module load openMPI/4.1.5/gnu/12.2.1

cd /u/dssc/s271711/HPC-final/exercise1

echo "time"

thin_bin/conway_serial_pgm.out -r -f big_grid -n 100 -s 0 -e 1

thin_bin/algo2_conway_serial_pgm.out -r -f big_grid -n 100 -s 0 -e 1

