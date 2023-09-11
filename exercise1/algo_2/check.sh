#!/bin/bash
make 

echo "checking coherency... with snark loop"

echo "checking results for static evolution are coherent"

./algo2_conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.serial

export OMP_NUM_THREADS=16
./algo2_conway_omp_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.omp

mpirun -np 8 algo2_conway_mpi_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.mpi

echo "comparing serial and omp"
cmp snapshot_00100.serial snapshot_00100.omp >&1
echo "comparing omp and mpi"
cmp snapshot_00100.omp snapshot_00100.mpi >&1

rm snapshot_*

echo "checking coherency... with glider gun"

echo "checking results for static evolution are coherent"

./algo2_conway_serial_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.serial

export OMP_NUM_THREADS=16
./algo2_conway_omp_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.omp

mpirun -np 8 algo2_conway_mpi_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.mpi

echo "comparing serial and omp"
cmp snapshot_00400.serial snapshot_00400.omp >&1
echo "comparing omp and mpi"
cmp snapshot_00400.omp snapshot_00400.mpi >&1

rm snapshot_*

make clean
