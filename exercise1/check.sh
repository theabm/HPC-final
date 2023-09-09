#!/bin/bash
make 

./conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.serial

./conway_serial_pgm.v2.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.serial.v2

./conway_omp_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.omp

mpirun -np 2 conway_mpi_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.mpi

cmp snapshot_00100.serial snapshot_00100.omp >&1
cmp snapshot_00100.omp snapshot_00100.mpi >&1
cmp snapshot_00100.serial snapshot_00100.serial.v2 >&1


rm snapshot_*

./conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.serial

./conway_serial_pgm.v2.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.serial.v2

./conway_omp_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.omp

mpirun -np 2 conway_mpi_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.mpi

cmp snapshot_00100.serial snapshot_00100.omp >&1
cmp snapshot_00100.omp snapshot_00100.mpi >&1

rm snapshot_*

make clean
