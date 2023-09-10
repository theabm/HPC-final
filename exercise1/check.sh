#!/bin/bash
make 

echo "checking coherency..."

echo "checking results for static evolution are coherent"

./conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.serial

./conway_serial_pgm.v2.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.serial.v2

export OMP_NUM_THREADS=16
./conway_omp_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.omp

./conway_omp_pgm.v2.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.omp.v2

mpirun -np 8 conway_mpi_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.mpi

export OMP_NUM_THREADS=4
export OMP_PLACES=cores
mpirun --map-by socket -np 2 conway_hybrid_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.hybrid

export OMP_NUM_THREADS=16

echo "comparing serial and omp"
cmp snapshot_00100.serial snapshot_00100.omp >&1
echo "comparing omp and mpi"
cmp snapshot_00100.omp snapshot_00100.mpi >&1
echo "comparing mpi and hybrid"
cmp snapshot_00100.mpi snapshot_00100.hybrid >&1

echo "starting checks for v2"
echo "comparing serial v1 and serial v2"
cmp snapshot_00100.serial snapshot_00100.serial.v2 >&1
echo "comparing omp v1 and omp v2 "
cmp snapshot_00100.omp snapshot_00100.omp.v2 >&1

rm snapshot_*

echo "checking results for ordered evolution are coherent"

./conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.serial

./conway_serial_pgm.v2.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.serial.v2

./conway_omp_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.omp

mpirun -np 2 conway_mpi_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.mpi

mpirun -np 8 conway_hybrid_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.hybrid

echo "comparing serial and omp"
cmp snapshot_00100.serial snapshot_00100.omp >&1
echo "comparing omp and mpi"
cmp snapshot_00100.omp snapshot_00100.mpi >&1
echo "comparing mpi and hybrid"
cmp snapshot_00100.mpi snapshot_00100.hybrid >&1

echo "starting checks for v2"
echo "comparing serial v1 and serial v2"
cmp snapshot_00100.serial snapshot_00100.serial.v2 >&1


rm snapshot_*

make clean
