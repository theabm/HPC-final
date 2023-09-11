#!/bin/bash
make 

echo "checking coherency... with snark loop"

echo "checking results for static evolution are coherent"

./conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.serial

export OMP_NUM_THREADS=16
./conway_omp_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.omp

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

rm snapshot_*

echo "checking results for ordered evolution are coherent"

./conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.serial

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

rm snapshot_*

echo "checking coherency... with glider gun"

echo "checking results for static evolution are coherent"

./conway_serial_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.serial

export OMP_NUM_THREADS=16
./conway_omp_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.omp

mpirun -np 8 conway_mpi_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.mpi

export OMP_NUM_THREADS=4
export OMP_PLACES=cores
mpirun --map-by socket -np 2 conway_hybrid_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.hybrid

export OMP_NUM_THREADS=16

echo "comparing serial and omp"
cmp snapshot_00400.serial snapshot_00400.omp >&1
echo "comparing omp and mpi"
cmp snapshot_00400.omp snapshot_00400.mpi >&1
echo "comparing mpi and hybrid"
cmp snapshot_00400.mpi snapshot_00400.hybrid >&1

rm snapshot_*

echo "checking results for ordered evolution are coherent"

./conway_serial_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 0
mv snapshot_00400 snapshot_00400.serial

./conway_omp_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 0
mv snapshot_00400 snapshot_00400.omp

mpirun -np 2 conway_mpi_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 0
mv snapshot_00400 snapshot_00400.mpi

mpirun -np 8 conway_hybrid_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 0
mv snapshot_00400 snapshot_00400.hybrid

echo "comparing serial and omp"
cmp snapshot_00400.serial snapshot_00400.omp >&1
echo "comparing omp and mpi"
cmp snapshot_00400.omp snapshot_00400.mpi >&1
echo "comparing mpi and hybrid"
cmp snapshot_00400.mpi snapshot_00400.hybrid >&1

rm snapshot_*


make clean
