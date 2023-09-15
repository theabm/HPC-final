#!/bin/bash
make 
cd algo_2
make 
cd ..

echo "checking coherency... static evolution"

echo "snark loop..."

./conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.serial

./algo_2/algo2_conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.serial.algo2

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

mpirun --map-by socket -np 2 conway_hybrid_pgm.v2.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 1
mv snapshot_00100 snapshot_00100.hybrid.v2

export OMP_NUM_THREADS=16

echo "comparing serial and omp"
cmp snapshot_00100.serial snapshot_00100.omp >&1
echo "comparing omp and mpi"
cmp snapshot_00100.omp snapshot_00100.mpi >&1
echo "comparing mpi and hybrid"
cmp snapshot_00100.mpi snapshot_00100.hybrid >&1
echo "comparing omp v1 and v2"
cmp snapshot_00100.omp snapshot_00100.omp.v2 >&1
echo "comparing hybrid v2 and v2"
cmp snapshot_00100.hybrid snapshot_00100.hybrid.v2 >&1

echo "comparing serial algo1 and serial algo 2"
cmp snapshot_00100.serial snapshot_00100.serial.algo2 >&1

rm snapshot_*

echo "glider gun..."

./conway_serial_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.serial

./algo_2/algo2_conway_serial_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.serial.algo2

export OMP_NUM_THREADS=16
./conway_omp_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.omp

./conway_omp_pgm.v2.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.omp.v2

mpirun -np 8 conway_mpi_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.mpi

export OMP_NUM_THREADS=4
export OMP_PLACES=cores
mpirun --map-by socket -np 2 conway_hybrid_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.hybrid

mpirun --map-by socket -np 2 conway_hybrid_pgm.v2.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 1
mv snapshot_00400 snapshot_00400.hybrid.v2

export OMP_NUM_THREADS=16

echo "comparing serial and omp"
cmp snapshot_00400.serial snapshot_00400.omp >&1
echo "comparing omp and mpi"
cmp snapshot_00400.omp snapshot_00400.mpi >&1
echo "comparing mpi and hybrid"
cmp snapshot_00400.mpi snapshot_00400.hybrid >&1
echo "comparing omp v1 and v2"
cmp snapshot_00400.omp snapshot_00400.omp.v2 >&1
echo "comparing hybrid v1 and v2"
cmp snapshot_00400.hybrid snapshot_00400.hybrid.v2 >&1

echo "comparing serial algo1 and serial algo 2"
cmp snapshot_00400.serial snapshot_00400.serial.algo2 >&1

rm snapshot_*

echo "checking coherency... ordered evolution"

echo "snark loop..."

./conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.serial

./algo_2/algo2_conway_serial_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.serial.algo2

./conway_omp_pgm.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.omp

./conway_omp_pgm.v2.out -r -f snark_loop_01.pgm -n 100 -s 100 -e 0
mv snapshot_00100 snapshot_00100.omp.v2

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
echo "comparing omp v1 and v2"
cmp snapshot_00100.omp snapshot_00100.omp.v2 >&1

echo "comparing serial algo1 and serial algo 2"
cmp snapshot_00100.serial snapshot_00100.serial.algo2 >&1

rm snapshot_*

echo "glider gun..."

./conway_serial_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 0
mv snapshot_00400 snapshot_00400.serial

./algo_2/algo2_conway_serial_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 0
mv snapshot_00400 snapshot_00400.serial.algo2

./conway_omp_pgm.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 0
mv snapshot_00400 snapshot_00400.omp

./conway_omp_pgm.v2.out -r -f gosper_glider_gun_01.pgm -n 400 -s 400 -e 0
mv snapshot_00400 snapshot_00400.omp.v2

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
echo "comparing omp v1 and v2"
cmp snapshot_00400.omp snapshot_00400.omp.v2 >&1

echo "comparing serial algo1 and serial algo 2"
cmp snapshot_00400.serial snapshot_00400.serial.algo2 >&1

rm snapshot_*

cd algo_2 

echo "running check for algo2"
bash check.sh

cd ..
make clean
