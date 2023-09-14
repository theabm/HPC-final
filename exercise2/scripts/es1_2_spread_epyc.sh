#!/bin/bash
#SBATCH --job-name=es1_2_spread_epyc
#SBATCH --output=es1_2_spread_epyc.out
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=EPYC
date;hostname;pwd
module load mkl/latest
module load openBLAS/0.3.23-omp
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export LD_LIBRARY_PATH=/u/dssc/s271711/blis_epyc/lib:$LD_LIBRARY_PATH

TESTTYPE=epyc_spread
cd /u/dssc/s271711/fast/exercise2/2-fixed_size
for DATATYPE in float double; do
OUTFILE_MKL=${TESTTYPE}/${DATATYPE}/results_mkl.csv
OUTFILE_OBLAS=${TESTTYPE}/${DATATYPE}/results_oblas.csv
OUTFILE_BLIS=${TESTTYPE}/${DATATYPE}/results_blis.csv
mkdir -p ${TESTTYPE}/${DATATYPE}
echo "cores,size,size,size,time,gflops" >> "${OUTFILE_MKL}"
echo "cores,size,size,size,time,gflops" >> "${OUTFILE_OBLAS}"
echo "cores,size,size,size,time,gflops" >> "${OUTFILE_BLIS}"
for j in 1 10 20 30 40 50 60 70 80 90 100 110 120 128; do
export OMP_NUM_THREADS=$j
export BLIS_NUM_THREADS=$j
for i in {0..9} ; do
	printf "${OMP_NUM_THREADS}," >> "${OUTFILE_MKL}"
	srun ./bins/epyc/dgemm_mkl_${DATATYPE}.x 10000 10000 10000 >> "${OUTFILE_MKL}"
	printf "${OMP_NUM_THREADS}," >> "${OUTFILE_OBLAS}"
	srun ./bins/epyc/dgemm_oblas_${DATATYPE}.x 10000 10000 10000 >> "${OUTFILE_OBLAS}"
	printf "${OMP_NUM_THREADS}," >> "${OUTFILE_BLIS}"
	srun ./bins/epyc/dgemm_blis_${DATATYPE}.x 10000 10000 10000 >> "${OUTFILE_BLIS}"
done
done
done
