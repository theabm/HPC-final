#!/bin/bash
#SBATCH --job-name=es1_1
#SBATCH --output=es1_1.out
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=2:00:00
#SBATCH --partition=THIN
date;hostname;pwd
module load mkl/latest
module load openBLAS/0.3.23-omp
export OMP_NUM_THREADS=24
export BLIS_NUM_THREADS=24
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/u/dssc/s271711/blis/lib:$LD_LIBRARY_PATH

TESTTYPE=thin
cd /u/dssc/s271711/fast/exercise2/1-fixed_cores
for DATATYPE in float double; do
OUTFILE_MKL=${TESTTYPE}/${DATATYPE}/results_mkl.csv
OUTFILE_OBLAS=${TESTTYPE}/${DATATYPE}/results_oblas.csv
OUTFILE_BLIS=${TESTTYPE}/${DATATYPE}/results_blis.csv
mkdir -p ${TESTTYPE}/${DATATYPE}
echo "size,size,size,time,gflops" > "${OUTFILE_MKL}"
echo "size,size,size,time,gflops" > "${OUTFILE_OBLAS}"
echo "size,size,size,time,gflops" > "${OUTFILE_BLIS}"
for j in {2000..20000..2000}; do
for i in {0..9} ; do
	srun ./bins/${TESTTYPE}/dgemm_mkl_${DATATYPE}.x $j $j $j >> "${OUTFILE_MKL}"
	srun ./bins/${TESTTYPE}/dgemm_oblas_${DATATYPE}.x $j $j $j >> "${OUTFILE_OBLAS}"
	srun ./bins/${TESTTYPE}/dgemm_blis_${DATATYPE}.x $j $j $j >> "${OUTFILE_BLIS}"
done
done
done
