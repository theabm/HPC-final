# Exercise 2 - Benchmarking MKL, OpenBLAS, and BLIS

Exercise 2 consists of benchmarking the level 3 math function `gemm` for matrix-matrix 
multiplication for three HPC math libraries, namely, MKL, OpenBLAS, and BLIS. 

We compare the obtained GFLOPS with the expected peak performance for the two types 
of architectures in Orfeo in two cases: 
1. Fixed cores: We fix the number of cores to use the full node. For THIN nodes, 
this corresponds to 24, while for EPYC nodes this corresponds to 128 cores. Then, 
using all three libraries, we compute the GFLOPS by running `gemm` using 
squared matrices of increasing sizes, from $2000$ to $20000$, with steps of 
$2000$.
2. Fixed size: We fix the dimensions of the squared matrices to $10000$ and we 
slowly increase the number of cores used for the matrix matrix multiplication. 
We start from 1 and end at 24 for THIN and 128 for EPYC. 

We repeat these measures for both double and single precision operations, and for 
THIN end EPYC nodes.

In both cases, we compare the GFLOPS obtained with the theoretical peak performance 
of the resources we use.

# Structure of this directory: 
```
1-fixed_cores/all the data and ipynb to generate graphs for this section.
2-fixed_size/all the data and ipynb to generate graphs for this section.
bins/all the binaries for THIN and EPYC, for all three libraries.
scripts/all scripts used to launch the jobs.
dgemm.c - the `gemm` program we used.
Makefile -the Makefile we used to compile `dgemm.c`
README.md - This README
```
# To compile the code 
Follow instructions [here](https://github.com/Foundations-of-HPC/Foundations_of_HPC_2022/blob/main/Assignment/README.MD) 
and [here](https://github.com/Foundations-of-HPC/Foundations_of_HPC_2022/blob/main/Assignment/exercise2/README.md) to compile 
the BLIS library. The MKL and OpenBLAS libraries were already installed in Orfeo.

# To run the code

To run the code with N OMP threads using matrices of size K, select the 
appropriate executable and run: 
```
$export OMP_NUM_THREADS=N 
$./selected_executable K K K 
```
