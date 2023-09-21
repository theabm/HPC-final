# Foundations of High Performance Computing - Final Exam 
## DSSC - A.A 2022-2023 

This repostory contains the implementation, data, scripts, and report 
for the final assignment for the course of High Performance Computing. 

The project consists of: 
- Exercise 1: Implementing a scalable, hybrid MPI + OMP version of Conway's game
of life and analyzing the scalability.
- Exercise 2: Benchmarking the three high performance math libraries: MKL, 
OpenBLAS, and BLIS.

The full problem description can be found in the [repository](https://github.com/Foundations-of-HPC/Foundations_of_HPC_2022/tree/main/Assignment) of the course.

All the problems were analyzed using the Orfeo super computer using: 
- "THIN" nodes: 24 cores of Intel Xeon Gold 6126 @ 2.6 GHz.
- "EPYC" nodes: 128 cores of AMD EPYC

# Structure of this github directory: 
```
README.md
report.pdf
report/all the files used for the report 
exercise1/README.md
exercise1/all the files used for exercise 1
exercise2/README.md
exercise2/all the files used for exercise 2
```

Where: 
- README.md briefly describes what was done and how to navigate the directory.
- report.pdf contains a detailed report of all the exercises. 
- exercise1/README.md contains all the necessary details for exercise 1 
- exercise2/README.md contains all the necessary details for exercise 2 


# Exercise 1
Exercise 1 consists in the implementation of a scalable, hybrid MPI + OMP 
version of Conway's game of life. Furthermore, an analysis of the strong OMP, 
strong MPI, and weak MPI scalability will be conducted in the report.

# Exercise 2 
Exercise 2 consists of benchmarking the level 3 math function `gemm` for matrix-matrix 
multiplication for three HPC math libraries, namely, MKL, OpenBLAS, and BLIS. 
We compare the obtained GFLOPS with the expected peak performance for the two types 
of architectures in Orfeo. 


