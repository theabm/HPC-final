all : algo2_conway_serial_pgm.out algo2_conway_omp_pgm.out algo2_conway_mpi_pgm.out

algo2_conway_serial_pgm.out : algo2_conway_serial_pgm.c
	gcc -O3 -march=native -Wall -Werror -Wextra $^ -o $@

algo2_conway_omp_pgm.out : algo2_conway_omp_pgm.c
	gcc -fopenmp -O3 -march=native -Wall -Werror -Wextra $^ -o $@

algo2_conway_mpi_pgm.out : algo2_conway_mpi_pgm.c
	mpicc -O3 -march=native -Wall -Werror -Wextra $^ -o $@

clean : 
	rm *.out
	rm snapshot_* a1.*

