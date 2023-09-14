all : conway_serial_pgm.out conway_omp_pgm.out conway_mpi_pgm.out conway_hybrid_pgm.out conway_omp_pgm.v2.out

conway_serial_pgm.out : conway_serial_pgm.c
	gcc -O3 -march=native -Wall -Werror -Wextra $^ -o $@

conway_omp_pgm.out : conway_omp_pgm.c
	gcc -fopenmp -O3 -march=native -Wall -Werror -Wextra $^ -o $@

conway_omp_pgm.v2.out : conway_omp_pgm.v2.c
	gcc -fopenmp -O3 -march=native -Wall -Werror -Wextra $^ -o $@

conway_mpi_pgm.out : conway_mpi_pgm.c
	mpicc -O3 -march=native -Wall -Werror -Wextra $^ -o $@

conway_hybrid_pgm.out : conway_hybrid_pgm.c
	mpicc -fopenmp -O3 -march=native -Wall -Werror -Wextra $^ -o $@

clean : 
	rm *.out

clean_output : 
	rm snapshot_* a1.*

