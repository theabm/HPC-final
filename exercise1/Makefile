all : conway_serial_pgm.out conway_omp_pgm.out conway_mpi_pgm.out

conway_serial_pgm.out : conway_serial_pgm.c
	gcc -S -O3 -march=native -Wall -Wextra $^ -o conway_serial_pgm.s
	gcc -O3 -march=native -Wall -Werror -Wextra $^ -o $@

conway_omp_pgm.out : conway_omp_pgm.c
	gcc -fopenmp -O3 -march=native -Wall -Werror -Wextra $^ -o $@

conway_mpi_pgm.out : conway_mpi_pgm.c
	mpicc -O3 -march=native -Wall -Werror -Wextra $^ -o $@

clean : 
	rm *.out
	rm snapshot* a1.mpg
#
# sr : 
# 	cp pattern_pgm/gosper_glider_gun.pgm .
# 	./conway_serial_pgm.out -r -f gosper_glider_gun -n 400 -s 1 
# 	ffmpeg -f image2 -i snapshot_%05d.pgm ./a1.mpg
# 	vlc a1.mpg
# 	rm snapshot* gosper_glider_gun.pgm
# 	rm a1.mpg
#
# mr : 
# 	cp pattern_pgm/gosper_glider_gun.pgm .
# 	mpirun -np 5 conway_serial_pgm.out -r -f gosper_glider_gun -n 400 -s 1 
# 	ffmpeg -f image2 -i snapshot_%05d.pgm ./a1.mpg
# 	vlc a1.mpg
# 	rm snapshot* gosper_glider_gun.pgm
# 	rm a1.mpg
#
#

