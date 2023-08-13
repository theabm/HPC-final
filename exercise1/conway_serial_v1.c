#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <mpi.h>

#define INIT 1
#define RUN  2

#define DEAD 0 
#define ALIVE 1

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1


// DEFAULT VALUES
//
char fname_deflt[] = "game_of_life.pgm";

int   action = 0;
int   k      = K_DFLT;
int   rows   = 0;
int   cols   = 0;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;


void get_args( int argc, char **argv )
{
    char *optstring = "irk:e:f:n:s:";

    int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {

            case 'i':
                action = INIT; break;

            case 'r':
                action = RUN; break;

            case 'k':
                k = atoi(optarg); rows = k; cols = k; break;

            case 'e':
                e = atoi(optarg); break;

            case 'f':
                fname = (char*)malloc( sizeof(optarg)+1 );
                sprintf(fname, "%s", optarg );
                break;

            case 'n':
                n = atoi(optarg); break;

            case 's':
                s = atoi(optarg); break;

            default :
                printf("argument -%c not known\n", c ); break;

        }
    }
}

void display_args(int rank, int size){
    printf("I am rank %d of %d.\naction (i : 1\tr : 2) -- %d\nk (size) -- %d\ne (0 : ORDERED\t1 : STATIC) -- %d\nf (filename) -- %s\nn (steps) -- %d\ns (save frequency) -- %d\n", rank, size, action, k, e, fname, n, s );
}

int get_my_rows(int total_rows, int rank, int size)
{
    // returns the number of rows that each process needs to handle

    // We divide the total rows by the size. In most cases this will not 
    // be a perfect division, so we must handle the remainder R. Since 
    // R < size, we simply give an extra row to all of the ranks that are 
    // less than R.
    // So, if we have 8 ranks, and the remainder is 3, only the first 3 ranks 
    // get an extra row to deal with.
    //
    int my_rows = total_rows/size + 1*(rank<(total_rows%size));

    return my_rows;
}

int get_my_offset(int total_rows, int rank, int size)
{
    // returns the offset from which the rank must get its rows 

    // for example: for rank 0, the offset will be zero as it will read 
    // my_rows(of rank 0) from the start. However, rank 1 will have to read 
    // after the rows that rank 0 handles (accounting for the remaineder), 
    // and so on.

    // nrows is the interger division between total_rows and size
    int nrows = total_rows/size;
    int remainder = total_rows%size;

    // The idea is the following, if the division was pefect, we would do 
    // nrows * rank. However, we still need to account for the remainder.
    // If the rank is one of the processors who got an extra row, we simply 
    // need to add rank to the offset. If we are dealing with a processor 
    // who is greater than the remainder, then we simply need to add the 
    // whole remainder.
    
    int my_offset = nrows*rank
        + rank*(rank < remainder) 
        + remainder*(rank > remainder);

    return my_offset;
}

int main(int argc, char **argv){


    // initialize MPI
    MPI_Init(&argc, &argv);

    // each rank gets info about itself
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // get arguments
    // MPI standard does not specify how arguments are processed. 
    // Ideally it would be better to have master process and then broad cast 
    // the arguments. However, for now this is ok.
    get_args(argc, argv);

    // setting up neighbors
    int prev = (size - 1)*(rank==0) + (rank - 1)*!(rank==0);
    int next = 0*(rank==(size-1)) + (rank + 1)*!(rank==(size-1));

    // for each rank, determine its part of the matrix
    int my_rows = get_my_rows(rows, rank, size);
    int my_offset = get_my_offset(rows, rank, size);

    // printf("I am rank %d of %d and I need to deal with %d rows out of %d. My offset is %d\n", rank, size, my_rows, rows, my_offset);

    char ** grid, ** grid_prev;
    char * data, * data_prev;

    // allocate memory for the data each rank will store
    // grid is an array of pointers which will be connected to the 1D structure 
    // that actually holds the data and will allow us to access this data as if 
    // it was a 2D matrix. This is a common trick taught in class to 
    // avoid allocating memory inefficiently. 
    
    grid = (char **) malloc((my_rows + 2) * sizeof(char *));
    grid_prev = (char **) malloc((my_rows + 2) * sizeof(char *));

    data = (char *) malloc( (my_rows + 2) * (cols + 2) * sizeof(char));
    data_prev = (char *) malloc( (my_rows + 2) * (cols + 2) * sizeof(char));

    // now we need to make sure that each of the pointers in grid point to the 
    // right place of the data.
    int augmented_cols = cols+2;

    for(int i = 0; i<my_rows+2; ++i){
        grid[i] = data + i*augmented_cols;
        grid_prev[i] = data_prev + i*augmented_cols;
    }

    // initialize the halo regions to being DEAD
    for(int j = 0; j<cols+2; ++j){
        grid[0][j] = grid[my_rows + 1][j] = grid_prev[0][j]
            = grid_prev[my_rows + 1][j] = DEAD;
    }
    for(int j = 0; j<my_rows+2; ++j){
        grid[j][0] = grid[j][cols + 1] = grid_prev[j][0]
            = grid_prev[j][cols + 1] = DEAD;
    }

    if(action == INIT){
        srand48(1*rank);
        // initialize the matrix randomly
        // possibility to optimize loop here by exploiting 
        // ILP
        for(int i = 0; i<my_rows; ++i){
            for(int j = 0; j<cols; ++j)
                grid[i][j] = drand48() > 0.5 ? ALIVE : DEAD ;
        }

        // MPI_File fh;
        // int err;
        // MPI_Info info = MPI_NULL;
        //
        // err = MPI_File_Open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, info, &fh);
        // if(err != MPI_SUCCESS){
        //     fprintf(stderr, "Error opening %s\n", fname);
        //     return err;
        // }
        //
        // err = MPI_File_write_at_all(fh, my_offset,)

    }




    



    // display_args(rank, size);

    if ( fname != NULL )
      free ( fname );

    MPI_Finalize();
    return 0;
}

