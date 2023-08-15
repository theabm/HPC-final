// THIS IS VERSION V1 WHICH USES THE PGM FILE FORMAT. THIS MEANS THAT 
// EVERYWHERE WE ARE USING MATRICES OF CHARS AND ENCODING DEAD OR ALIVE 
// AS 0 AND 255 RESPECTIVELY. SOME KEY DIFFERENCES ARE:
// THE HEADER USES P5 INSTEAD OF P4 FOR PBM.
// THE HEADER INCLUDES AN EXTRA MAX VAL PARAMETER. 
// IT IS NOT NECESSARY TO ENCODE EVERYTHING AT BIT LEVEL.


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <mpi.h>

#define INIT 1
#define RUN  2

// for pgm files, 0 is black and 255(MAXVAL) is white
// since we want black squares to denote alive, we define 
// alive to be 0 (black) and dead to be 255 (white)
#define DEAD 255
#define ALIVE 0
#define MAX_VAL 255

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1

#define HEADER_FORMAT_STRING "P5 %d %d %d\n"


// DEFAULT VALUES
char fname_deflt[] = "game_of_life.pgm";

int   action = 0;
int   k      = K_DFLT;
int   rows   = 0;
int   cols   = 0;
int   e      = ORDERED;
int   n      = 10000;
int   s      = 1;
char *fname  = NULL;
char *f_prefix = ".pgm";


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
                fname = (char*)malloc( sizeof(optarg)+sizeof(f_prefix)+1 );
                sprintf(fname, "%s%s", optarg, f_prefix);
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

int get_my_row_offset(int total_rows, int rank, int size)
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
void save_grid(char * fname, MPI_Comm comm, int rank, char * header, int header_size, MPI_Offset offset, char * data, int my_rows, int cols, MPI_File fh){

    // Opening the file in MPI_MODE_WRITEONLY or MPI_MODE_CREATE_ONLY
    const int err = MPI_File_open(comm, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    // Check that the file was opened correctly
    if(err != MPI_SUCCESS){
        fprintf(stderr, "Error opening %s\n", fname);
        MPI_Abort(MPI_COMM_WORLD, err);
    }
    
    if(rank == 0){
        MPI_File_write_at(fh, 0, header, header_size, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    MPI_File_write_at_all(fh, offset, data + cols, my_rows*cols, MPI_CHAR, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
}

int main(int argc, char **argv){
    // initialize MPI
    MPI_Init(&argc, &argv);

    // each rank gets info about itself
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char ** grid, ** grid_prev;
    char * data, * data_prev;

    // get arguments
    // MPI standard does not specify how arguments are processed. 
    // Ideally it would be better to have master process and then broad cast 
    // the arguments. However, for now this is ok.
    get_args(argc, argv);

    // setting up neighbors
    const int prev = (size - 1)*(rank==0) + (rank - 1)*!(rank==0);
    const int next = 0*(rank==(size-1)) + (rank + 1)*!(rank==(size-1));


    if(action == INIT){

        // for each rank, determine its part of the matrix
        const int my_rows = get_my_rows(rows, rank, size);
        const int my_row_offset = get_my_row_offset(rows, rank, size);
        const MPI_Offset my_file_offset = my_row_offset * cols * sizeof(char);

        printf("I am rank %d of %d and I need to deal with %d rows out of %d. My offset is %d\n", rank, size, my_rows, rows, my_row_offset);

        const int augmented_rows = my_rows + 2;

        grid = (char **) malloc(augmented_rows * sizeof(char *));
        grid_prev = (char **) malloc(augmented_rows * sizeof(char *));

        data = (char *) malloc( augmented_rows * cols * sizeof(char));
        data_prev = (char *) malloc( augmented_rows * cols * sizeof(char));

        // now we need to make sure that each of the pointers in grid point to the 
        // right place of the data.
        for(int i = 0; i<augmented_rows; ++i){
            grid[i] = data + i*cols;
            grid_prev[i] = data_prev + i*cols;
        }

        // initialize the halo regions to being DEAD
        for(int j = 0; j<cols; ++j){
            grid[0][j] = grid[my_rows + 1][j] = grid_prev[0][j]
                = grid_prev[my_rows + 1][j] = DEAD;
        }

        MPI_File fh;

        // HEADER INFO CALCULATION

        // intuition: https://stackoverflow.com/questions/3919995/determining-sprintf-buffer-size-whats-the-standard
        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);
        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        const MPI_Offset header_offset = header_size * sizeof(char);
        const MPI_Offset my_total_file_offset = my_file_offset + header_offset;


        // INITIALIZE MATRIX RANDOMLY
        
        // for reproducible results. If we want truly random then we also 
        // need to include some changing information such as time
        srand48(1*rank); 

        // possibility to optimize loop here by exploiting 
        // ILP
        // HERE IS A POSSIBILITY TO USE OPENMP
        for(int i = 1; i<my_rows+1; ++i){
            for(int j = 0; j<cols; ++j)
                grid[i][j] = drand48() > 0.5 ? ALIVE : DEAD ;
                // grid[i][j] = ALIVE;
        }

        save_grid(fname, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data, my_rows, cols, fh);
        // // Opening the file in MPI_MODE_WRITEONLY or MPI_MODE_CREATE_ONLY
        // const int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        //
        // // Check that the file was opened correctly
        // if(err != MPI_SUCCESS){
        //     fprintf(stderr, "Error opening %s\n", fname);
        //     MPI_Abort(MPI_COMM_WORLD, err);
        // }
        // 
        // if(rank == 0){
        //     MPI_File_write_at(fh, 0, header, header_size, MPI_CHAR, MPI_STATUS_IGNORE);
        // }
        //
        // MPI_File_write_at_all(fh, my_total_file_offset, data + cols, my_rows*cols, MPI_CHAR, MPI_STATUS_IGNORE);
        //
        // MPI_File_close(&fh);

    }
    else if (action == RUN){

        MPI_File fh;

        // read file 
        // in general, we cannot assume we will use the same number 
        // of processors as the initialization phase
        // so we need to first read the header, broadcast to the other 
        // processes, distribute an area of the matrix for each process 
        // and read from the file
        
        const int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

        if(err != MPI_SUCCESS){
            fprintf(stderr, "Error opening %s\n", fname);
            MPI_Abort(MPI_COMM_WORLD, err);
        }

        int opt_args[2] = {0,0};

        if(rank == 0){

            FILE * fh_posix = fopen(fname, "r");

            // we know that the magic number is P5 so we set the offset as the  
            // length of P5 * sizeof(char)
            if(fh_posix == NULL){
                fprintf(stderr, "Error opening %s\n", fname);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);

            }

            fscanf(fh_posix, "P5 %d %d 255\n",opt_args, opt_args+1 );
        }

        MPI_Bcast(opt_args, 2, MPI_INT, 0, MPI_COMM_WORLD);

        rows = opt_args[0];
        cols = opt_args[1];

        // printf("I am rank %d and I got %d rows and %d cols\n", rank, rows, cols);
        
        const int my_rows = get_my_rows(rows, rank, size);
        const int my_row_offset = get_my_row_offset(rows, rank, size);
        const MPI_Offset my_file_offset = my_row_offset * cols * sizeof(char);

        const int augmented_rows = my_rows + 2;

        grid = (char **) malloc(augmented_rows * sizeof(char *));
        grid_prev = (char **) malloc(augmented_rows * sizeof(char *));

        data = (char *) malloc( augmented_rows * cols * sizeof(char));
        data_prev = (char *) malloc( augmented_rows * cols * sizeof(char));

        // now we need to make sure that each of the pointers in grid point to the 
        // right place of the data.
        for(int i = 0; i<augmented_rows; ++i){
            grid[i] = data + i*cols;
            grid_prev[i] = data_prev + i*cols;
        }

        // initialize the halo regions to being DEAD
        for(int j = 0; j<cols; ++j){
            grid[0][j] = grid[my_rows + 1][j] = grid_prev[0][j]
                = grid_prev[my_rows + 1][j] = DEAD;
        }

        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);
        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        const MPI_Offset header_offset = header_size * sizeof(char);
        const MPI_Offset my_total_file_offset = my_file_offset + header_offset;

        MPI_File_read_at_all(fh, my_total_file_offset, data + cols, my_rows*cols, MPI_CHAR, MPI_STATUS_IGNORE);

        MPI_File_close(&fh);

        char * fname2 = "lala2.pgm";

        save_grid(fname2, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data, my_rows, cols, fh);

    
    }


    // display_args(rank, size);

    if ( fname != NULL )
      free ( fname );

    MPI_Finalize();
    return 0;
}

