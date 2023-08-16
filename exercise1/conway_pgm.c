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
                fname = (char*)malloc(sizeof(optarg)+sizeof(f_prefix)+1 );
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
void save_grid(char * fname, MPI_Comm comm, int rank, char * header, int header_size, MPI_Offset offset, char * data, int my_rows, int cols)
{

    MPI_File fh;

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

void upgrade_cell(char ** grid_prev, char** grid, int i, int j)
{

    // by construction, the index for rows will work as is, since we take 
    // care of halo cells explicitly.
    // However, we need to deal with the fact that a column may be out of bounds
    int jm1 = (j-1)%(int)cols + cols*(j-1<0);
    int jp1 = (j+1)%(int)cols + cols*(j+1<0);

    char n_alive_cells = grid_prev[i-1][jm1]
        + grid_prev[i-1][j]
        + grid_prev[i-1][jp1]
        + grid_prev[i][jm1]
        + grid_prev[i][jp1]
        + grid_prev[i+1][jm1]
        + grid_prev[i+1][j]
        + grid_prev[i+1][jp1];
    
    if (n_alive_cells >= 2 && n_alive_cells <=3){
        grid[i][j] = ALIVE;
    }
    else{
        grid[i][j] = DEAD;
    }
}

int main(int argc, char **argv){
    // initialize MPI
    MPI_Init(&argc, &argv);

    // each rank gets info about itself
    int rank;
    int size;

    char ** grid, ** grid_prev;
    char * data, * data_prev;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // get arguments
    // MPI standard does not specify how arguments are processed. 
    // Ideally it would be better to have master process the args and then 
    // broadcast the arguments. However, for now this is ok.
    get_args(argc, argv);

    // setting up neighbors for 1D splitting
    // Note that with 1D splitting, we have a limitation on how many processes 
    // this program will work with. 
    // If the processes are more than the number of rows, then the program will 
    // not use some of these processes.
    const int prev = (size - 1)*(rank==0) + (rank - 1)*!(rank==0);
    const int next = 0*(rank==(size-1)) + (rank + 1)*!(rank==(size-1));

    if(action == INIT){

        // for each rank, we want to determine how many rows it will need to 
        // handle.
        // Since we are running the program with -i option (INIT), 
        // we also specify the number of rows and columns. 
        // So this operation is well defined.
        
        const int my_rows = get_my_rows(rows, rank, size);
        const int my_row_offset = get_my_row_offset(rows, rank, size);
        const int augmented_rows = my_rows + 2;

        const MPI_Offset my_file_offset = my_row_offset * cols * sizeof(char);

        grid = (char **) malloc(augmented_rows * sizeof(char *));
        grid_prev = (char **) malloc(augmented_rows * sizeof(char *));

        data = (char *) malloc( augmented_rows * cols * sizeof(char));
        data_prev = (char *) malloc( augmented_rows * cols * sizeof(char));

        if(
                grid == NULL
                || grid_prev == NULL
                || data == NULL
                || data_prev == NULL
                )
        {
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);

        }

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
                grid_prev[i][j] = drand48() > 0.5 ? ALIVE : DEAD ;
                // grid[i][j] = ALIVE;
        }

        save_grid(fname, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data_prev, my_rows, cols);

    }
    else if (action == RUN){

        // read file 
        // in general, we cannot assume we will use the same number 
        // of processors as the initialization phase
        // so we need to first read the header, broadcast to the other 
        // processes, distribute an area of the matrix for each process 
        // and read from the file
        
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
            fclose(fh_posix);
        }

        MPI_Bcast(opt_args, 2, MPI_INT, 0, MPI_COMM_WORLD);

        rows = opt_args[0];
        cols = opt_args[1];

        const int my_rows = get_my_rows(rows, rank, size);
        const int my_row_offset = get_my_row_offset(rows, rank, size);
        const int augmented_rows = my_rows + 2;

        const MPI_Offset my_file_offset = my_row_offset * cols * sizeof(char);

        grid = (char **) malloc(augmented_rows * sizeof(char *));
        grid_prev = (char **) malloc(augmented_rows * sizeof(char *));

        data = (char *) malloc( augmented_rows * cols * sizeof(char));
        data_prev = (char *) malloc( augmented_rows * cols * sizeof(char));

        if(
                grid == NULL
                || grid_prev == NULL
                || data == NULL
                || data_prev == NULL
                )
        {
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);

        }
        
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

        MPI_File fh;

        const int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

        if(err != MPI_SUCCESS){
            fprintf(stderr, "Error opening %s\n", fname);
            MPI_Abort(MPI_COMM_WORLD, err);
        }

        MPI_File_read_at_all(fh, my_total_file_offset, data_prev + cols, my_rows*cols, MPI_CHAR, MPI_STATUS_IGNORE);

        MPI_File_close(&fh);

        // At this point, all the processes have read their portion of the 
        // matrix. The halo regions have been set to zero, and we are ready 
        // to start the game of life.
        
        // file name (it will always be the same length so we only need it once)
        char * file_name = malloc(snprintf(NULL, 0, "snapshot_%05d.pgm", 0)+1);
        if(file_name == NULL){
            printf("Not enough space.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }

        MPI_Request prev_send_request, next_send_request;
        MPI_Request prev_recv_request, next_recv_request;

        const int prev_tag = 0; 
        const int next_tag = 1;

        char *tmp = NULL;

        for(int t = 0; t < n; ++t){

            // no longer needed since it is taken cafe of by MPI_Wait 
            // and MPI_Request_free().
            // prev_send_request = MPI_REQUEST_NULL;
            // next_send_request = MPI_REQUEST_NULL;
            //
            // prev_recv_request = MPI_REQUEST_NULL;
            // next_recv_request = MPI_REQUEST_NULL;

            // Step 1. Start non blocking exchange of halo cells.

            // we post send and recv requests for two halo regions.
            // matrix: 
            // ---------------------  row 0 (HALO)
            // 010100011001100100100  row 1
            // 010010000100101001111  row 2             --
            // ...                                        | these rows dont need 
            // 001010111110011000001  row my_rows - 1   __| halo regions
            // 010010010000100101000  row my_rows
            // ---------------------  row my_rows + 1 (HALO)
            
            // send row: 1 to the previous rank
            MPI_Isend(data_prev + cols, cols, MPI_CHAR, prev, prev_tag, MPI_COMM_WORLD, &prev_send_request);
            // send row: my_rows to next rank
            MPI_Isend(data_prev + my_rows*cols, cols, MPI_CHAR, next, next_tag, MPI_COMM_WORLD, &next_send_request);

            // receive from prev and put in row 0 (halo region)
            MPI_Irecv(data_prev, cols, MPI_CHAR, prev, next_tag, MPI_COMM_WORLD, &prev_recv_request);
            // receive from next and put in row my_rows + 1 (halo region)
            MPI_Irecv(data_prev + cols*my_rows + cols, cols, MPI_CHAR, next, prev_tag, MPI_COMM_WORLD, &next_recv_request);

            // Once the send and receive have completed, each process should 
            // have the halo regions. So we can update the entire grid. 

            // However, in the meantime, we can process all of the internal rows 
            // which don't need halo regions. 
            // This allows us to hide the latency which occurs in messagge 
            // passing.
            // These rows that dont need the halo regions are 
            // rows 2 and my_rows - 1 (look at diagram above)

            // Step 2. Process internal cells to hide latency

            // POSSIBILITY TO PARARELLIZE WITH OMP
            for(int row = 2; row < my_rows; ++row){
                for(int col = 0; col < cols; ++cols){
                    // watch for access pattern of memory. 
                    // could be necessary to optimize and access differently
                    upgrade_cell(grid_prev, grid, row, col);
                }
            }

            // At this point, all internal cells have been processed
            // and we can check if the recv operation has completed
            
            // Step 3. Check that recv has been completed.
            // (theoretically, we could also use MPI_Test to see if its 
            // completed and if it isnt, do some other things, but all that's 
            // really left is the computation of the borders)
            
            // look at https://stackoverflow.com/questions/10882581/mpi-isend-request-parameter 
            // and 
            // https://stackoverflow.com/questions/22410827/mpi-reuse-mpi-request
            // for more information on why this is necessary.
            // However, the main idea is that these operations free the request 
            // handles and set them to MPI_REQUEST_NULL which can then be reused.
            MPI_Wait(&prev_recv_request, MPI_STATUS_IGNORE);
            MPI_Wait(&next_recv_request, MPI_STATUS_IGNORE);
            MPI_Request_free(&prev_send_request);
            MPI_Request_free(&next_send_request);

            // At this point, the halo has been received by all processes.
            // So we need to update the 2 rows that require halo regions.
            // These are row 1 and row my_rows

            // Step 4. Update limiting rows (row 1 and row my_rows)
            for(int col=0; col<cols; ++col){
                upgrade_cell(grid_prev, grid, 1, col);
                upgrade_cell(grid_prev, grid, my_rows, col);
            }

            // At this point, data_prev contains the data at t-1, 
            // while data contains the data for time t.
            // Now we need to check if we need to save

            // Step 5. Check if need to save, and if we do, save grid to pgm
            
            if(t%s == 0){
                sprintf(file_name, "snapshot_%05d", t);
                save_grid(file_name, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data, my_rows, cols);
            }

            // Step 6. Swap data_prev and data 
            tmp = data;
            data = data_prev;
            data_prev = tmp;
            // it is good practice to avoid dangling pointers.
            tmp = NULL;
        }
    
        free(file_name);
    }
    else{
        printf("Unknown action. Abort");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    free(grid);
    free(grid_prev);
    free(data);
    free(data_prev);

    if ( fname != NULL )
      free ( fname );

    MPI_Finalize();
    return 0;
}

