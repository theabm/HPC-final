#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <mpi.h>
#include <unistd.h>
#include <omp.h>

#define INIT 1
#define RUN  2

#define DATA(i,j) (data[(i)*cols + (j)])
#define DATA_PREV(i,j) (data_prev[(i)*cols + (j)])

// for pgm files, 0 is black and MAXVAL is white.
// So dead is 0 and alive is 1
#define DEAD 0
#define ALIVE 1
#define MAX_VAL 1

#define CACHE_LINE_SIZE 64

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1

#define HEADER_FORMAT_STRING "P5 %d %d %d\n"


// DEFAULT VALUES
char fname_deflt[] = "game_of_life.pgm";

int   action = 0;
int   k      = K_DFLT;
int   rows   = K_DFLT;
int   cols   = K_DFLT;
int   e      = STATIC;
int   n      = 10000;
int   s      = 0;
char *fname  = NULL;

void get_args( int argc, char **argv )
{
    char *optstring = "irk:e:f:n:s:";

    int c;
    while ((c = getopt(argc, argv, optstring)) != -1)
    {
        switch(c) 
        {

            case 'i':
                action = INIT; break;

            case 'r':
                action = RUN; break;

            case 'k':
                k = atoi(optarg); rows = cols = k; break;

            case 'e':
                e = atoi(optarg); break;

            case 'f':
                size_t str_len = strlen(optarg) + 1;
                fname = (char*)malloc(str_len * sizeof(char));
                sprintf(fname, "%s", optarg);
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

void display_args(int rank, int size)
{
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
        + rank*(rank <= remainder) 
        + remainder*(rank > remainder);

    return my_offset;
}
void save_grid(char * restrict fname, MPI_Comm comm, int rank, char * restrict header, int header_size, MPI_Offset offset, unsigned char * restrict data, int my_rows, int cols)
{

    MPI_File fh;

    // Opening the file in MPI_MODE_WRITEONLY or MPI_MODE_CREATE_ONLY
    const int err = MPI_File_open(comm, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    // Check that the file was opened correctly
    if(err != MPI_SUCCESS)
    {
        fprintf(stderr, "Error opening %s\n", fname);
        MPI_Abort(MPI_COMM_WORLD, err);
    }
    
    if(rank == 0)
    {
        MPI_File_write_at(fh, 0, header, header_size, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    MPI_File_write_at_all(fh, offset, data + cols, my_rows*cols, MPI_CHAR, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
}

void upgrade_cell_static(unsigned char * restrict data_prev, unsigned char * restrict data, int i, int j)
{
    DATA(i,j) = DEAD;

    register int jm1 = j==0 ? cols-1 : j-1;
    register int jp1 = j==(cols-1) ? 0 : j+1;
    register int im1 = i-1;
    register int ip1 = i+1;
    
    unsigned char tmp0=DATA_PREV(im1, jm1);
    unsigned char tmp3=DATA_PREV(i,jm1);
    unsigned char tmp5=DATA_PREV(ip1,jm1);

    unsigned char tmp1=DATA_PREV(im1,j);
    unsigned char tmp2=DATA_PREV(im1,jp1);

    unsigned char tmp4=DATA_PREV(i,jp1);

    unsigned char tmp6=DATA_PREV(ip1,j);
    unsigned char tmp7=DATA_PREV(ip1,jp1);

    unsigned char n_alive_cells = tmp0+tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7;

    DATA(i,j) = ALIVE*(n_alive_cells==3) + DATA_PREV(i,j)*(n_alive_cells==2);
}

void upgrade_cell_ordered(unsigned char * data, int i, int j)
{
    register int jm1 = j==0 ? cols-1 : j-1;
    register int jp1 = j==(cols-1) ? 0 : j+1;
    register int im1 = i-1;
    register int ip1 = i+1;

    unsigned char tmp0=DATA(im1, jm1);
    unsigned char tmp3=DATA(i,jm1);
    unsigned char tmp5=DATA(ip1,jm1);

    unsigned char tmp1=DATA(im1,j);
    unsigned char tmp2=DATA(im1,jp1);

    unsigned char tmp4=DATA(i,jp1);

    unsigned char tmp6=DATA(ip1,j);
    unsigned char tmp7=DATA(ip1,jp1);

    register unsigned char n_alive_cells = tmp0+tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7;

    DATA(i,j) = ALIVE*(n_alive_cells==3) + DATA(i,j)*(n_alive_cells==2);
}

int main(int argc, char **argv)
{

    // initialize MPI
    MPI_Init(&argc, &argv);

    // each rank gets info about itself
    int rank;
    int size;

    unsigned char * data, * data_prev;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // get arguments
    // MPI standard does not specify how arguments are processed. 
    // Ideally it would be better to have master process the args and then 
    // broadcast the arguments. However, for now this is ok.
    get_args(argc, argv);

    // if(n>99999)
    // {
    //     printf("n cannot be greater than 99999. Using this value");
    //     n = 99999;
    // }
    //
    // if(s>99999)
    // {
    //     printf("n cannot be greater than 99999. Using this value");
    //     s = 99999;
    // }

    // setting up neighbors for 1D splitting
    // Note that with 1D splitting, we have a limitation on how many processes 
    // this program will work with. 
    // If the processes are more than the number of rows, then the program will 
    // not use some of these processes.
    const int prev = (size - 1)*(rank==0) + (rank - 1)*!(rank==0);
    const int next = 0*(rank==(size-1)) + (rank + 1)*!(rank==(size-1));

    if(action == INIT)
    {

        // for each rank, we want to determine how many rows it will need to 
        // handle.
        // Since we are running the program with -i option (INIT), 
        // we also specify the number of rows and columns. 
        // So this operation is well defined.
        
        const int my_rows = get_my_rows(rows, rank, size);
        const int my_row_offset = get_my_row_offset(rows, rank, size);
        const int augmented_rows = my_rows + 2;

        // printf("I am rank %d and I am getting %d rows at offset %d\n", rank, my_rows, my_row_offset);

        const MPI_Offset my_file_offset = my_row_offset * cols * sizeof(unsigned char);

        data = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));

        if(!data)
        {
            printf("Allocation failed.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }

        // initialize the halo regions to being DEAD
        for(int j = 0; j<cols; ++j)
        {
            DATA(0,j) = DATA(my_rows + 1,j) = DEAD;
        }

        // HEADER INFO CALCULATION

        // intuition: https://stackoverflow.com/questions/3919995/determining-sprintf-buffer-size-whats-the-standard
        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);

        if(!header)
        {
            printf("Allocation failed.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        const MPI_Offset header_offset = header_size * sizeof(char);
        const MPI_Offset my_total_file_offset = my_file_offset + header_offset;


        // INITIALIZE MATRIX RANDOMLY
        
        // for reproducible results. If we want truly random then we also 
        // need to include some changing information such as time
        srand48(10*rank); 

        for(int i = 1; i<my_rows+1; ++i)
        {
            for(int j = 0; j<cols; ++j)
            {
                DATA(i,j) = drand48() > 0.5 ? ALIVE : DEAD ;
            }
        }

        save_grid(fname, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data, my_rows, cols);
        free(header);
        free(data);

    }
    else if (e == STATIC && action == RUN)
    {

        // read file 
        // in general, we cannot assume we will use the same number 
        // of processors as the initialization phase
        // so we need to first read the header, broadcast to the other 
        // processes, distribute an area of the matrix for each process 
        // and read from the file
        
        int opt_args[2] = {0,0};

        if(rank == 0)
        {

            FILE * fh_posix = fopen(fname, "r");

            // we know that the magic number is P5 so we set the offset as the  
            // length of P5 * sizeof(char)
            if(!fh_posix)
            {
                fprintf(stderr, "Error opening %s\n", fname);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);
            }

            int args_scanned = fscanf(fh_posix, "P5 %d %d 1\n",opt_args, opt_args+1 );
            if(args_scanned != 2)
            {
                printf("fscanf failed.");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);
            }

            fclose(fh_posix);
            // printf("I am rank 0 and I have received %d rows %d cols\n", *opt_args, *(opt_args+1));
        }

        MPI_Bcast(opt_args, 2, MPI_INT, 0, MPI_COMM_WORLD);

        // printf("I am rank %d and I have received %d rows %d cols\n", rank, *opt_args, *(opt_args+1));

        rows = opt_args[0];
        cols = opt_args[1];

        if(size > rows)
        {
            printf("This program cannot handle more processes than rows. Make sure P <= R");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);
        }

        const int my_rows = get_my_rows(rows, rank, size);
        const int my_row_offset = get_my_row_offset(rows, rank, size);

        // printf("I am rank %d and I am getting %d rows at offset %d\n", rank, my_rows, my_row_offset);

        // since we need to deal with halos we need to allocate two more 
        // rows than what I need.
        
        const int augmented_rows = my_rows + 2;

        const MPI_Offset my_file_offset = my_row_offset * cols * sizeof(char);

        data = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));
        data_prev = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));

        if(!data || !data_prev)
        {
            printf("Allocation failed.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }
        
        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);

        if(!header)
        {
            printf("Allocation failed.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        const MPI_Offset header_offset = header_size * sizeof(char);
        const MPI_Offset my_total_file_offset = my_file_offset + header_offset;

        MPI_File fh;

        const int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

        if(err != MPI_SUCCESS)
        {
            fprintf(stderr, "Error opening %s\n", fname);
            MPI_Abort(MPI_COMM_WORLD, err);
        }

        MPI_File_read_at_all(fh, my_total_file_offset, data_prev + cols, my_rows*cols, MPI_CHAR, MPI_STATUS_IGNORE);

        MPI_File_close(&fh);

        // At this point, all the processes have read their portion of the 
        // matrix. The halo regions have been set to zero, and we are ready 
        // to start the game of life.
        
        // file name (it will always be the same length so we only need it once)
        char * snapshot_name = malloc(32);
        if(!snapshot_name)
        {
            printf("Not enough space.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }

        MPI_Request prev_send_request, next_send_request;
        MPI_Request prev_recv_request, next_recv_request;

        const int prev_tag = 0; 
        const int next_tag = 1;

        unsigned char *tmp_data = NULL;

        const int MAX_THREADS = omp_get_max_threads();

        int chunk = (my_rows-2)*cols/MAX_THREADS;
        int remainder = chunk%CACHE_LINE_SIZE;
        // ensure chunk is a multiple of CACHE_LINE_SIZE bytes 
        chunk = remainder==0 ? chunk : chunk+CACHE_LINE_SIZE-remainder;

        int small_chunk = cols/MAX_THREADS;
        remainder = small_chunk%CACHE_LINE_SIZE;
        // ensure chunk is a multiple of CACHE_LINE_SIZE bytes 
        small_chunk = remainder==0 ? small_chunk : small_chunk+CACHE_LINE_SIZE-remainder;

        int save_counter = 0;
        const unsigned int my_rows_x_cols = my_rows*cols;
        const unsigned int my_rows_x_cols_p_cols = my_rows_x_cols + cols;

        double start_time = MPI_Wtime();
        for(int t = 1; t < n+1; ++t)
        {

            MPI_Isend(data_prev + cols, cols, MPI_CHAR, prev, prev_tag, MPI_COMM_WORLD, &prev_send_request);
            MPI_Isend(data_prev + my_rows_x_cols, cols, MPI_CHAR, next, next_tag, MPI_COMM_WORLD, &next_send_request);

            MPI_Irecv(data_prev, cols, MPI_CHAR, prev, next_tag, MPI_COMM_WORLD, &prev_recv_request);
            MPI_Irecv(data_prev + my_rows_x_cols_p_cols, cols, MPI_CHAR, next, prev_tag, MPI_COMM_WORLD, &next_recv_request);

            // we use the for loop over cells instead of the double for loop 
            // over columns and rows
            #pragma omp parallel for schedule(dynamic, chunk)
            for(int cell = 2*cols; cell< my_rows_x_cols; ++cell)
            {
                    int i = cell/cols;
                    int tmp2 = -i*cols;
                    int j = cell + tmp2;
                    upgrade_cell_static(data_prev, data, i, j);
            }

            MPI_Wait(&prev_recv_request, MPI_STATUS_IGNORE);
            MPI_Wait(&next_recv_request, MPI_STATUS_IGNORE);
            MPI_Request_free(&prev_send_request);
            MPI_Request_free(&next_send_request);

            #pragma omp parallel for schedule(dynamic, small_chunk)
            for(int col=0; col<cols; ++col)
            {
                upgrade_cell_static(data_prev, data, 1, col);
            }

            #pragma omp parallel for schedule(dynamic, small_chunk)
            for(int col=0; col<cols*(my_rows>1); ++col)
            {
                upgrade_cell_static(data_prev, data, my_rows, col);
            }

            ++save_counter;
            if(s>0 && save_counter == s && t<100000)
            {
                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data, my_rows, cols);
                save_counter = 0;
            }

            tmp_data = data;
            data = data_prev;
            data_prev = tmp_data;
            tmp_data = NULL;
        }
        double end_time = MPI_Wtime();
        double local_time = end_time-start_time;
        double global_time_avg;
        MPI_Reduce(&local_time, &global_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(rank==0)
        {
            global_time_avg/=size;
            printf("%d,%d,%d,%lf\n",size,MAX_THREADS,size*MAX_THREADS,global_time_avg);
        }
        if(s==0)
        {
            sprintf(snapshot_name, "snapshot_%05d", n);
            save_grid(snapshot_name, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data_prev, my_rows, cols);
        }
    
        free(snapshot_name);
        free(header);
        free(data);
        free(data_prev);
    }
    else if(e == ORDERED && action == RUN)
    {

        int opt_args[2] = {0,0};

        if(rank == 0)
        {

            FILE * fh_posix = fopen(fname, "r");

            if(!fh_posix)
            {
                fprintf(stderr, "Error opening %s\n", fname);
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);
            }

            int args_scanned = fscanf(fh_posix, "P5 %d %d 1\n",opt_args, opt_args+1 );
            if(args_scanned != 2)
            {
                printf("fscanf failed.");
                MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);
            }

            fclose(fh_posix);
        }

        MPI_Bcast(opt_args, 2, MPI_INT, 0, MPI_COMM_WORLD);

        rows = opt_args[0];
        cols = opt_args[1];

        if(size > rows)
        {
            printf("This program cannot handle more processes than rows. Make sure P <= R");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);
        }

        const int my_rows = get_my_rows(rows, rank, size);
        const int my_row_offset = get_my_row_offset(rows, rank, size);

        const int augmented_rows = my_rows + 2;

        const MPI_Offset my_file_offset = my_row_offset * cols * sizeof(char);

        data = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));

        if(!data)
        {
            printf("Allocation failed.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }
        
        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);

        if(!header)
        {
            printf("Allocation failed.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        const MPI_Offset header_offset = header_size * sizeof(char);
        const MPI_Offset my_total_file_offset = my_file_offset + header_offset;

        MPI_File fh;

        const int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

        if(err != MPI_SUCCESS)
        {
            fprintf(stderr, "Error opening %s\n", fname);
            MPI_Abort(MPI_COMM_WORLD, err);
        }

        MPI_File_read_at_all(fh, my_total_file_offset, data + cols, my_rows*cols, MPI_CHAR, MPI_STATUS_IGNORE);

        MPI_File_close(&fh);

        char * snapshot_name = malloc(32);
        if(!snapshot_name)
        {
            printf("Not enough space.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }

        MPI_Request prev_send_request, next_send_request;
        MPI_Request prev_recv_request, next_recv_request;

        const int prev_tag = 0; 
        const int next_tag = 1;

        int save_counter = 0;
        const unsigned int my_rows_x_cols = my_rows*cols;
        const unsigned int my_rows_x_cols_p_cols = my_rows_x_cols + cols;

        double start_time = MPI_Wtime();
        for(int t = 1; t < n+1; ++t)
        {
            // for ordered evolution, we cannot parallelize, and each process 
            // must work on their part of the grid in serial order. 
            
            // first, each process posts a non blocking receive operation for 
            // the top halo row and the bottom halo row. 
            // Since we have ordered evolution each process will then wait 
            // until the top halo row has been received to start working on 
            // its grid.
            MPI_Irecv(data, cols, MPI_CHAR, prev, next_tag, MPI_COMM_WORLD, &prev_recv_request);
            MPI_Irecv(data + my_rows_x_cols_p_cols, cols, MPI_CHAR, next, prev_tag, MPI_COMM_WORLD, &next_recv_request);

            // however, each process will also have to send its top row 
            // (the bottom halo row for the previous process) which is necessary 
            // for the previous process to calculate the grid.
            // However, this excludes rank 0, because it needs to send its first
            // row to rank size-1 row only AFTER it has been updated 
            if(rank != 0){ MPI_Isend(data + cols, cols, MPI_CHAR, prev, prev_tag, MPI_COMM_WORLD, &prev_send_request);}

            // to get things started for rank 0, this is the first forward 
            // message of the top halo row.
            if(rank == size-1) { MPI_Isend(data + my_rows_x_cols, cols, MPI_CHAR, next, next_tag, MPI_COMM_WORLD, &next_send_request); }

            // All processes will wait until receive of top halo row is complete 
            // (since they need it to start ordered evolution)
            //
            // the first time, only process 0 should pass.
            MPI_Wait(&prev_recv_request, MPI_STATUS_IGNORE);
            // printf("I am rank %d and I am in\n", rank);

            // Since we have received the top halo row, we can process all 
            // data EXCEPT the last row, since we are not sure that we have 
            // received the bottom halo yet.
            // If my_rows = 1 this will be skipped (as it should)
            for(int row = 1; row < my_rows; ++row)
            {
                for(int col = 0; col < cols; ++col)
                {
                    upgrade_cell_ordered(data, row, col);
                }
            }
            if(my_rows>1)
            {
                // in this case, the for loop above has been executed and 
                // now that the first row has been processed, rank zero can send it 
                // to rank size-1
                if(rank == 0){ MPI_Isend(data + cols, cols, MPI_CHAR, prev, prev_tag, MPI_COMM_WORLD, &prev_send_request);}

                // We have finished processing all our grid except for last row. 
                // For this, we need to make sure that bottom halo has been received.

                // wait for receive of bottom halo
                MPI_Wait(&next_recv_request, MPI_STATUS_IGNORE);

                // process last row
                for(int col=0; col<cols; ++col)
                {
                    upgrade_cell_ordered(data, my_rows, col);
                }

                // now that we have computed the last row, we need to send it 
                // to the next process. This will be the top halo for the next 
                // process, which in turn will finally pass the wait statement 
                // above.
                
                // all except process n-1 because this will be taken care of in the 
                // next generation
                if(rank<(size-1)){ MPI_Isend(data + my_rows_x_cols, cols, MPI_CHAR, next, next_tag, MPI_COMM_WORLD, &next_send_request);}

            }
            else
            {
                // in this case we have skipped the for loop
                // we need to wait to receive the bottom halo from next
                MPI_Wait(&next_recv_request, MPI_STATUS_IGNORE);
                
                // we process row 1
                for(int col=0; col<cols; ++col)
                {
                    upgrade_cell_ordered(data, my_rows, col);
                }

                // rank zero sends this row to rank n-1
                if(rank == 0){ MPI_Isend(data + cols, cols, MPI_CHAR, prev, prev_tag, MPI_COMM_WORLD, &prev_send_request);}

                // all ranks except last, send their last row (the same one )
                // to the next rank
                if(rank<(size-1)){ MPI_Isend(data + my_rows_x_cols, cols, MPI_CHAR, next, next_tag, MPI_COMM_WORLD, &next_send_request);}
            }

            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Request_free(&prev_send_request);
            MPI_Request_free(&next_send_request);

            ++save_counter;
            if(s>0 && save_counter == s && t<100000)
            {
                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data, my_rows, cols);
                save_counter = 0;
            }

        }
        double end_time = MPI_Wtime();
        double local_time = end_time-start_time;
        double global_time_avg;
        MPI_Reduce(&local_time, &global_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(rank==0)
        {
            global_time_avg/=size;
            printf("%d,%d,%d,%lf\n",size,1,size,global_time_avg);
        }
        if(s==0)
        {
            sprintf(snapshot_name, "snapshot_%05d", n);
            save_grid(snapshot_name, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data, my_rows, cols);
        }
    
        free(snapshot_name);
        free(header);
        free(data);
    }
    else
    {
        printf("Unknown action. Abort");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    if ( fname != NULL )
      free ( fname );

    MPI_Finalize();
    return 0;
}

