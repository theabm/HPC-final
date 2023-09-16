#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <mpi.h>
#include <unistd.h>
#include <immintrin.h>

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
    unsigned char count_of_neighbors;

    register int jm1 = j==0 ? cols-1 : j-1;
    register int jp1 = j==(cols-1) ? 0 : j+1;
    register int im1 = i-1;
    register int ip1 = i+1;
    
    // if we get here we have non zero element
    // so it is either 
    // 1. a cell with some neighbors 
    // 2. a live cell 
    // 3. both
    //
    // However, since we already copied the contents 
    // of data prev, we only need to keep track of 
    // events that *change* the status of the cell 
    // i.e births and deaths
    
    count_of_neighbors = DATA_PREV(i,j) >> 1;

    // we first deal with deaths, however, we need 
    // to make sure that we only update in the case 
    // that the cell was alive and died 

    // cell was alive, so last bit is 1
    if(DATA_PREV(i,j) & 0x01)
    {
        // if greater than 3 or less than 2, we kill it. 
        // otherwise we dont do anything.
        if( count_of_neighbors>3 || count_of_neighbors<2 )
        {

            // to set cell to dead we need to
            // 1. flip its last bit. We can do this by doing 
            // a bit wise and with ~0x01. i.e 11111110
            // since 1 & 0 = 0 and 1 & 1 = 1, the first 7 bits 
            // remain unchanged. while the last one is turned to 
            // zero.
            
            DATA(i,j) &= ~0x01;

            // 2. get neighbors and decrease the neighbor count 
            // of each one. 
            // This can be achieved by subtracting 00000010 (or 0x02)
            // from each cell.
            // lets check possible limiting cases: 
            // 1. will it ever be negative? We need a value smaller 
            // than 0x02 for this to happen, namely 0x01.
            // However, by construction, this neighbors had at least 
            // one live "on" cell (the one we are killing right now) 
            // So the smallest possible value for any neighboring cell 
            // is 0x02 (i.e 00000010).
            // So we will never obtain a negative value causing 
            // strange behavior. Thats good.
            //
            // 2. will the subtraction ever change the state of the 
            // cell? No, because we are guaranteed in this way that 
            // the last bit will remain the same
            
            DATA(im1, jm1)-=0x02;
            DATA(im1,j)-=0x02;
            DATA(im1,jp1)-=0x02;
            DATA(i,jm1)-=0x02;
            DATA(i,jp1)-=0x02;
            DATA(ip1,jm1)-=0x02;
            DATA(ip1,j)-=0x02;
            DATA(ip1,jp1)-=0x02;
        }   
    }
    else
    { 
        // the cell is dead.
        if(count_of_neighbors == 3) 
        {

            // set to living and increase counter of neighbors

            // we set the last bit to 1 with a bit wise or
            DATA(i,j) |= 0x01; 

            // 2. get neighbors and increase the neighbor count 
            // of each one. 
            // To do this, we need to sum 0x02

            DATA(im1, jm1)+=0x02;
            DATA(im1,j)+=0x02;
            DATA(im1,jp1)+=0x02;
            DATA(i,jm1)+=0x02;
            DATA(i,jp1)+=0x02;
            DATA(ip1,jm1)+=0x02;
            DATA(ip1,j)+=0x02;
            DATA(ip1,jp1)+=0x02;
        }
    }
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

void vectorized_bitwise_and(unsigned char * data_prev, int rows, int cols){
    // we create an alias of data_prev with vectorized data type 
    // __m256 is a 256 long bit register, aka it holds 32 chars (cells)
    unsigned char * start = data_prev+cols;

    __m256i * d_p = (__m256i*)start;

    __m256i result;

    __m256i bits = _mm256_set1_epi8(0x01);

    const int n = (rows*cols)/32;
    const int remainder = (rows*cols)%32;

    for(int k = 0; k<n; ++k){
        result = _mm256_and_si256(d_p[k], bits);
        _mm256_storeu_si256(d_p+k, result);
    }

    // we need to find the last address that was processed 
    // so at the end, I have processed 32*n elements of the array.
    const int ofs = n*32;

    for(int k=0; k<remainder; ++k){
        *(data_prev+ofs+k) &= 0x01;
    }
}

void bitwise_and(unsigned char * data_prev, int start, int end)
{
    for(int k=start; k<end; ++k)
    {
        *(data_prev+k) &= 0x01;
    }
}

void preprocess_cell(unsigned char * data_array, int i, int j, int cols)
{
    // calculate j+1, j-1 with wrapping
    register int jm1 = j==0 ? cols-1 : j-1;
    register int jp1 = j==(cols-1) ? 0 : j+1;
    register int im1 = i-1;
    register int ip1 = i+1;

    // get neighbors state
    unsigned char tmp0=data_array[im1*cols+jm1] & 0x01;
    unsigned char tmp3=data_array[i*cols+jm1] & 0x01;
    unsigned char tmp5=data_array[ip1*cols+jm1] & 0x01;

    unsigned char tmp1=data_array[im1*cols+j] & 0x01;
    unsigned char tmp2=data_array[im1*cols+jp1] & 0x01;

    unsigned char tmp4=data_array[i*cols+jp1] & 0x01;

    unsigned char tmp6=data_array[ip1*cols+j] & 0x01;
    unsigned char tmp7=data_array[ip1*cols+jp1] & 0x01;

    // this will be a number between 8 and 0
    unsigned char n_alive_cells = tmp0+tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7;

    // so now I need to shift everything left by one
    n_alive_cells = n_alive_cells<<1;

    // lastly, if the cell is dead, we leave the last bit as is 
    // and if the cell is alive, we need to set the last bit to 1
    // keeping everything else the same. 
    // so if cells is dead we have                      00000000
    // and we have if we do bitwise or with nalivecells 000xxxx0 
    // (x denotes any possible values)
    // we obtain the again n_alive_cells
    // if the cells is alive 00000001
    // and we OR with        000xxxx0
    // we obtain             000xxxx1
    // to do this, we need to do bitwise OR with DATA
    data_array[i*cols+j] |= n_alive_cells;
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

        // this is blocking (maybe use non blocking later)
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

        if(my_rows<5)
        {
            printf("this program logic cannot handle the case where each process gets less than 5 rows. Either use less processes or use algo 1");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_REQUEST);
        }

        const int my_end = (my_rows+1)*cols;
        const int augmented_rows = my_rows + 2;

        const MPI_Offset my_file_offset = my_row_offset * cols * sizeof(char);

        data = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));
        data_prev = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));

        if(!data || !data_prev)
        {
            printf("Allocation failed.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }
        
        // initialize the halo regions to being DEAD
        for(int j = 0; j<cols; ++j)
        {
            DATA(0,j) = DATA(my_rows + 1,j) = DATA_PREV(0,j)
                = DATA_PREV(my_rows + 1,j) = DEAD;
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

        // each process has its part of the grid, but must first 
        // preprocess it to get it into the right format

        MPI_Request prev_send_request, next_send_request;
        MPI_Request prev_recv_request, next_recv_request;
        MPI_Request prev_send_request1, next_recv_request1;

        const int prev_tag = 0; 
        const int next_tag = 1;

        // we initiate exchange of halo cells
        // send row 1 and row n
        MPI_Isend(data_prev + cols, cols, MPI_CHAR, prev, prev_tag, MPI_COMM_WORLD, &prev_send_request);
        MPI_Isend(data_prev + my_rows*cols, cols, MPI_CHAR, next, next_tag, MPI_COMM_WORLD, &next_send_request);

        // receive in row 0 and row n+1
        MPI_Irecv(data_prev, cols, MPI_CHAR, prev, next_tag, MPI_COMM_WORLD, &prev_recv_request);
        MPI_Irecv(data_prev + cols*my_rows + cols, cols, MPI_CHAR, next, prev_tag, MPI_COMM_WORLD, &next_recv_request);

        // since we havent received yet, we process row 2 to my_rows-1 in 
        // meantime 
        for(int i = 2; i < my_rows; ++i)
        {
            for(int j = 0; j < cols; ++j)
            {
                preprocess_cell(data_prev, i, j, cols);
            }
        }

        // we wait for receive to complete
        MPI_Wait(&prev_recv_request, MPI_STATUS_IGNORE);
        MPI_Wait(&next_recv_request, MPI_STATUS_IGNORE);
        MPI_Request_free(&prev_send_request);
        MPI_Request_free(&next_send_request);

        // process row 0 and my_rows
        for(int j = 0; j < cols; ++j)
        {
            preprocess_cell(data_prev, 1, j, cols);
            preprocess_cell(data_prev, my_rows, j, cols);
        }

        // grid been preprocessed
        
        // file name (it will always be the same length so we only need it once)
        char * snapshot_name = malloc(32);
        if(!snapshot_name)
        {
            printf("Not enough space.");
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_NO_SPACE);
        }

        unsigned char *tmp_data = NULL;

        int save_counter = 0;
        const unsigned int my_rows_x_cols = my_rows*cols;
        const unsigned int my_rows_x_cols_p_cols = my_rows_x_cols + cols;
        const int my_grid_size_bytes = my_rows*cols*sizeof(unsigned char);

        double start_time = MPI_Wtime();
        for(int t = 1; t < n+1; ++t)
        {

            // the idea is the same as omp.v2, however, whenever we 
            // say "copy row 0 to row n" etc, we mean make a communication
            //
            // the steps can be summarized as follows:
            // 0. copy data_prev into data
            // 1. send row my_rows to rank next (non blocking) and put it in row 0
            // 2. process rows 3 to my_rows - 2 (to hide latency)
            // 3. wait until receive is complete
            // 4. process rows 1 and 2
            // 5. send row 0 to prev and put it in row n (non blocking)
            // 6. send rows 1 to prev (non blocking)
            // 7. wait for receive of row n
            // 8. process row n-1
            // 9. wait for receive of row n+1
            // 10. process row n
            // 11. send back row n+1 to next and put it in row 1 
            
            memcpy(data+cols, data_prev+cols, my_grid_size_bytes);

            // send row: my_rows to next rank
            MPI_Isend(data + my_rows_x_cols, cols, MPI_CHAR, next, next_tag, MPI_COMM_WORLD, &next_send_request);

            // receive from prev and put in row 0 (halo region)
            MPI_Irecv(data, cols, MPI_CHAR, prev, next_tag, MPI_COMM_WORLD, &prev_recv_request);

            // process rows 3 until my_rows - 2
            for(int i = 3; i < my_rows-1; ++i)
            {
                for(int j = 0; j < cols; ++j)
                {
                    if(DATA_PREV(i,j)==0)
                    {
                        continue;
                    }

                    upgrade_cell_static(data_prev, data, i, j);
                }

            }

            // wait to receive row 0
            MPI_Wait(&prev_recv_request, MPI_STATUS_IGNORE);
            MPI_Request_free(&next_send_request);

            // process rows 1 and 2
            for(int j = 0; j < cols; ++j)
            {
                if(DATA_PREV(1,j)==0 && DATA_PREV(2,j)==0)
                {
                    continue;
                }
                else if(DATA_PREV(1,j) !=0 && DATA_PREV(2,j)!=0)
                {
                    upgrade_cell_static(data_prev, data, 1, j);
                    upgrade_cell_static(data_prev, data, 2, j);
                }
                else if(DATA_PREV(1,j)!=0)
                {
                    upgrade_cell_static(data_prev, data, 1, j);
                }
                else
                {
                    upgrade_cell_static(data_prev, data, 2, j);
                }
            }
                
            // send row 0 to prev into row n
            MPI_Isend(data, cols, MPI_CHAR, prev, prev_tag, MPI_COMM_WORLD, &prev_send_request);
            MPI_Irecv(data + my_rows_x_cols, cols, MPI_CHAR, next, prev_tag, MPI_COMM_WORLD, &next_recv_request);

            // send row 1 to prev into row n+1
            MPI_Isend(data + cols, cols, MPI_CHAR, prev, prev_tag+1, MPI_COMM_WORLD, &prev_send_request1);
            MPI_Irecv(data + my_rows_x_cols_p_cols, cols, MPI_CHAR, next, prev_tag+1, MPI_COMM_WORLD, &next_recv_request1);

            // wait for receive of row n
            MPI_Wait(&next_recv_request, MPI_STATUS_IGNORE);
            MPI_Request_free(&prev_send_request);

            // process row n-1
            for(int j = 0; j < cols; ++j)
            {
                if(DATA_PREV(my_rows-1,j)==0)
                {
                    continue;
                }

                upgrade_cell_static(data_prev, data, my_rows-1, j);
            }

            // wait until receive of row n+1
            MPI_Wait(&next_recv_request1, MPI_STATUS_IGNORE);
            MPI_Request_free(&prev_send_request1);

            // process row n
            for(int j = 0; j < cols; ++j)
            {
                if(DATA_PREV(my_rows,j)==0)
                {
                    continue;
                }

                upgrade_cell_static(data_prev, data, my_rows, j);
            }
            
            // send row n+1 to next into row 1
            MPI_Isend(data + my_rows_x_cols + cols, cols, MPI_CHAR, next, next_tag, MPI_COMM_WORLD, &next_send_request);
            // blocking receive
            MPI_Recv(data + cols, cols, MPI_CHAR, prev, next_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Request_free(&next_send_request);

            // Check if need to save, and if we do, save grid to pgm
            
            ++save_counter;

            if(s>0 && save_counter == s && t<100000)
            {
                memcpy(data_prev + cols, data + cols, my_grid_size_bytes);

                bitwise_and(data_prev, cols, my_end);

                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data_prev, my_rows, cols);
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
            printf("%d,%lf\n",size,global_time_avg);
        }
        if(s==0)
        {
            memcpy(data + cols, data_prev + cols, my_grid_size_bytes);

            bitwise_and(data, cols, my_end);

            sprintf(snapshot_name, "snapshot_%05d", n);
            save_grid(snapshot_name, MPI_COMM_WORLD, rank, header, header_size, my_total_file_offset, data, my_rows, cols);
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
            printf("%d,%lf\n",size,global_time_avg);
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

