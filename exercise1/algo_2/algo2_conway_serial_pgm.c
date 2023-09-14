// this version improves upon the other one in an algorithmic sense. 
// before, we were using one byte to store data for one single cell. 
// And for each cell, we would do 8 updates. Naturally, we spend the most 
// time calling the function upgrade_cell. 
//
// However, it is worth noting, that the likelyhood of a cell surviving, or 
// becoming alive, is much lower than the likelihood of it being dead. 
// Suppose we look at a living cell:
// We have 2^8 combinations of alive and dead for its neighbors. 
// Out of those, only 8 choose 3 + 8 choose 2 combinations keep the cell alive. 
// 8!/(3! * 5!) + 8!/(2! * 6!) = 84 / 2^8 = 32.8% chance of survival.
// 
// For a dead cell to come to life instead, we need 3 live neighbors, which is 
// 56/2^8 = 21.8% 
//
// With an initial uniform distribution, a cell is alive or dead with 50% 
// probability. 
// So the probability of becoming/staying alive is .5*.32 + .5*.21 = 26%. 
//
// This means that on average, most cells will be dead and have dead neighbors, 
// so their status wont change! which means, that if we encode information better, 
// we can avoid making repeated calls to upgrade_cell for cells that dont change.
//
// This idea (taken from : source) is as follows:
// use a byte to encode a cell state, however, 
// the last bit will encode if the cells is alive of dead 
// the 4 bits before that will encode the sum of the neighbors (it can only 
// go up to 8 so we only need 4 bits), and the first 3 bits will be empty.
//
// In this way, by looking at the three bytes, we can already tell the situation 
// of the neighbors of the cell, and we need only to scan for non zero bytes. 
// of course, once we change the status of a cell, we have to change all of its 
// 8 neighbors (the counts). However, as said before, we will likely not do this 
// for all cells, so its still a win.
//
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <immintrin.h>
#include <time.h>


#define DATA(i,j) (data[(i)*cols + (j)])
#define DATA_PREV(i,j) (data_prev[(i)*cols + (j)])

// for pgm files, 0 is black and MAXVAL is white.
// So dead is 0 and alive is 1
#define DEAD 0
#define ALIVE 1
#define MAX_VAL 1

#define CACHE_LINE_SIZE 64

#define K_DFLT 100

#define INIT 1
#define RUN  2

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

void save_grid(char * restrict fname, char * restrict header, int header_size, unsigned char * restrict data, int rows, int cols)
{

    FILE * fh = fopen(fname, "wb");

    if(fh == NULL)
    {
        fprintf(stderr, "Error opening %s\n", fname);
        exit(0);
    }

    fwrite(header,1,header_size,fh);
    
    fwrite(data+cols, sizeof(unsigned char), rows*cols,fh);

    fclose(fh);

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
    unsigned char count_of_neighbors;

    register int jm1 = j==0 ? cols-1 : j-1;
    register int jp1 = j==(cols-1) ? 0 : j+1;
    register int im1 = i-1;
    register int ip1 = i+1;

    count_of_neighbors = DATA(i,j) >> 1;

    if(DATA(i,j) & 0x01)
    {
        if( count_of_neighbors>3 || count_of_neighbors<2 )
        {
            
            DATA(i,j) &= ~0x01;
            
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
        if(count_of_neighbors == 3) 
        {

            DATA(i,j) |= 0x01; 

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

void display_args(){
    printf("action (i : 1\tr : 2) -- %d\nk (size) -- %d\ne (0 : ORDERED\t1 : STATIC) -- %d\nf (filename) -- %s\nn (steps) -- %d\ns (save frequency) -- %d\n", action, k, e, fname, n, s );
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

struct timespec diff(struct timespec start, struct timespec end)
{
        struct timespec temp;
        if ((end.tv_nsec-start.tv_nsec)<0) {
                temp.tv_sec = end.tv_sec-start.tv_sec-1;
                temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
        } else {
                temp.tv_sec = end.tv_sec-start.tv_sec;
                temp.tv_nsec = end.tv_nsec-start.tv_nsec;
        }
        return temp;
}

int main(int argc, char **argv)
{

    unsigned char * data, * data_prev;
    struct timespec start_time, end_time;
    double elapsed;

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

    if(action == INIT)
    {
        // init is done in the same way since we need to save a pgm file

        // we add two rows for halo regions
        const int augmented_rows = rows + 2;

        data = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));

        if(!data)
        {
            printf("Allocation failed.");
            exit(0);
        }

        // initialize the halo regions to being DEAD
        for(int j = 0; j<cols; ++j)
        {
            DATA(0,j) = DATA(rows + 1,j) = DEAD;
        }

        // HEADER INFO CALCULATION

        // intuition: https://stackoverflow.com/questions/3919995/determining-sprintf-buffer-size-whats-the-standard
        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = (char *) malloc(header_size + 1);

        if(!header)
        {

            printf("Allocation failed.");
            exit(0);

        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);

        srand48(10); 

        for(int i = 1; i<rows+1; ++i)
        {
            for(int j = 0; j<cols; ++j)
            {
                DATA(i,j) = drand48() > 0.5 ? ALIVE : DEAD ;
            }
        }

        save_grid(fname, header, header_size, data, rows, cols);
        free(header);
        free(data);
    }
    else if (e == STATIC && action == RUN)
    {
        // reading is done in the same way.

        int opt_args[2] = {0,0};


        FILE * fh_posix = fopen(fname, "r");

        // we know that the magic number is P5 so we set the offset as the  
        // length of P5 * sizeof(char)
        if(fh_posix == NULL)
        {
            fprintf(stderr, "Error opening %s\n", fname);
            exit(0);
        }

        int args_scanned = fscanf(fh_posix, "P5 %d %d 1\n",opt_args, opt_args+1 );
        if(args_scanned != 2)
        {
            printf("fscanf failed.");
            exit(0);
        }

        fclose(fh_posix);

        rows = opt_args[0];
        cols = opt_args[1];

        if(rows<3 || cols<3)
        {
            // if cols = 3, we have problem that wrapping around updates 
            // neighbors twice
            printf("Matrix is too small. Use normal serial version.");
            exit(0);
        }

        const int end = (rows+1)*cols;
        const int augmented_rows = rows + 2;

        data = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));
        data_prev = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));

        if(!data || !data_prev)
        {
            printf("Allocation failed.");
            exit(0);
        }

        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);

        if(!header)
        {
            printf("Allocation failed.");
            exit(0);
        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);

        FILE * fh = fopen(fname, "rb");

        if(!fh)
        {
            fprintf(stderr, "Error opening %s\n", fname);
            exit(0);
        }

        fseek(fh, header_size, 1);

        size_t args_read = fread(data_prev+cols, sizeof(unsigned char), rows*cols, fh);

        if(args_read != (size_t)rows*cols)
        {
            printf("fread failed.");
            exit(0);
        }

        fclose(fh);


        // at this point, we have read the matrix, however, we need to do 
        // some preprocessing to get into the format we described. 
        // As is the case with most optimizations, there is a price to pay, 
        // i.e there is no free lunch.
        
        // step one is to swap the halo rows

        memcpy(data_prev, data_prev + rows*cols, cols*sizeof(char));
        memcpy(data_prev+rows*cols+cols, data_prev+cols, cols*sizeof(char));

        for(int i = 1; i < rows+1; ++i)
        {
            for(int j = 0; j < cols; ++j)
            {
                preprocess_cell(data_prev, i, j, cols);
            }
        }

        // so now, data_prev contains our cells in the desired format. Apart 
        // from the halo rows (which need to be swapped again in the new format)

        // now we have preprocessed the grid to our desired format.
        char * snapshot_name = malloc(32);
        if(!snapshot_name)
        {
            printf("Not enough space.");
            exit(0);
        }

        unsigned char *tmp_data = NULL;

        const int row_len_bytes = cols*sizeof(unsigned char);
        int save_counter = 0;
        const unsigned int rows_x_cols = rows*cols;
        const unsigned int rows_x_cols_p_cols = rows_x_cols + cols;
        const int grid_size_bytes = rows*cols*sizeof(unsigned char);

        clock_gettime(CLOCK_MONOTONIC, &start_time);
        for(int t = 1; t < n+1; ++t)
        {
            
            // in this case, we need to handle the update with care. 
            // We can't simply process each cell individually as before, because 
            // we have some interdependency between data, i.e modifying the state 
            // of row i can modify row i-i and row i+1.

            // Before, I simply swapped the halo cells, copied the array into 
            // data (for ease of use), and then updated data while looking at 
            // data_prev. 
            // However, I realized this doesnt work! why? 
            //
            // Once we swap the halo rows, row 0 will be a copy of row n, and row 
            // n+1 will be a copy of row 1. 
            // When we process row 1, we could update row 0, and when we update 
            // row n we could update row n+1.
            // These changes need to get propagated to the corresponding row 
            // somehow, otherwise we get wrong results. 
            //
            // The way I chose to do this is: 
            // we copy all of data_prev into data. This is to make the update 
            // simpler since otherwise we would need to keep track of which  
            // neighbors are in one matrix and which are in another. Also, by 
            // the rules of conway, in most cases we do nothing only births and 
            // deaths affect change the state of cells.
            
            // Then we work only in data. We copy row n into halo 0, 
            // we process row 1 (and perhaps it affects halo 0),
            // so we copy halo 0 *back* into row n so it is updated. 
            // Then we process row 2 (which maybe changes row 1). 
            // So now, we send the updated version of row 1 as a halo n+1.
            // This will maybe get updated when we process row n, so at the end 
            // of processing the rest of the array (row 3 to n)
            // we copy row n+1 back into row 1.
            //
            // Note that for the mpi case, this will correspond to 4 
            // communications in total, some of them blocking.


            // copy data_prev into data (only internal columns)
            memcpy(data+cols, data_prev+cols, grid_size_bytes);

            // copy row n into halo 0
            memcpy(data, data + rows_x_cols, row_len_bytes);

            // update row 1 and 2
            // remember, we *look* at cell(i,j) in data_prev 
            // (only for internal cells)
            // but we update cell(i,j) and its neighbors in data
            for(int j = 0; j < cols; ++j)
            {
                if(DATA_PREV(1,j)==0)
                {
                    continue;
                }

                upgrade_cell_static(data_prev, data, 1, j);

            }

            for(int j = 0; j < cols; ++j)
            {
                if(DATA_PREV(2,j)==0)
                {
                    continue;
                }

                upgrade_cell_static(data_prev, data, 2, j);

            }

            // copy halo 0 back into row n
            memcpy(data + rows_x_cols, data, row_len_bytes);

            // copy row 1 into halo n+1 
            memcpy(data+rows_x_cols_p_cols, data+cols, row_len_bytes);

            // we go through the array (from row 3 to rows)
            for(int i = 3; i < rows+1; ++i)
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
            
            
            // we copy halo n+1 into row 1
            memcpy(data+cols, data+rows_x_cols_p_cols, row_len_bytes);

            ++save_counter;

            // now comes the (ugly) of writing to file.
            // While this process is faster because of the reasons cited above, 
            // I/O becomes much more complicated. For each cell, we need to 
            // extract the information contained in the last bit, and the write 
            // it to file.

            // there are surely numerous ways to do this, however, we choose the 
            // simplest. Instead of allocating even more memory, we simply copy 
            // all of the contents of data into data_prev (since its an unused 
            // region at the moment), then, we will do bitwise AND with 0x01 
            // with each cell to extract the last bit. Finally, we write data_prev
            
            if(s>0 && save_counter == s && t<100000)
            {
                memcpy(data_prev + cols, data + cols, grid_size_bytes);

                bitwise_and(data_prev, cols, end);

                // at this point we have extracted the info of the cells and 
                // we are now ready to print as usual.

                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, header, header_size, data_prev, rows, cols);
                save_counter = 0;
            }

            // Finally, we need to swap data and data_prev.
            tmp_data = data;
            data = data_prev;
            data_prev = tmp_data;
            tmp_data = NULL;

            // data_prev will contain the updated cells (ready for a new 
            // generation) and will begin the next iteration in the for loop 
            // by swapping the halo cells
        }
        clock_gettime(CLOCK_MONOTONIC, &end_time);
        elapsed = (double)diff(start_time,end_time).tv_sec + (double)diff(start_time,end_time).tv_nsec / 1000000000.0;
        printf("time: %lf\n", elapsed);
    
        free(snapshot_name);
        free(header);
        free(data);
        free(data_prev);
    }
    else if(e == ORDERED && action == RUN)
    {
        // the ordered version is virtually the same as the above, 
        // except we directly operate on data directly. 

        int opt_args[2] = {0,0};


        FILE * fh_posix = fopen(fname, "r");

        // we know that the magic number is P5 so we set the offset as the  
        // length of P5 * sizeof(char)
        if(!fh_posix)
        {
            fprintf(stderr, "Error opening %s\n", fname);
            exit(0);
        }

        int args_scanned = fscanf(fh_posix, "P5 %d %d 1\n",opt_args, opt_args+1 );
        if(args_scanned != 2)
        {
            printf("fscanf failed.");
            exit(0);
        }

        fclose(fh_posix);

        rows = opt_args[0];
        cols = opt_args[1];

        if(rows<3 || cols<3)
        {
            // if cols = 3, we have problem that wrapping around updates 
            // neighbors twice
            printf("Matrix is too small. Use normal serial version.");
            exit(0);
        }

        const int end = (rows+1)*cols;
        const int augmented_rows = rows + 2;

        data = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));
        data_prev = (unsigned char *) aligned_alloc(CACHE_LINE_SIZE, augmented_rows * cols * sizeof(unsigned char));

        if(!data || !data_prev)
        {
            printf("Allocation failed.");
            exit(0);
        }

        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);

        if(!header)
        {
            printf("Allocation failed.");
            exit(0);
        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);

        FILE * fh = fopen(fname, "rb");

        if(!fh)
        {
            fprintf(stderr, "Error opening %s\n", fname);
            exit(0);
        }

        fseek(fh, header_size, 1);

        size_t args_read = fread(data+cols, sizeof(unsigned char), rows*cols, fh);

        if(args_read != (size_t)rows*cols)
        {
            printf("fread failed.");
            exit(0);
        }

        fclose(fh);


        // at this point, we have read the matrix, however, we need to do 
        // some preprocessing to get into the format we described. 
        // As is the case with most optimizations, there is a price to pay, 
        // i.e there is no free lunch.
        
        // step one is to swap the halo rows
        memcpy(data, data + rows*cols, cols*sizeof(char));
        memcpy(data+rows*cols+cols, data+cols, cols*sizeof(char));

        for(int i = 1; i < rows+1; ++i)
        {
            for(int j = 0; j < cols; ++j)
            {
                preprocess_cell(data, i, j, cols);
            }
        }
        // so now, data_prev contains our cells in the desired format. Apart 
        // from the halo rows (which need to be swapped again in the new format)

        // now we have preprocessed the grid to our desired format.
        char * snapshot_name = malloc(32);
        if(!snapshot_name)
        {
            printf("Not enough space.");
            exit(0);
        }

        const int row_len_bytes = cols*sizeof(unsigned char);
        int save_counter = 0;
        const unsigned int rows_x_cols = rows*cols;
        const unsigned int rows_x_cols_p_cols = rows_x_cols + cols;
        const int grid_size_bytes = rows*cols*sizeof(unsigned char);

        clock_gettime(CLOCK_MONOTONIC, &start_time);
        for(int t = 1; t < n+1; ++t)
        {
            
            // We copy the bottom row into the top halo cell. 
            memcpy(data, data + rows_x_cols, row_len_bytes);

            // we can't do the same thing that we did in v1 where we process 
            // the whole internal portion except the last row, and then copy 
            // the first row into the bottom halo, and process the last row. 
            // Why? Because when we process the first row, it might happen that 
            // it changes the results of the halo cells above it (i.e. the last 
            // row). If we proceed to update all the rows, when we process row 
            // n-1, that row may influence row n. So then we have to merge 
            // the results for the top halo and row n. And then process row n.
            // To avoid this complication, we first process row 1, we copy the 
            // top halo back into row n, and then we proceed as usual. (i.e 
            // row 2 to n-1, copy row 1 into row n+1, and then process row n)
            // lastly processing row n may have changed row n+1, so we need to 
            // copy it back into row 1
            
            // process row 1
            for(int j = 0; j < cols; ++j)
            {
                if(DATA(1,j)==0)
                {
                    continue;
                }
                upgrade_cell_ordered(data, 1, j);
            }

            // copy top halo (possibly modified) back into last row
            memcpy(data + rows_x_cols, data, row_len_bytes);

            // proceed with rows 2 until n-1
            for(int i = 2; i < rows; ++i)
            {
                for(int j = 0; j < cols; ++j)
                {
                    if(DATA(i,j)==0)
                    {
                        continue;
                    }
                    upgrade_cell_ordered(data, i, j);
                }
            }

            // then we copy the first row into the last halo cell
            memcpy(data+rows_x_cols_p_cols, data+cols, row_len_bytes);

            // we update the last row now that we have the updated information.
            for(int j = 0; j < cols; ++j)
            {
                if(DATA(rows,j)==0)
                {
                    continue;
                }
                upgrade_cell_ordered(data, rows, j);
            }

            // also, processing the last row may have affected the bottom halo 
            // so once we are done, we have to copy this back into the first row
            memcpy(data + cols, data+rows_x_cols_p_cols, row_len_bytes);

            ++save_counter;

            if(s>0 && save_counter == s && t<100000)
            {
                memcpy(data_prev + cols, data + cols, grid_size_bytes);

                bitwise_and(data_prev, cols, end);

                // at this point we have extracted the info of the cells and 
                // we are now ready to print as usual.

                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, header, header_size, data_prev, rows, cols);
                save_counter = 0;
            }

            // data_prev will contain the updated cells (ready for a new 
            // generation) and will begin the next iteration in the for loop 
            // by swapping the halo cells
        }
        clock_gettime(CLOCK_MONOTONIC, &end_time);
        elapsed = (double)diff(start_time,end_time).tv_sec + (double)diff(start_time,end_time).tv_nsec / 1000000000.0;
        printf("time: %lf\n", elapsed);
    
        free(snapshot_name);
        free(header);
        free(data);
        free(data_prev);
    }
    else
    {
        printf("Unknown action. Abort");
        exit(0);
    }

    if (fname != NULL )
      free(fname);

    return 0;
}

