// this is an omp version of serial v2
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <immintrin.h>
#include <omp.h>


#define DATA(i,j) (data[(i)*cols + (j)])
#define DATA_PREV(i,j) (data_prev[(i)*cols + (j)])

// for pgm files, 0 is black and MAXVAL is white.
// So dead is 0 and alive is 1
#define DEAD 0
#define ALIVE 1
#define MAX_VAL 1

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

    if(fh == NULL){
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

    register int jm1 = (j-1)%(int)cols + cols*((j-1)<0);
    register int jp1 = (j+1)%(int)cols + cols*((j+1)<0);
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

    register int jm1 = (j-1)%(int)cols + cols*((j-1)<0);
    register int jp1 = (j+1)%(int)cols + cols*((j+1)<0);
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
    register int jm1 = (j-1)%(int)cols + cols*((j-1)<0);
    register int jp1 = (j+1)%(int)cols + cols*((j+1)<0);
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

int main(int argc, char **argv){

    unsigned char * data, * data_prev;

    get_args(argc, argv);
    
    if(action == INIT)
    {
        // init is done in the same way since we need to save a pgm file

        // we add two rows for halo regions
        const int augmented_rows = rows + 2;

        data = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));

        if(
                !data
                )
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

        const int end = (rows+1)*cols;
        const int augmented_rows = rows + 2;

        data = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));
        data_prev = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));

        if(
                data == NULL
                || data_prev == NULL
                )
        {
            printf("Allocation failed.");
            exit(0);
        }

        // initialize the halo regions to being DEAD
        for(int j = 0; j<cols; ++j)
        {
            DATA(0,j) = DATA(rows + 1,j) = DATA_PREV(0,j)
                = DATA_PREV(rows + 1,j) = DEAD;
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

        if(fh == NULL)
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
        char * snapshot_name = malloc(snprintf(NULL, 0, "snapshot_%05d", 0)+1);
        if(snapshot_name == NULL)
        {
            printf("Not enough space.");
            exit(0);
        }

        unsigned char *tmp_data = NULL;

        for(int t = 1; t < n+1; ++t)
        {
            
            // As before, we need to swap the two halo cells
            memcpy(data_prev, data_prev + rows*cols, cols*sizeof(unsigned char));
            memcpy(data_prev+rows*cols+cols, data_prev+cols, cols*sizeof(unsigned char));

            // we copy all the contents of data_prev into data
            memcpy(data+cols, data_prev+cols, cols*rows*sizeof(unsigned char));

            // to parallelize here, we cannot simply put a #pragma omp parallel for 
            // as we did for v1. 
            // Lets analyze why. 
            //
            // Suppose we have a 15x15 grid. After all is said and done, we have 
            // two 17x17 grids, data and dataprev
            // 
            // Suppose we have three threads, so each one will deal with a 5x5 
            // block (ignoring halo rows). 
            // 
            // If thread 0 were to process row 4 while thread 1 starts processing 
            // row 6 (its first row), we would have a race condition since 
            // both threads can potentially modify row 5. 
            //
            // Can this scenario happen? After all, both threads should start with 
            // their first row at the same time. Unfortunately, OMP only schedules 
            // the threads for us, however, the actual management of threads 
            // is done by the OS. So, if, by misfortune, thread 1 was delayed because 
            // that core was extremely busy, and it starts to process its first row 
            // when thread 0 is going through row 4, we may get wrong results. 
            //
            // Of course, this is unlikely, but it is preferable not to leave 
            // the correctness of the algorithm to chance. 
            //
            // Furthermore, note that this is also true for rows 5 and 7, since 
            // they both can modify row 6. 
            //
            // So what can we do? One solution is:
            //
            // Since the issue is processing rows that concurrently write in 
            // the same place, we can simply put space between us. 
            //
            // In other words, we can process the grid 3 at a time, and 
            // paralellize that for loop.
            // In this way, at the boundary, where the problem used to be, 
            // the threads can process with any delay, since between these rows 
            // there are 3 other rows of data in reality. 
            //
            // This may be a little worse for cache friendliness, but it is 
            // the easiest solution.
            
            #pragma omp parallel for schedule(dynamic)
            for(int i = 1; i < rows+1; i+=3)
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

            #pragma omp parallel for schedule(dynamic)
            for(int i = 2; i < rows+1; i+=3)
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

            #pragma omp parallel for schedule(dynamic)
            for(int i = 3; i < rows+1; i+=3)
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

            if(t%s == 0)
            {
                memcpy(data_prev + cols, data + cols, cols*rows*sizeof(unsigned char));

                bitwise_and(data_prev, cols, end);

                // at this point we have extracted the info of the cells and 
                // we are now ready to print as usual.

                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, header, header_size, data_prev, rows, cols);
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

        const int end = (rows+1)*cols;
        const int augmented_rows = rows + 2;

        data = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));
        data_prev = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));

        if(
                data == NULL
                || data_prev == NULL
                )
        {
            printf("Allocation failed.");
            exit(0);
        }

        // initialize the halo regions to being DEAD
        for(int j = 0; j<cols; ++j)
        {
            DATA(0,j) = DATA(rows + 1,j) = DEAD;
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

        if(fh == NULL)
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
        char * snapshot_name = malloc(snprintf(NULL, 0, "snapshot_%05d", 0)+1);
        if(snapshot_name == NULL)
        {
            printf("Not enough space.");
            exit(0);
        }

        for(int t = 1; t < n+1; ++t)
        {
            
            // We copy the bottom row into the top halo cell. 
            memcpy(data, data + rows*cols, cols*sizeof(unsigned char));

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
            // row 2 to n-1, swap, and row n)
            
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
            memcpy(data + rows*cols, data, cols*sizeof(unsigned char));

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
            memcpy(data+rows*cols+cols, data+cols, cols*sizeof(unsigned char));

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
            memcpy(data + cols, data+rows*cols+cols, cols*sizeof(unsigned char));

            if(t%s == 0)
            {
                memcpy(data_prev + cols, data + cols, cols*rows*sizeof(unsigned char));

                bitwise_and(data_prev, cols, end);

                // at this point we have extracted the info of the cells and 
                // we are now ready to print as usual.

                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, header, header_size, data_prev, rows, cols);
            }

            // data_prev will contain the updated cells (ready for a new 
            // generation) and will begin the next iteration in the for loop 
            // by swapping the halo cells
        }
    
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

