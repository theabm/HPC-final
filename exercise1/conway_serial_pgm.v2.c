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
    // suggest prefetching
    DATA(i,j) = DEAD;

    // this makes the column index periodic (without branches) 
    // for example, if j = 0, j-1 = -1, and -1%cols = -1. The last term 
    // is true, so we obtain cols - 1 
    // if j = col - 1 then j+1 = col (out of bounds). col%col = 0 and the 
    // second term is 0 so we obtain that the new coordinate is 0.
    register int jm1 = (j-1)%(int)cols + cols*((j-1)<0);
    register int jp1 = (j+1)%(int)cols + cols*((j+1)<0);
    register int im1 = i-1;
    register int ip1 = i+1;

    // since we are processing elements sequentially, apart from boundary cells, 
    // the indexes with j - 1 are in potentially loaded caches already since 
    // we processed (i, j-1) previously. 
    // So we have two options
    // 1. (i, j-1) are in the same line as (i,j). In this case, all the elements 
    // with j-1 are already in cache (a cache is surely longer than 2 bytes)
    // So accessing these elements before is beneficial.
    // This is the best case scenario. If we assume a cache line is 64 bytes, 
    // then, we load 64 cells in a line. For 62 cells, the access pattern is HHHHHHHH = 100%. 
    // while for the first cell, the access pattern is MMMHHHHH = 62%
    // 2. (i, j-1) and (i, j) are in different cache lines. Then, since (i, j-1) 
    // was processed previously, we are sure that (i-1, j-1), (i, j-1), and (i+1, j-1)
    // are in the cache. So, in case the cache is full, and we need to flush these 
    // lines out, we access them before.
    // This is the worst case scenario.
    // The access pattern is HHHMMHMH which is 5/8 = 62%
    //
    // So, every 64 bytes, we (62*100% + 2*62%)/64 = 98% on average.
    //
    //
    // Unfortunately, corner cases and boundary cells bring this number down.
    // 
    // In the four corners, the access pattern is MMHMMMMH = 25%. 
    // In the border cases, the situation is similar as explained above. 
    // i.e either 100% or 62%.
    
    unsigned char tmp0=DATA_PREV(im1, jm1);
    unsigned char tmp3=DATA_PREV(i,jm1);
    unsigned char tmp5=DATA_PREV(ip1,jm1);

    unsigned char tmp1=DATA_PREV(im1,j);
    unsigned char tmp2=DATA_PREV(im1,jp1);

    unsigned char tmp4=DATA_PREV(i,jp1);

    unsigned char tmp6=DATA_PREV(ip1,j);
    unsigned char tmp7=DATA_PREV(ip1,jp1);

    register unsigned char n_alive_cells = tmp0+tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7;

    // the majority of cells will be dead and will stay dead 
    // so by reordering the conditions, we enhance branch prediction
    DATA(i,j) = ALIVE*(n_alive_cells==3) + DATA_PREV(i,j)*(n_alive_cells==2);

}

void upgrade_cell_ordered(unsigned char * data, int i, int j)
{
    register int jm1 = (j-1)%(int)cols + cols*((j-1)<0);
    register int jp1 = (j+1)%(int)cols + cols*((j+1)<0);
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

void display_args(){
    printf("action (i : 1\tr : 2) -- %d\nk (size) -- %d\ne (0 : ORDERED\t1 : STATIC) -- %d\nf (filename) -- %s\nn (steps) -- %d\ns (save frequency) -- %d\n", action, k, e, fname, n, s );
}

int main(int argc, char **argv){

    unsigned char * data, * data_prev;

    get_args(argc, argv);
    
    if(action == INIT){
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
        for(int j = 0; j<cols; ++j){
            DATA(0,j) = DATA(rows + 1,j) = DEAD;
        }

        // HEADER INFO CALCULATION

        // intuition: https://stackoverflow.com/questions/3919995/determining-sprintf-buffer-size-whats-the-standard
        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = (char *) malloc(header_size + 1);

        if(!header){

            printf("Allocation failed.");
            exit(0);

        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);

        srand48(10); 

        for(int i = 1; i<rows+1; ++i){
            for(int j = 0; j<cols; ++j)
                DATA(i,j) = drand48() > 0.5 ? ALIVE : DEAD ;
        }

        save_grid(fname, header, header_size, data, rows, cols);
        free(header);
        free(data);
    }
    else if (e == STATIC && action == RUN){
        // reading is done in the same way.

        int opt_args[2] = {0,0};


        FILE * fh_posix = fopen(fname, "r");

        // we know that the magic number is P5 so we set the offset as the  
        // length of P5 * sizeof(char)
        if(fh_posix == NULL){
            fprintf(stderr, "Error opening %s\n", fname);
            exit(0);
        }

        int args_scanned = fscanf(fh_posix, "P5 %d %d 1\n",opt_args, opt_args+1 );
        if(args_scanned != 2){
            printf("fscanf failed.");
            exit(0);
        }

        fclose(fh_posix);

        rows = opt_args[0];
        cols = opt_args[1];

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
        for(int j = 0; j<cols; ++j){
            DATA(0,j) = DATA(rows + 1,j) = DATA_PREV(0,j)
                = DATA_PREV(rows + 1,j) = DEAD;
        }

        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);

        if(!header){

            printf("Allocation failed.");
            exit(0);

        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);

        FILE * fh = fopen(fname, "rb");

        if(fh == NULL){
            fprintf(stderr, "Error opening %s\n", fname);
            exit(0);
        }

        fseek(fh, header_size, 1);

        size_t args_read = fread(data_prev+cols, sizeof(unsigned char), rows*cols, fh);

        if(args_read != (size_t)rows*cols){
            printf("fread failed.");
            exit(0);
        }

        fclose(fh);


        // at this point, we have read the matrix, however, we need to do 
        // some preprocessing to get into the format we described. 
        // As is the case with most optimizations, there is a price to pay, 
        // i.e there is no free lunch.
        
        // step one is to copy the halo rows
        for(int col=0; col<cols;++col){
            DATA_PREV(0,col) = DATA_PREV(rows,col);
            DATA_PREV(rows+1,col) = DATA_PREV(1,col);
        }

        // at the moment, for each cell we have to count the 8 neighbors 
        // encode them in bits, and then do some masking.
        register int jm1, jp1, im1, ip1;
        unsigned char tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
        unsigned char n_alive_cells;
        for(int i = 1; i < rows+1; ++i){
            for(int j = 0; j < cols; ++j){
                jm1 = (j-1)%(int)cols + cols*((j-1)<0);
                jp1 = (j+1)%(int)cols + cols*((j+1)<0);
                im1 = i-1;
                ip1 = i+1;
                tmp0=DATA_PREV(im1, jm1);
                tmp3=DATA_PREV(i,jm1);
                tmp5=DATA_PREV(ip1,jm1);

                tmp1=DATA_PREV(im1,j);
                tmp2=DATA_PREV(im1,jp1);

                tmp4=DATA_PREV(i,jp1);

                tmp6=DATA_PREV(ip1,j);
                tmp7=DATA_PREV(ip1,jp1);

                // this will be a number between 8 and 0
                n_alive_cells = tmp0+tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7;

                // so now I need to shift everything left by one
                n_alive_cells = n_alive_cells<<1;
                
                // lastly, if the cell is dead, we leave the last bit as is 
                // so we do nothing, 
                // and if the cell is alive, we need to set the last bit to 1
                // keeping everything else the same. 
                // to do this, we need to do bitwise OR with DATA_PREV
                DATA_PREV(i,j) |= n_alive_cells;
            }
        }

        // now we have preprocessed the grid to our desired format.
        char * snapshot_name = malloc(snprintf(NULL, 0, "snapshot_%05d", 0)+1);
        if(snapshot_name == NULL){
            printf("Not enough space.");
            exit(0);
        }

        unsigned char *tmp_data = NULL;
        unsigned char *cell_ptr = NULL;


        unsigned char count_of_neighbors = 0;

        register unsigned char val;
        register unsigned char count; 

        for(int t = 1; t < n+1; ++t){
            
            // I have done this already in preprocessing phase
            //
            // for(int col=0; col<cols;++col){
            //     DATA_PREV(0,col) = DATA_PREV(rows,col);
            //     DATA_PREV(rows+1,col) = DATA_PREV(1,col);
            // }
            //
            memcpy(data, data_prev, cols*augmented_rows*sizeof(unsigned char));
            
            cell_ptr = data;

            for(int i = 1; i < rows+1; ++i)
            {
                for(int j = 0; j < cols; ++j)
                {
                    if(*cell_ptr==0)
                    {
                        continue;
                    }
                    // if we get here we have non zero element
                    // we only need to keep track of events that 
                    // change the status of the cell 
                    // i.e births and deaths
                    
                    count_of_neighbors = *cell_ptr >> 1;
                    

                    // we first deal with deaths, however, we need 
                    // to make sure that we only update in the case 
                    // that the cell was alive and died 
                    if(count_of_neighbors!=3 && count_of_neighbors!=2 && *cell_ptr & 0x01)
                    {

                        // to set cell to dead we need to
                        // 1. set cell of data to contents of neighbor count 
                        // with last bit set to zero
                        *cell_ptr = count_of_neighbors<<1; 

                        // 2. get neighbors and decrease the neighbor count 
                        // of each one. We need to be careful, because some 
                        // values have already been copied in the array.
                        
                        jm1 = (j-1)%(int)cols + cols*((j-1)<0);
                        jp1 = (j+1)%(int)cols + cols*((j+1)<0);
                        im1 = i-1;
                        ip1 = i+1;

                        tmp0=DATA(im1, jm1);

                        // we need the state of the last bit.
                        // to get this, we do bit wise and with 0x01
                        val = tmp0 & 0x01;
                        // we then need to decrease neighbor count by one
                        // so we shift the contents right, subtract one, 
                        // and reshift the result one bit left
                        count = ((tmp0>>1) - 0x01)<<1;
                        // finally, we need to join the results. This is done by 
                        // bitwise OR
                        DATA(im1, jm1) = count | val;

                        tmp1=DATA(im1,j);

                        val = tmp1 & 0x01;
                        count = ((tmp1>>1) - 0x01)<<1;
                        DATA(im1, j) = count | val;

                        tmp2=DATA(im1,jp1);

                        val = tmp2 & 0x01;
                        count = ((tmp2>>1) - 0x01)<<1;
                        DATA(im1, jp1) = count | val;

                        tmp3=DATA(i,jm1);

                        val = tmp3 & 0x01;
                        count = ((tmp3>>1) - 0x01)<<1;
                        DATA(i, jm1) = count | val;

                        tmp4=DATA(i,jp1);

                        val = tmp4 & 0x01;
                        count = ((tmp4>>1) - 0x01)<<1;
                        DATA(i, jp1) = count | val;

                        tmp5=DATA(ip1,jm1);

                        val = tmp5 & 0x01;
                        count = ((tmp5>>1) - 0x01)<<1;
                        DATA(ip1, jm1) = count | val;

                        tmp6=DATA(ip1,j);

                        val = tmp6 & 0x01;
                        count = ((tmp6>>1) - 0x01)<<1;
                        DATA(ip1, j) = count | val;

                        tmp7=DATA(ip1,jp1);

                        val = tmp7 & 0x01;
                        count = ((tmp7>>1) - 0x01)<<1;
                        DATA(ip1, jp1) = count | val;
                    }
                    else if(count_of_neighbors == 3 && !(*cell_ptr & 0x01) )
                    {
                        // set to alive

                        // First, we need to set the state of the cell in 
                        // data to alive and have the correct neighbors
                        // this can be done with bitwise or between the 
                        // previous state, and 0x01
                        *cell_ptr |= 0x01; 

                        // 2. get neighbors and increase the neighbor count 
                        // of each one. 

                        jm1 = (j-1)%(int)cols + cols*((j-1)<0);
                        jp1 = (j+1)%(int)cols + cols*((j+1)<0);
                        im1 = i-1;
                        ip1 = i+1;

                        tmp0=DATA(im1, jm1);

                        val = tmp0 & 0x01;
                        count = ((tmp0>>1) + 0x01)<<1;
                        DATA(im1, jm1) = count | val;

                        tmp1=DATA(im1,j);

                        val = tmp1 & 0x01;
                        count = ((tmp1>>1) + 0x01)<<1;
                        DATA(im1, j) = count | val;

                        tmp2=DATA(im1,jp1);

                        val = tmp2 & 0x01;
                        count = ((tmp2>>1) + 0x01)<<1;
                        DATA(im1, jp1) = count | val;

                        tmp3=DATA(i,jm1);

                        val = tmp3 & 0x01;
                        count = ((tmp3>>1) + 0x01)<<1;
                        DATA(i, jm1) = count | val;

                        tmp4=DATA(i,jp1);

                        val = tmp4 & 0x01;
                        count = ((tmp4>>1) + 0x01)<<1;
                        DATA(i, jp1) = count | val;

                        tmp5=DATA(ip1,jm1);

                        val = tmp5 & 0x01;
                        count = ((tmp5>>1) + 0x01)<<1;
                        DATA(ip1, jm1) = count | val;

                        tmp6=DATA(ip1,j);

                        val = tmp6 & 0x01;
                        count = ((tmp6>>1) + 0x01)<<1;
                        DATA(ip1, j) = count | val;

                        tmp7=DATA(ip1,jp1);

                        val = tmp7 & 0x01;
                        count = ((tmp7>>1) + 0x01)<<1;
                        DATA(ip1, jp1) = count | val;
                        
                    }

                }
                cell_ptr++;
            }

            // at the end of this, data contains the new generation.
            // except for the halo rows, which are from the old generation.
            
            // now comes the (ugly) of writing to file.
            
            // To do this, we will need to copy all the contents of 
            // data into data_prev


            if(t%s == 0){
                memcpy(data_prev + cols, data + cols, cols*rows*sizeof(unsigned char));

                // now we can exploit vectorization to do bit AND with 0x01 to 
                // get only last bit.

                // we create an alias of data_prev with vectorized data type 
                // __m256 is a 256 long bit register, aka it holds 32 chars (cells)
                // __m256 *d_p = (__m256*)data_prev
                // const __m256 bits = _mm256_set_epi8(0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01)

                // const int n = (rows*cols)/vector_length;
                // const int remainder = (rows*cols)%vector_length;

                for(int k=0; k<rows*cols; ++k){
                    *(data_prev+k) &= 0x01;
                }
                // at this point we have extracted the info of the cells and 
                // we are now ready to print as usual.

                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, header, header_size, data, rows, cols);
            }

            // finally, after printing (or not) the state is as follows: 
            // data has updated new generation (except halo cells), 
            // data_prev is in an unknown state (either prev gen or has gotten 
            // bit wise AND'ed).
            //
            // so we swap halos (to get ready for the next generation). 


            for(int col=0; col<cols;++col){
                DATA(0,col) = DATA(rows,col);
                DATA(rows+1,col) = DATA(1,col);
            }
            // now data has all the new halos, and is ready for a new generation. 
            // first, we need to swap data and data prev
            tmp_data = data;
            data = data_prev;
            data_prev = tmp_data;
            tmp_data = NULL;
        }
    
        free(snapshot_name);
        free(header);
        free(data);
        free(data_prev);
    }
    else if(e == ORDERED && action == RUN){

        int opt_args[2] = {0,0};


        FILE * fh_posix = fopen(fname, "r");

        // we know that the magic number is P5 so we set the offset as the  
        // length of P5 * sizeof(char)
        if(fh_posix == NULL){
            fprintf(stderr, "Error opening %s\n", fname);
            exit(0);
        }

        int args_scanned = fscanf(fh_posix, "P5 %d %d 1\n",opt_args, opt_args+1 );
        if(args_scanned != 2){
            printf("fscanf failed.");
            exit(0);
        }

        fclose(fh_posix);

        rows = opt_args[0];
        cols = opt_args[1];

        const int augmented_rows = rows + 2;

        data = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));

        if(
                data == NULL
                )
        {
            printf("Allocation failed.");
            exit(0);
        }

        // initialize the halo regions to being DEAD
        for(int j = 0; j<cols; ++j){
            DATA(0,j) = DATA(rows + 1,j) = DEAD;
        }

        int header_size = snprintf(NULL, 0, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);
        char * header = malloc(header_size + 1);

        if(!header){

            printf("Allocation failed.");
            exit(0);

        }

        sprintf(header, HEADER_FORMAT_STRING, rows, cols, MAX_VAL);

        FILE * fh = fopen(fname, "rb");

        if(fh == NULL){
            fprintf(stderr, "Error opening %s\n", fname);
            exit(0);
        }

        fseek(fh, header_size, 1);

        size_t args_read = fread(data+cols, sizeof(unsigned char), rows*cols, fh);

        if(args_read != (size_t)rows*cols){
            printf("fread failed.");
            exit(0);
        }

        fclose(fh);

        char * snapshot_name = malloc(snprintf(NULL, 0, "snapshot_%05d", 0)+1);
        if(snapshot_name == NULL){
            printf("Not enough space.");
            exit(0);
        }

        for(int t = 1; t < n+1; ++t){

            // first we copy the bottom row into the top halo cell
            for(int col=0; col<cols;++col){
                DATA(0,col) = DATA(rows,col);
            }

            // then we process all cells starting from row 1 to row 
            // rows - 1. 
            // We cant process row rows because we need to copy the updated 
            // row 1 into the bottom halo.
            for(int row = 1; row < rows; ++row){
                for(int col = 0; col < cols; ++col){
                    upgrade_cell_ordered(data, row, col);
                }
            }

            // we copy row 1 into the bottom halo
            for(int col=0; col<cols;++col){
                DATA(rows+1,col) = DATA(1,col);
            }

            // we update the last row now that we have the updated information.
            for(int col = 0; col < cols; ++col){
                upgrade_cell_ordered(data, rows, col);
            }

            if(t%s == 0){
                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, header, header_size, data, rows, cols);
            }

        }
    
        free(snapshot_name);
        free(header);
        free(data);
    }
    else {
        printf("Unknown action. Abort");
        exit(0);
    }

    if (fname != NULL )
      free ( fname );

    return 0;
}

