#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>


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

    if(!fh)
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

void display_args()
{
    printf("action (i : 1\tr : 2) -- %d\nk (size) -- %d\ne (0 : ORDERED\t1 : STATIC) -- %d\nf (filename) -- %s\nn (steps) -- %d\ns (save frequency) -- %d\n", action, k, e, fname, n, s );
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

        // we add two rows for halo regions
        const int augmented_rows = rows + 2;

        data = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));

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

        const int augmented_rows = rows + 2;

        data = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));
        data_prev = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));

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

        // even though we have limited n and s, compiler 
        // complains about the size. To prevent this, we allocate extra space
        // char * snapshot_name = malloc(snprintf(NULL, 0, "snapshot_%05d", 0)+1);
        // "snapshot_" requires 9 bytes. 
        // an int goes up to billions plus negative sign so it needs 11 bytes 
        // we also need one byte for end of line character. 
        // in total, that is 21 bytes
        // lets make it an even 32

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

        clock_gettime(CLOCK_MONOTONIC, &start_time);
        for(int t = 1; t < n+1; ++t)
        {
            memcpy(data_prev, data_prev + rows_x_cols, row_len_bytes);
            memcpy(data_prev + rows_x_cols_p_cols, data_prev+cols, row_len_bytes);

            for(int row = 1; row < rows+1; ++row)
            {
                for(int col = 0; col < cols; ++col)
                {
                    upgrade_cell_static(data_prev, data, row, col);
                }
            }
            ++save_counter;

            if(s>0 && save_counter == s && t<100000)
            {
                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, header, header_size, data, rows, cols);
                save_counter = 0;
            }

            tmp_data = data;
            data = data_prev;
            data_prev = tmp_data;
            tmp_data = NULL;
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

        const int augmented_rows = rows + 2;

        data = (unsigned char *) malloc( augmented_rows * cols * sizeof(unsigned char));

        if(!data)
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

        clock_gettime(CLOCK_MONOTONIC, &start_time);
        for(int t = 1; t < n+1; ++t)
        {

            // first we copy the bottom row into the top halo cell
            memcpy(data, data + rows_x_cols, row_len_bytes);

            // then we process all cells starting from row 1 to row 
            // rows - 1. 
            // We cant process row rows because we need to copy the updated 
            // row 1 into the bottom halo.
            for(int row = 1; row < rows; ++row)
            {
                for(int col = 0; col < cols; ++col)
                {
                    upgrade_cell_ordered(data, row, col);
                }
            }

            // we copy row 1 into the bottom halo
            memcpy(data + rows_x_cols_p_cols, data+cols, row_len_bytes);

            // we update the last row now that we have the updated information.
            for(int col = 0; col < cols; ++col)
            {
                upgrade_cell_ordered(data, rows, col);
            }
            ++save_counter;

            if(s>0 && save_counter == s && t<100000)
            {
                sprintf(snapshot_name, "snapshot_%05d", t);
                save_grid(snapshot_name, header, header_size, data, rows, cols);
                save_counter = 0;
            }

        }
        clock_gettime(CLOCK_MONOTONIC, &end_time);
        elapsed = (double)diff(start_time,end_time).tv_sec + (double)diff(start_time,end_time).tv_nsec / 1000000000.0;
        printf("time: %lf\n", elapsed);
    
        free(snapshot_name);
        free(header);
        free(data);
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

