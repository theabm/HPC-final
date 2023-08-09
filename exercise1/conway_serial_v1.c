#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>


#define INIT 1
#define RUN  2

#define K_DFLT 100

#define ORDERED 0
#define STATIC  1


// DEFAULT VALUES
//
char fname_deflt[] = "game_of_life.pgm";

int   action = 0;
int   k      = K_DFLT;
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
                k = atoi(optarg); break;

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
    if ( fname == NULL ){
        fname = (char*)malloc( sizeof(fname_deflt)+1 );
        sprintf(fname, "%s", fname_deflt );
    }
}

void display_args(){
    printf("action (i : 1\tr : 2) -- %d\nk (size) -- %d\ne (0 : ORDERED\t1 : STATIC) -- %d\nf (filename) -- %s\nn (steps) -- %d\ns (save frequency) -- %d\n", action, k, e, fname, n, s );
}

int main(int argc, char **argv){

    get_args(argc, argv);
    display_args();

    if ( fname != NULL )
      free ( fname );

    return 0;
}
