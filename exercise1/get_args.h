#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>


int get_args ( int argc, char **argv )
{
  int action = 0;
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

  return 0;
}

void check_args(){
    printf("action: %s\t k: %d\t e:%d\t f:%s\t n:%d\t s:%d\t", action, k, e, fname, n, s);
    return 0;
}
