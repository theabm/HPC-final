#include <stdio.h>
const unsigned int char_size = sizeof(char)*8;
const unsigned int rows = 5;
const unsigned int cols = 5;


// https://stackoverflow.com/questions/20494771/sweeping-through-a-2d-arrays-using-pointers-with-boundary-conditions

char check_if_nth_bit_is_set(char* g, int pos, int n)
{
    // Pretty self explanatory, returns 1 if bit is set and zero otherwise
    unsigned char flag = 1;   // flag = 0000.....00001
    flag = flag << n;
    return g[pos] & flag;      
}

char getElement(char* g, int i, int j)
{
    // This function takes as input the coordinates of the neighboring
    // cells. Therefore, this function accounts for boundary conditions in 
    // case that coordinates are negative or too big. 
    // To illustrate, an example where rows = 5
    // If the coordinates are negative, -1%rows = -1 (if rows is uint it is 0 )
    // so we need to cast to int for correct behavior)
    // -1 + 5 = 4 which is at the boundary on the other side
    // for positive coordinates, the second term is zero and we just use a mod 
    // to go to proper place
    i = i%(int)rows + rows*(i<0);
    j = j%(int)cols + cols*(j<0);

    return check_if_nth_bit_is_set(g, i*rows + j, 0);
}

void set_nth_bit(char* g, int pos, int n)
{
    // Sets the n-th of the array element at pos to 1.
    // Therefore, if we have a sequence 10110001 and we want 
    // to set the 1st (counting from zero) bit to 1, we 
    // shift 00000001 left by one -> 00000010 and do a bitwise
    // or which outputs 10110011 
    //
    unsigned char flag = 1;   // flag = 0000.....00001
    flag = flag << n;
    g[pos] = g[pos] | flag;

    return ;
}

void clear_nth_bit(char* g, int i, int n)
{
    // sets the nth bit of array element i to zero. 
    // as a quick example: we want to set the 4th bit of 
    // 10110001 to zero (counting from zero on the right)
    // we shift 1 to the left 4 times and negate to obtain
    // 11101111 and then do bitwise and to obtain
    // 10100001 which is just what we want.

    unsigned char flag = 1;
    flag = ~(flag << n);
    g[i] = g[i] & flag;

    return ;
}
char get_status_of_neighbors(char *g, int i, int j)
{
    char s1 = getElement(g, i-1, j-1);
    char s2 = getElement(g, i-1, j  );
    char s3 = getElement(g, i-1, j+1);
    char s4 = getElement(g, i  , j-1);
    char s5 = getElement(g, i  , j+1);
    char s6 = getElement(g, i+1, j-1);
    char s7 = getElement(g, i+1, j  );
    char s8 = getElement(g, i+1, j+1);

    return s1+s2+s3+s4+s5+s6+s7+s8;
}

void upgrade_cell(char *g1, int i, int j)
{
    int n_alive_cells = get_status_of_neighbors(g1, i, j);

    if (n_alive_cells >= 2 && n_alive_cells <=3){
        set_nth_bit(g1, i*rows + j, 1);
    }
    else{
        clear_nth_bit(g1, i*rows + j, 1);
    }
}
void display_grid(char *g)
{
    for(int i=0; i<rows; i++){
        printf("[ ");
        for(int j=0; j<cols; j++){
            for(int k=0; k<char_size; k++){
                printf("%d",!!((g[i*rows + j]<<k) & 0x80));
            }
            printf(" , ");
        }
        printf(" ]\n");
    }
    return;
}

int main()
{
    char g1[] = {
        1,0,0,1,0,
        0,1,0,1,1,
        0,0,1,0,0,
        1,0,1,1,1,
        0,0,0,0,1
    };
    printf("Grid of life before:\n");
    display_grid(g1);

    for (int i=0; i<rows; i++){
        for (int j=0; j<cols; j++){
            upgrade_cell(g1,i,j);
        }
    }
    
    printf("\nGrid of life after:\n");
    display_grid(g1);
    return 0;
}
