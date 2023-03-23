#include <stdio.h>
const unsigned int char_size = sizeof(char)*8;
const unsigned int rows = 5;
const unsigned int cols = 5;


// https://stackoverflow.com/questions/20494771/sweeping-through-a-2d-arrays-using-pointers-with-boundary-conditions

void SetBit(char A[],  int k )
{
    int i = k/32;        //gives the corresponding index in the array A
    int pos = k%32;      //gives the corresponding bit position in A[i]

    unsigned int flag = 1;   // flag = 0000.....00001

    flag = flag << pos;      // flag = 0000...010...000   (shifted k positions)

    A[i] = A[i] | flag;      // Set the bit at the k-th position in A[i]
}

char getElement(int* g, int i, int j){
    i = i%(int)rows + rows*(i<0);
    j = j%(int)cols + cols*(j<0);

    return g[i*rows + j];
}

void setElement(int* g, int i, int j, int val){
    g[i*rows + j] = val;

    return ;
}
int get_status_of_neighbors(int *g, int i, int j){
    // need to handle period conditions still
    int s1 = getElement(g, i-1, j-1);
    int s2 = getElement(g, i-1, j  );
    int s3 = getElement(g, i-1, j+1);
    int s4 = getElement(g, i  , j-1);
    int s5 = getElement(g, i  , j+1);
    int s6 = getElement(g, i+1, j-1);
    int s7 = getElement(g, i+1, j  );
    int s8 = getElement(g, i+1, j+1);

    int sum = s1+s2+s3+s4+s5+s6+s7+s8;


    return sum;

}

void upgrade_cell(int *g1, int *g2, int i, int j){
    
    int n_alive_cells = get_status_of_neighbors(g1, i, j);

    if (n_alive_cells >= 2 && n_alive_cells <=3){
        printf("so im alive! ");
        printf("g2 before: %i", g2[i*rows+j]);
        setElement(g2,i,j,1);
        printf("g2 before: %i\n", g2[i*rows+j]);
    }
    else{
        printf("so ops, im dead. ");
        printf("g2 before: %i", g2[i*rows+j]);
        setElement(g2,i,j,0);
        printf("g2 before: %i\n", g2[i*rows+j]);

    }

}

int main(){
    int g1[] = {
        1,0,0,1,0,
        0,1,0,1,1,
        0,0,1,0,0,
        1,0,1,1,1,
        0,0,0,0,1
    };
    int g2[rows*cols];

    for (int i=0; i<rows; i++){
        for (int j=0; j<cols; j++){
            upgrade_cell(g1,g2,i,j);
        }
    }
    printf("\n");
    for (int i=0; i<rows; i++){
        printf("[");
        for (int j=0; j<cols; j++){
            printf("%i,",g2[i*rows+j]);
        }
        printf("]\n");

    }
    
    return 0;
}
