/*
 * int matrix
 * Spring 2015
 * Author: William Willie Wells
 */

#include "include/matrix.h"
#include <iostream>

using namespace std;// new,delete

// matrix constructor
template <class t>
matrix<t>::matrix(int r, int c){
   initM(r,c);
}

/*template <class t>
void matrix<t>::initNest(int r, int c){
    mat=(t *) new t [r*c];
}*/

// matrix destructor
template <class t>
matrix<t>::~matrix(){
   delete []mat;
}

// matrix initialization
template <class t>
void matrix<t>::initM(int r,int c){
    mat=(t *) new t [r*c];// allocate memory
   
    for(int i=0;i<r*c;i++){
        mat[i]=(t)-1;
    }
    nc=c;
}

//
template <class t>
t &matrix<t>::operator()(int r, int c)const{
    return mat[r*nc+c];
}

