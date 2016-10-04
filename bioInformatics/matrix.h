/*
 * int matrix header file
 * Spring 2015
 * Author: William Willie Wells
 */

#ifndef __MATRIX_H__
#define __MATRIX_H__

// int matrix class 
template<class t>
class matrix{
   public:
      matrix(int row=1,int col=1);// constructor
      ~matrix();// destructor
      void initM(int row=1,int col=1);// initializer
      //void initNest(int row=1,int col=1);

      t &operator()(int row=1, int col=1) const;// overloaded operator
   private:
      t *mat;// 
      int nc;// about of columns 
};

#endif
