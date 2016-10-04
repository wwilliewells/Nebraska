/*
 * Neighbor Joining header file
 * Bioinformatics, Spring 2015
 * Author: William Willie Wells
 */

#ifndef _JOINNEIGHBOR_H_
#define _JOINNEIGHBOR_H_

#include <cstdio>
#include "matrix.h"

// neighbor joining class
class joinNear{
   public:
      joinNear(int,int); // constructor
      ~joinNear(); // destructor

      void joinTree(matrix<char>&);// perform neighbor joining
   private:
      void initJoinNear(int,int);// initializer
      // methods used by joinTree
      void initDistance(matrix<char>&);// start distance
      void minDist();// find leaves with minimum distance
      void join();// join neighbors
      void depth();// compute depth of each leaf/branch
      void sort();// sort "tree"
      int search(int);// search unsorted tree
      void print(matrix<char>&);// print in Newick format

      // variables used to facilitate neighbor joining
      matrix<int> tree; // output tree
      matrix<double> distance;// main distance matrix
      double *branchLen;// output branch weights
      double *residue;// normalized distance to remaining leaves
      int *leafID;// leaf identification, number is name (*,*)
      int leaves;// number of unjoined sequences, decreasing row value
      int neighbor[2];// neighbors to be joined
      int branch;// internal and terminal leaf node
      int *sorted;// sorted "tree" for printing
      int uroot;// unrooted parent value of final leaves
      int nodes;// number of internal + terminal leaves in the "tree"
      int aminoacids;// static column value
      bool newSubTree;// subtree flag
      int scount;
};

#endif
