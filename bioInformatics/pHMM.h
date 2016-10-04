/*
 * Profile Hidden Markov Model header file
 * Bioinformatics, Spring 2015
 * Author: William Willie Wells
 */

#ifndef _PHMM_H_
#define _PHMM_H_

#include <cstdio>
#include "matrix.h"

// profile hidden Markov Model class
class pHMM{
   public:
      pHMM(FILE*); // constructor
      ~pHMM(); // destructor
      void initPHMM(FILE*); // initializer

      void search(matrix<char>&,int,int);// searches for related sequences
   private:
      // methods used by initializer
      void initCMfromMSF(FILE*);
      void matchStates();
      void initModel();  
      void updateParameters();
      void countToProb();
      // methods used by search
      double forward(matrix<char>&,int,int);
      void viterbi(matrix<char>&,int,int,FILE*);
      void initPE(int,bool);

      // variables used to facilitate alignment
      // training data
      matrix<char> trainData; // rc[0] * rc[1]
      // probabilities (parameters)
      matrix<double> matchEmit;// 20 * matches
      matrix<double> matchTransit; // 3 * matches +1
      matrix<double> insertEmit; // 20 * matches +1
      matrix<double> insertTransit; // 3 * matches +1
      matrix<double> deleteTransit; // 3 * matches
      matrix<double> existence;// forward
      matrix<double> path;// Viterbi
      bool *matchCol;
      int matches;
      int count;
      int rc[2];
      char aminoAcid;
};

#endif
