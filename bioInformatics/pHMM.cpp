/*
 * Profile Hidden Markov Model -- pHMM.cpp
 * Bioinformatics, Spring 2015
 * Author: William Willie Wells
 */

#include <cstdio> // using fopen,fclose,fscanf,getc,rewind,sprintf,fputs,perror,NULL, EOF, FILE
#include <iostream> // using cout
#include "pHMM.h"
//#include "matrix.h"
#include <cmath> // using exp, log

using namespace std;

const int aa = 20;

// lookup table for score matrix indices 
int aaToN(char c){
   switch(c){
      case 'A':return 0;
      case 'R':return 1;
      case 'N':return 2;
      case 'D':return 3;
      case 'C':return 4;
      case 'Q':return 5;
      case 'E':return 6;
      case 'G':return 7;
      case 'H':return 8;
      case 'I':return 9;
      case 'L':return 10;
      case 'K':return 11;
      case 'M':return 12;
      case 'F':return 13;
      case 'P':return 14;
      case 'S':return 15;
      case 'T':return 16;
      case 'W':return 17;
      case 'Y':return 18;
      case 'V':return 19;
      default: return -1;
   }
}

// matrix constructor
template <class t>
matrix<t>::matrix(int r, int c){
   initM(r,c);
}

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

// overloaded operator ()
template <class t>
t &matrix<t>::operator()(int r, int c)const{
    return mat[r*nc+c];
}

// initialize char matrix from a FASTA formatted file
void initCMfromFASTA(FILE *f, matrix<char> &cm, int &r,int &c){
   char x;
   int mcx=0;
   int cx=0;
   int rx=0;
   bool flag=0;

   // count total characters
   do{
      x=getc(f);
      // if '>' is encountered stop count, increment rows, calculate max column
      if(x=='>'){
         if(cx > mcx){ mcx = cx; }
         flag=0;
         rx++;cx=0;
      }// set a flag to start counting after both '>' and '\n' are encountered
      if(x=='\n' && flag==0){
         flag=1;
      }// count charcters until '>' is encountered
      if(flag==1){
         if(((x >= 'a' && x <= 'z') || (x >= 'A' && x <= 'Z'))){
      	   cx++;
         }
      }
   }while(x != EOF);
   // assign referenced row and maximum column values
   r=rx;
   c=mcx;
   // allocate memory for char array of count size and rewind file
   cm.initM(rx,mcx);
   rewind(f);
   // initialize char matrix to characters in FASTA file
   rx=-1;
   flag=0;
   do{
      x=getc(f);
      // if '>' is encountered do not count
      if(x=='>'){
         flag=0;
         rx++;
         cx=0;
      }// set a flag to start counting after both '>' and '\n' are encountered
      if(x=='\n' && flag==0){
         flag=1;
      }// count characters until '>' is encountered
      if(flag==1){
         if(((x >= 'a' && x <= 'z') || (x >= 'A' && x <= 'Z'))){
            cm(rx,cx)=x;
      	    cx++;
         }
      }
   }while(x != EOF);
}

// convert default initialization character to indel character
void convertC(matrix<char>& mc, int r,int c){
    //
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            if(mc(i,j)<'A' || mc(i,j) > 'Z'){
                mc(i,j)='~';
            }
        }
    }
}

// initialize char matrix from a MSF formatted file 
void pHMM::initCMfromMSF(FILE *f){
   char x;
   int nx=0;
   int rx=0;
   int cx=0;
   int ix=0;
   int nc=0;
   bool flag[3]={false,false,false};
   char y=' ';

   // count total characters
   do{
      x=getc(f);
      // set a flag to start counting after both '/'
      if(x=='/'){
         flag[0]=true;
      }// count characters 
      if(flag[0]==true){
         if(x=='.' || x=='~' || (x >= 'A' && x <= 'Z')){
      	    ix++;// total count
         }
         if(x >= 'a' && x <= 'z'){flag[1]=true;}
         if(flag[1]==true && x=='\n' && x!=y && flag[2]==false){
            nx++;// number of rows
         }
         if(x=='\n' && x==y){
            flag[1]=false;
            if(nx>2){ flag[2]=true; }
         }
      }
      y=x;
   }while(x != EOF);
   // allocate memory for char array of count size and rewind file
   trainData.initM(nx,ix/nx);
   rewind(f);
   flag[0]=false;
   flag[1]=false;
   // initialize char matrix to characters in MSF file
   do{
      x=getc(f);
      // set a flag to start counting after both '/'
      if(x=='/'){
         flag[0]=true;
      }// count characters 
      if(flag[0]){
         if(rx==0 && flag[1]==false){ cx=nc; flag[1]=true; }
         if(x=='.' || x=='~' || (x >= 'A' && x <= 'Z')){
            trainData(rx,cx)=x;
            cx++;
         }else if(x=='\n' && cx>nc){
             rx++;
             if(rx<nx){ cx=nc; }
         }
         if(rx==nx){ rx=0; flag[1]=false; nc=cx; }
      }
   }while(x != EOF);
   rc[0]=nx;// assign row and column values of trianing data
   rc[1]=ix/nx;
}

// constructor
pHMM::pHMM(FILE *f){
    initPHMM(f);
}

// destructor
pHMM::~pHMM(){
    delete []matchCol;
}

// initializer
void pHMM::initPHMM(FILE *f){
    // initialize training data
    initCMfromMSF(f);
    // determine amount of match states
    matchStates();
    // create states in model
    initModel();
    // initialize state emmission and transition parameters
    updateParameters();
}

// determine amount of match states
void pHMM::matchStates(){
    // allocate memory for pHMM class boolean array
    matches=0;
    matchCol = (bool *) new bool [rc[1]];
    for(int j=0;j<rc[1];j++){
        count=0;// reset count
        for(int i=0;i<rc[0];i++){
            if(trainData(i,j) >='A' && trainData(i,j) <= 'Z'){
                count++;// vote for match column
            }
        }// if votes > half amount of rows increment matches
        if(count > rc[0]/2){ matches++; matchCol[j]=1; }
        else{ matchCol[j]=0; }
    }
}

// determine how to define the model's architecture
void pHMM::initModel(){
    // initialize emission and transition parameters
    int transit = 3;
    matchEmit.initM(matches,aa);
    insertEmit.initM(matches+1,aa);
    matchTransit.initM(matches+1,transit);
    insertTransit.initM(matches+1,transit);
    deleteTransit.initM(matches+1,transit);

    // use Laplace's rule for prior distribution 
    double prior = 1.0;
    for(int i=0;i<matches+1;i++){
        // initialize emission parameters
        for(int j=0;j<aa;j++){
            if(i<matches){ matchEmit(i,j)=prior; }
            insertEmit(i,j)=prior/aa;
        }
        // initialize transition parameters
        for(int k=0;k<transit;k++){
            if(i == matches && k == 2){
               prior = 0.0;
            }
            matchTransit(i,k)=prior;
            insertTransit(i,k)=prior;
            if(i>0){ deleteTransit(i,k)=prior; }
            if(i==matches && k == 2 || i==0){ deleteTransit(i,k)=0.0; }
        }
    }
}

// convert counts to probabilities
void pHMM::countToProb(){
    // sum total counts and divide individual counts by total
    double totals[4];
    for(int i=0;i<matches+1;i++){
        // reset totals
        totals[1]=0.0;
        totals[2]=0.0;
        totals[3]=0.0;
        if(i<matches){
            totals[0]=0.0;
            // calculate match emission probabilities
            for(int j=0;j<aa;j++){
                totals[0] = totals[0] + matchEmit(i,j);
            }
            for(int m=0;m<aa;m++){
                matchEmit(i,m) = matchEmit(i,m) / totals[0];
            }
        }
        // calculate transition probabilities
        for(int k=0;k<3;k++){
            totals[1] = totals[1] + matchTransit(i,k);
            totals[2] = totals[2] + insertTransit(i,k);
            if(i>0){ totals[3] = totals[3] + deleteTransit(i,k); }
        }
        for(int r=0;r<3;r++){
            matchTransit(i,r) = matchTransit(i,r) / totals[1];
            insertTransit(i,r) = insertTransit(i,r) / totals[2];
            if(i>0){ deleteTransit(i,r) = deleteTransit(i,r) / totals[3]; }
        }
    }
}

// count amino acid emissions and state transitions
void pHMM::updateParameters(){
    // count
    int pstate;
    int cstate; 
    for(int i=0;i<rc[0];i++){
        count=0;
        pstate=-1;
        for(int j=0;j<rc[1];j++){
            aminoAcid=trainData(i,j);
            // count match state emissions and set transition flag
            if(matchCol[j]==1){
                count++;
     	        if(aminoAcid >= 'A' && aminoAcid <= 'Z'){
                    matchEmit(count-1,aaToN(aminoAcid))++;
                    cstate=0;
                }else{
                    cstate=2;
                }
            }else{
                if(aminoAcid >= 'A' && aminoAcid <= 'Z'){ cstate=1; }
                else{ cstate=-1; } // deal with leading and trailing gaps
            }
            // count transitions
            if(cstate > -1){
                if(pstate==0){
                    if(cstate==1){
                        matchTransit(count,cstate)++; 
                    }else{ 
                        matchTransit(count-1,cstate)++; 
                    }
                }else if(pstate==1){
                    insertTransit(count,cstate)++;
                }else if(pstate==2){
                    if(cstate==1){ deleteTransit(count,cstate)++; }
                    else{ deleteTransit(count-1,cstate)++; }
                }else{
                    if(cstate==1){ matchTransit(count,cstate)++; }
                    else{ matchTransit(count-1,cstate)++; }
                }
            }else{
                if(count == matches){
                    if(pstate==0){
                        matchTransit(count,0)++;
                    }else if(pstate==1){
                        insertTransit(count,0)++;
                    }else if(pstate==2){
                        deleteTransit(count,0)++;
                    }
                }
            }
            // set previous state to current state
            pstate=cstate;
        }
    }
    // convert counts to probabilities
    countToProb();
}

// initialize path (Viterbi) or existence (forward) matrix
void pHMM::initPE(int timesteps, bool flag){
    if(flag==0){ existence.initM(timesteps,3); }
    else{ path.initM(timesteps,3); }
    for(int i=0;i<timesteps;i++){
        for(int j=0;j<3;j++){
            if(flag==0){ existence(i,j)=0.0; }
            else{ path(i,j)=0.0; }
        }
    }
    if(flag==0){ existence(0,0)=1.0; }
    else{ path(0,0)=1.0; }
}

// compute profile Hidden Markov Model forward algorithm
double pHMM::forward(matrix<char> &db,int x,int c){
    double pNew;
    int index;
    count=0;
    // initialize existence matrix
    initPE(c,false);
    // check entire sequence
    while(count<c){
        // check for valid amino acid 
        aminoAcid = db(x,count);
        if(aminoAcid=='~'){ break; } // end of FASTA sequence reached
        index = aaToN(aminoAcid);
        if(index==-1){ aminoAcid=db(x,count+1); index=aaToN(aminoAcid); }
        if(count==0){
            // first match state 
            pNew=exp(existence(count,0))*matchTransit(count,0);
            existence(count+1,0)=log(pNew)+log(aa*matchEmit(count,index));
            // first insert state
            pNew=exp(existence(count,0))*matchTransit(count,1);
            existence(count,1)=log(pNew)+log(aa*insertEmit(count,index));
            // first delete state
            pNew=exp(existence(count,0))*matchTransit(count,2);
            existence(count+1,2)=log(pNew);
        }else if(count>0 && count<matches){
            // next match state
            pNew=exp(existence(count,0))*matchTransit(count,0);
            pNew=pNew+exp(existence(count,1))*insertTransit(count,0);
            if(count > 1){ // columns 0 and 1 have no delete state inputs
                pNew=pNew+exp(existence(count,2))*deleteTransit(count,0);
            }
            existence(count+1,0)=log(pNew)+log(aa*matchEmit(count,index));
            // next insert state
            pNew=exp(existence(count,0))*matchTransit(count,1);
            pNew=pNew+exp(existence(count-1,1))*insertTransit(count,1);
            if(count > 1){
                pNew=pNew+exp(existence(count,2))*deleteTransit(count,1);
            }
            existence(count,1)=log(pNew)+log(aa*insertEmit(count,index));
            // next delete state
            pNew=exp(existence(count,0))*matchTransit(count,2);
            pNew=pNew+exp(existence(count,1))*insertTransit(count,2);
            if(count > 1){    
                pNew=pNew+exp(existence(count,2))*deleteTransit(count,2);
            }    
            existence(count+1,2)=log(pNew);
        }else{ // entering looped final insert state
            if(count==matches){ // first loop
                pNew=exp(existence(count,0))*matchTransit(count,1);
                pNew=pNew+exp(existence(count-1,1))*insertTransit(count,1);
                pNew=pNew+exp(existence(count,2))*deleteTransit(count,1);
                existence(count,1)=log(pNew)+log(aa*insertEmit(count,index));
            }else{ // consecutive loops
                pNew=exp(existence(count-1,1))*insertTransit(matches,1);
                existence(count,1)=log(pNew)+log(aa*insertEmit(matches,index));
            }
        }
        count++;// get next amino acid
    }
    // P(sequence) = sum(0:k,f_k(L)P(0|k), transition to end state
    pNew=exp(existence(count,0))*matchTransit(matches,0);
    pNew=pNew+exp(existence(count,1))*insertTransit(matches,0);
    pNew=pNew+exp(existence(count,2))*deleteTransit(matches,0);
    existence(count+1,0)=log(pNew);

    return existence(count+1,0);// end state
}

// compute most likely path through profile Hidden Markov Model
void pHMM::viterbi(matrix<char>& db,int r,int c,FILE* f){
    double tmp;
    double p;
    int cstate=-1;
    int index;
    int m=0;
    int iter=0;
    bool flag=0;

    // initialize path matrix
    initPE(c,true);
    count=0;
    while(count<c && iter < rc[1]){
        // check for valid amino acid 
        aminoAcid = db(r,count);
        if(aminoAcid=='~'){break;}
        index = aaToN(aminoAcid);
        if(index==-1){count++; flag=1; continue;}
        if(count==0 || flag){
            // first match state 
            p=log(aa*matchEmit(m,index));
            path(count+1,0)=p+path(count,0)+log(matchTransit(m,0));
            tmp=path(count+1,0); 
            cstate=0;
            // first insert state
            p=log(aa*insertEmit(m,index));
            path(count,1)=p+path(count,0)+log(matchTransit(m,1));
            if(path(count,1)>tmp){ tmp=path(count,1); cstate=1; }
            // first delete state
            path(count+1,2)=path(count,0)+log(matchTransit(m,2));
            // print path, start in an emitting state
            if(cstate==0){ m++; count++; fputc(aminoAcid,f); }
            else if(cstate==1){ count++; fputc(aminoAcid,f); }
        }else if(count>0 && m<matches){
            // next match state
            tmp=path(count,0)+log(matchTransit(m,0));
            p=path(count,1)+log(insertTransit(m,0));
            if(p > tmp){ tmp=p; pstate=1; }
            if(m > 1){ 
                p=path(count,2)+log(deleteTransit(m,0));
                if(p > tmp){ tmp=p; pstate=2; }
            }
            // add maximum transition probability to emission probability
            path(count+1,0)=tmp+log(aa*matchEmit(m,index));
            // next insert state
            tmp=path(count-1,0)+log(matchTransit(m,1));
            p=path(count-1,1)+log(insertTransit(m,1));
            if(p > tmp){ tmp=p; pstate=1; }
            if(m > 1){
                p=path(count-1,2)+log(deleteTransit(m,1));
                if(p > tmp){ tmp=p; pstate=2; }
            }
            // add maximum transition probability to emission probability
            path(count,1)=tmp+log(aa*insertEmit(m,index));

            // next delete state
            tmp=path(count,0)+log(matchTransit(m,2));
            p=path(count,1)+log(insertTransit(m,2));
            if(p > tmp){ tmp=p; pstate=1; }
            if(m > 1){
                p=path(count,2)+log(deleteTransit(m,2));
                if(p > tmp){ tmp=p; pstate=2; }  
            }
            path(count+1,2)=tmp;

            // determine current maximum
            tmp=path(count+1,0);
            cstate=0;
            if(path(count,1) > tmp){ tmp=path(count,1); cstate=1; }
            if(path(count+1,2) > tmp){ tmp=path(count+1,2); cstate=2; }
            // print path
            if(cstate==0){ m++; count++; fputc(aminoAcid,f); }
            else if(cstate==1){ count++; fputc(aminoAcid,f); }
            else if(cstate==2){ fputc('~',f); iter++; }
        }else{ // trailing insert states
            if(m==matches){ // first trailing insert state
                tmp=path(count-1,0) + log(matchTransit(m,1));
                p=path(count-1,1)+log(insertTransit(m,1));
                if(p > tmp){ tmp=p; }
                p=path(count-1,2)+log(deleteTransit(m,1));
                if(p > tmp){ tmp=p; }
                path(count,1)=tmp+log(aa*insertEmit(m,index));
            }else{ // consecutive insert states
                p=path(count-1,1)+log(insertTransit(m,1));
                path(count,1)=p+log(aa*insertEmit(m,index));
            }
            fputc(aminoAcid,f);
            count++;
        }
    }// space output
    fputc('\n',f);
    fputc('\n',f);
    // P(sequence) = sum(0:k,v_k(L)P(0|k), end state value
    tmp=path(count,0)+log(matchTransit(matches,0)); 
    p=path(count,1)+log(insertTransit(matches,0));
    if(p > tmp){ tmp=p; }
    p=path(count,2)+log(deleteTransit(matches,0));
    if(p > tmp){ tmp=p; }
    path(count+1,0)=tmp;
}

// search database for proteins related to the model
void pHMM::search(matrix<char>& db,int r,int c){
    double p1,p2;
    FILE *f;
    
    // perform forward algorithm on each sequence
    f=fopen("phmm.txt","w");
    for(int i=0;i<r;i++){
        p1=forward(db,i,c);
        // if > threshold perform Viterbi on sequence
        if(p1 > 0){
            viterbi(db,i,c,f);
        }
    }
    fclose(f);
}

// Main
int main(int argc, char** argv){
   // declare a reused file variable
   FILE *file1;
   
   // parse training file (multiple alignment), and create model
   file1=fopen("train.msf","r");
   if(file1==NULL){perror("file does not exist\n");}
   pHMM phmm(file1);
   fclose(file1);
      
   // read in database file
   file1=fopen("database3.txt","r");
   if(file1==NULL)perror("file does not exist\n");
   matrix<char> basedata;
   int bdRow,bdCol;
   initCMfromFASTA(file1,basedata,bdRow,bdCol);
   convertC(basedata,bdRow,bdCol);
   fclose(file1);

   // search database for sequences related to Profile Hidden Markov Model
   phmm.search(basedata,bdRow,bdCol);
 
   return 0; 
}
