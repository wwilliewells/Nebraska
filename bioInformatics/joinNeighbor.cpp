/*
 * Neighbor Joining -- joinNeighbor.cpp
 * Bioinformatics, Spring 2015
 * Author: William Willie Wells
 */

#include <cstdio> // using fopen,fclose,getc,rewind,fprintf,perror,NULL,EOF,FILE
#include <iostream> // using cout
#include "include/joinNeighbor.h"
#include "include/matrix.h"
#include <cmath> // using log
#include <cstdlib>//rand
#include <ctime>// using srand. time

using namespace std;

const double aa = 20.0;

// lookup table for random number generator
char nToAa(int i){
   switch(i){
      case 0: return 'A';
      case 1: return 'R';
      case 2: return 'N';
      case 3: return 'D';
      case 4: return 'C';
      case 5: return 'Q';
      case 6: return 'E';
      case 7: return 'G';
      case 8: return 'H';
      case 9: return 'I';
      case 10: return 'L';
      case 11: return 'K';
      case 12: return 'M';
      case 13: return 'F';
      case 14: return 'P';
      case 15: return 'S';
      case 16: return 'T';
      case 17: return 'W';
      case 18: return 'Y';
      case 19: return 'V';
      default: return 'z';
   }
}/**/

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

// constructor
joinNear::joinNear(int r, int c){
    initJoinNear(r,c);
}

// destructor
joinNear::~joinNear(){
    delete []residue;
    delete []branchLen;
    delete []leafID;
    delete []sorted;
}

// initializer
void joinNear::initJoinNear(int r, int c){
    // initialize static (memory)
    leaves=r;
    nodes=2*leaves-2;
    aminoacids=c;
    branch=-1;
    uroot=1;
    newSubTree=true;
    neighbor[0]=-1;
    neighbor[1]=-1;

    // initialize dynamic (memory)
    tree.initM(4,nodes);
    distance.initM(leaves,leaves);
    branchLen = (double*) new double [nodes-1];
    residue = (double*) new double [leaves];
    sorted = (int*) new int [nodes];
    leafID = (int*) new int [leaves];
    for(int i=0;i<leaves;i++){ leafID[i]=i; }
}

// initialize distance using modified Jukes-Cantor
void joinNear::initDistance(matrix<char>& database){ 
    double iFrac=(aa-1)/aa;
    double frac=1/((double)aminoacids + 0.00001);// log(1)=0
    double count;

    // base distance off of amount of similar elements
    for(int i=0;i<leaves;i++){
        for(int j=0;j<leaves;j++){
            count=0.0; // reset count
            for(int k=0;k<aminoacids;k++){
                if(database(i,k) != database(j,k)){
                    count++;// fraction of elements that differ
                }
            }// removed element remaining fraction from log for stability
            distance(i,j) = (-1.0*iFrac)*log(1-count*frac);
        }
    }
}

// find leaves that are closest, minimum distance apart
void joinNear::minDist(){ 
    // calculate residue for each leaf
    double remaining=(double)leaves - 2.0;
    for(int i=0;i<leaves;i++){
        residue[i]=0.0;
        for(int j=0;j<leaves;j++){
            if(i==j){distance(i,j)=0.0; }
            residue[i]+=distance(i,j);
        }
        residue[i]/=remaining;
     }

     // find leaves that are closest
     double D,min=1000.0;
     for(int m=0;m<leaves;m++){
        for(int k=m+1;k<leaves;k++){
            D = distance(m,k) - (residue[m] + residue[k]);
            if(D < min){ min=D; neighbor[0]=m; neighbor[1]=k; }
        }
    }
}

// join leaves that are closest 
void joinNear::join(){
    // copy non-selected leaves into temporary matrix and copy identification
    matrix<double> dNew;
    dNew.initM(leaves-1,leaves-1);
    int ID[leaves];// leaf identification
    int n;
    int m=-1;
    for(int i=0;i<leaves;i++){
        // copy if row is not one of the selected leaves
        if(i != neighbor[0] && i != neighbor[1]){
            n=-1;
            m++;
            ID[m]=leafID[i];
            // copy if column is not one of the selected leaves
            for(int j=0;j<leaves;j++){
                if(j != neighbor[0] && j != neighbor[1]){
                    n++;
                    dNew(m,n)=distance(i,j);
                }
            }
        }
    }
    // create new leaf from selected leaves compute distances and add to matrix 
    n=0;
    for(int k=0;k<leaves;k++){
        if(k != neighbor[0] && k != neighbor[1]){
            dNew(m+1,n) = distance(neighbor[0],k) + distance(neighbor[1],k);
            dNew(m+1,n) = dNew(m+1,n) - distance(neighbor[0],neighbor[1]);
            dNew(m+1,n) = dNew(m+1,n)/2.0;
            dNew(n,m+1) = dNew(m+1,n);
            n++;
        }
    }
    dNew(m+1,m+1)=0.0;
    // zeroing out old matrix for stability (not necessary)
    for(m=0;m<leaves;m++){
        for(n=0;n<leaves;n++){
            distance(m,n)=0.0;
        }
    }
    // write new matrix into old matrix and add identification for the new leaf
    for(m=0;m<leaves-1;m++){
        for(n=0;n<leaves-1;n++){
            distance(m,n) = dNew(m,n);
        }
        if(m < leaves - 2){ leafID[m]=ID[m]; }
        else{ leafID[m] = branch; }
    }
}

// search unsorted "tree" for leaf, direction gives forward and backward   
int joinNear::search(int direction){
    for(int i=0;i<nodes;i++){
        if(tree(direction,i) == branch){
            return i; 
        }
    }
    return -1;// error case, not found
}

// compute depth of leaves for printing '(' characters
void joinNear::depth(){
    for(int i=0;i<nodes;i=i+2){
        branch=tree(1,i);
        tree(2,i)=0;
        while(branch != uroot){
            if(branch != uroot){
                tree(2,i)++;
                branch=tree(1,search(0));
            }
        }cout<<tree(2,i)<<"  ";
        tree(2,i+1)=tree(2,i);// siblings have the same parent
    }cout<<endl;
}

// sort "tree", into ordered progressive "subtrees" 
void joinNear::sort(){
    //int ti;
    int index=0;
    while(scount<nodes){
        branch=tree(0,index);
        if(branch > -1){
            sorted[scount]=branch;
            tree(3,search(0))=1;
            scount++;
            if(tree(0,search(0)+1) > -1){
                sorted[scount]=tree(0,search(0)+1);
                tree(3,search(0)+1)=1;
                scount++;
                branch=tree(1,search(0));
            }else{
                branch=tree(0,search(0)+1);
                if(tree(3,search(1))==1 && tree(3,search(1)+1)==1){
                    sorted[scount]=branch;
                    tree(3,search(0))=1;
                    scount++;
                    branch=tree(1,search(0));
                }else{
                    index=search(1);
                }
            }
        }else{
            if(tree(3,search(1))==1 && tree(3,search(1)+1)==1){
                sorted[scount]=branch;
                tree(3,search(0))=1;
                scount++;
                branch=tree(1,search(0));
            }else{
                index=search(1);
            }
        } 
    }
}
                
// print tree in Newick format
void joinNear::print(matrix<char>& database){
    int count=0;
    int i,j,index=0;
    int open=0;
    char buf[aminoacids+1];
    FILE* f;
    f=fopen("tree1.txt","w");
    fprintf(f,"(");
    // iterate through sorted "tree"
    while(count < nodes-1){cout<<sorted[count]<<' ';
        // get current leaf/branch value and position in unsorted "tree"
        branch = sorted[count];
        index=search(0);
        // determine if a new "subtree" is required and print '(' characters
        if(sorted[count] > -1 && sorted[count+1] > -1){ newSubTree = true; }
        if(newSubTree==true){
            for(j=0;j<tree(2,index)-open;j++){ fprintf(f,"("); }
            open=tree(2,index);
            newSubTree=false;
        }
        // if leaf is internal print ')' otherwise print applicable sequence
        if(branch > -1){
            for(i=0;i<aminoacids;i++){
                buf[i]=database(branch,i);
            }
            buf[aminoacids]='\n';
            fprintf(f,"%s",buf);
        }else{
            fprintf(f,")");
            open--;
        }
        // print corresponding distance value
        if(tree(1,index) != uroot){
            fprintf(f,":%f",branchLen[index]);
        }
        // print commas
        if(count < nodes-2 && sorted[count+1] > -1){ fprintf(f,",\n"); }

        count++;
    }cout<<sorted[count]<<endl;
    // print last unrooted "tree" distance and close file
    fprintf(f,"):%f)\n",branchLen[count-1]);
    fclose(f);
}

// perform neighbor joining and output tree in Newick format
void joinNear::joinTree(matrix<char>& database){
    // calculating starting distance measure
    initDistance(database);
    // join leaves 
    FILE *f;
    f=fopen("tree2.txt","w");
    char buf1[aminoacids];
    char buf2[aminoacids];
    int count=0; 
    while(leaves > 2){
        // determine leaves with the minimal distance
        minDist();

        // add chosen leaves to tree, assign parent node (branch)
        tree(0,count+1)=leafID[neighbor[1]];
        tree(0,count)=leafID[neighbor[0]];
        tree(1,count+1)=branch;
        tree(1,count)=branch;
        // calculate distances that are associated with chosen leaves
        branchLen[count]=distance(neighbor[0],neighbor[1]);
        branchLen[count]=branchLen[count]+residue[neighbor[0]];
        branchLen[count]=branchLen[count]-residue[neighbor[1]];
        branchLen[count]=branchLen[count]/2.0;
        branchLen[count+1]=distance(neighbor[0],neighbor[1]) - branchLen[count];

        // perform joining, "pop off" chosen leaves "push on" branch 
        join();
        sprintf(buf1,"%d",tree(0,count));
        sprintf(buf2,"%d",tree(0,count+1));
        for(int i=0;i<aminoacids;i++){
            if(tree(0,count) > -1){ buf1[i]=database(tree(0,count),i); }
            if(tree(0,count+1) > -1){ buf2[i]=database(tree(0,count+1),i); }
        }
        fprintf(f,"(%s:%f,\n %s:%f):%d, \n",buf1,branchLen[count],buf2,branchLen[count+1],branch);
        // update next: branch indentification, tree position, number of leaves
        branch--;
        count+=2;
        leaves--;
    }
    // add last two leaves to the tree, assign branch which is the anti-root
    tree(0,count)=leafID[0];
    tree(0,count+1)=leafID[1];
    tree(1,count)=branch;
    tree(1,count+1)=branch;
    uroot=branch;// unrooted value, present as a parent but not as a node
    // assign final distance between unchosen leaves
    branchLen[count]=distance(0,1);
    sprintf(buf1,"%d",tree(0,count));
    sprintf(buf2,"%d",tree(0,count+1));
    for(int i=0;i<aminoacids;i++){
        if(tree(0,count) > -1){ buf1[i]=database(tree(0,count),i); }
        if(tree(0,count+1) > -1){ buf2[i]=database(tree(0,count+1),i); }
    }
    fprintf(f,"(%s:%f,\n %s:%f):%d, \n",buf1,branchLen[count],buf2,branchLen[count+1],branch);
        
    // compute depth of each leaf/branch, then create a sorted tree, and print
    depth();
    branch=uroot;
    //scount=0;
    //sort();
    //print(database);
}
        
// Main
int main(int argc, char** argv){
   // declare a reused file variable
   FILE *file1;
      
   // read in database file
   file1=fopen("treedata.txt","r");
   if(file1==NULL)perror("file does not exist\n");
   matrix<char> basedata;
   int bdRow,bdCol;
   initCMfromFASTA(file1,basedata,bdRow,bdCol);
   convertC(basedata,bdRow,bdCol);
   fclose(file1);

   /*// randomly generated test data
   srand(time(NULL));
   int rr = rand() % 20 + 10;
   int rc = rand() % 20 + 50;
   matrix<char> testdata;
   testdata.initM(rr,rc);
 
   for(int i=0;i<rr;i++){
       for(int j=0;j<rc;j++){
           testdata(i,j)=nToAa(rand() % 20);
           cout<<testdata(i,j);
       }cout<<endl;
   }cout<<endl;*/
   
   // create a neighbor joining class instance
   joinNear JN(bdRow,bdCol);
   //joinNear test(rr,rc);

   // join neighbors, create tree, and print tree in Newick format.
   JN.joinTree(basedata);
   //test.joinTree(testdata);

   return 0; 
}
