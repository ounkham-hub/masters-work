#include<stdio.h>
#include<stdlib.h>

//***************************************************************************
//Program created to count the number of O1-H* and O2-H* weights among
//Zundels and since these Zundels are not undergoing reactivity symmetric
//weights need to be accounted for. In other words a O1-H* and O2-H* of 
//1.000 and 0.031 is the same as an  O1-H* and O2-H* weight of 0.031 and
//1.000.
//
//Input: total-WN-Zundel-shared-weights-nocutoffs.txt
//	 potential-weights.txt
//Output:
//
//Created  by: Lelee Ounkham
//Date last modified: November 2, 2018
//
//***************************************************************************

int main(){
   int O1,O2,H1,H2,lines,mlines,y,z,row,mrow,snap;
   float wO1H1,wO2H2,weight1,weight2;
   FILE *input,*winput;

   int indexmatrix[536000][3];
   int maincount[716];
   float weightmatrix[536000][3];
   float mainmatrix[716][3];

   for(y=0; y<536000; y++){
      for(z=0; z<3; z++){
	 indexmatrix[y][z]=0;
         weightmatrix[y][z]=0.000;
      }
   }
   for(y=0; y<716; y++){
      maincount[y]=0;
      mainmatrix[y][0]=0.000;
      mainmatrix[y][1]=0.000;
   } 
   lines=0;
   mlines=0;

//Read in file and store into matrix 

   input=fopen("total-WN-Zundel-shared-weights-nocutoffs.txt","r");
   winput=fopen("potential-weights.txt","r");

   while(fscanf(input,"%d %d %d %d %d %f %f\n",&snap,&O1,&O2,&H1,&H2,&wO1H1,&wO2H2)==7){
      indexmatrix[lines][0]=O1;
      indexmatrix[lines][1]=O2;
      indexmatrix[lines][2]=H1; //should be the same as H2
      weightmatrix[lines][1]=wO1H1;
      weightmatrix[lines][2]=wO2H2;
      lines++;

   }//end while-loop
   fclose(input);

   while(fscanf(winput,"%f %f\n",&weight1,&weight2)==2){
      mainmatrix[mlines][1]=weight1;
      mainmatrix[mlines][2]=weight2;
      mlines++;

//      printf("%0.3f %0.3f\n",weight1,weight2);
   }//end while-loop
   fclose(winput);

//printf("%d %d\n",lines,mlines);
//***************************************************************************
//Section to identify and count weight pairings. O and H indices can be 
//outputted as a crosscheck.
//***************************************************************************

   for(mrow=0; mrow<mlines; mrow++){
      for(row=0; row<lines; row++){
	 if(weightmatrix[row][1]!=0.000 && weightmatrix[row][2]!=0.000){
	    if(mainmatrix[mrow][1]==weightmatrix[row][1] && mainmatrix[mrow][2]==weightmatrix[row][2]){
//If weights matches in same order as listed
	       maincount[mrow]=maincount[mrow]+1; //increment count;
	       
	       weightmatrix[row][1]=0.000;
	       weightmatrix[row][2]=0.000;


	       for(y=0; y<lines; y++){
		  if(weightmatrix[y][1]!=0.000 && weightmatrix[y][2]!=0.000){
	             if(mainmatrix[mrow][1]==weightmatrix[y][1] && mainmatrix[mrow][2]==weightmatrix[y][2]){
		    	maincount[mrow]=maincount[mrow]+1; 
		        weightmatrix[y][1]=0.000;
		        weightmatrix[y][2]=0.000;

		     }
		     else if(mainmatrix[mrow][1]==weightmatrix[y][2] && mainmatrix[mrow][2]==weightmatrix[y][1]){
			maincount[mrow]=maincount[mrow]+1; 
		        weightmatrix[y][1]=0.000;
		        weightmatrix[y][2]=0.000;
		     }
	          }//end y-loop 
	       }//end y-loop
//printf("%d %0.3f %0.3f %d\n",mrow,mainmatrix[mrow][1],mainmatrix[mrow][2],maincount[mrow]);	
	    }
	 }
      }//end row-loop
   }//end mrow-loop

//***************************************************************************
//Section to print maincount matrix with all the counts and weights
//***************************************************************************

   for(y=0; y<mlines; y++){
     if(maincount[y]!=0){
       printf("%0.3f %0.3f %d\n",mainmatrix[y][1],mainmatrix[y][2],maincount[y]);
     }
   }

}//end main-loop
   
