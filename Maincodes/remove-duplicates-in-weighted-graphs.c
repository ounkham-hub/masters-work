#include<stdio.h>
#include<stdlib.h>

//**********************************************************************
//Program created to remove duplicate lines in weighted graphs.
//
//Input: weighted-O#.H#.Graph 
//Output:
//
//Created by Lelee Ounkham
//Last modified August 31, 2018
//**********************************************************************

int main(){

   int index1,type1,rindex1,index2,type2,rindex2;
   int snap,lines,y,z,row,nrow,bead;
   float weight;
   FILE *input,*output;
   char ifile[100],ofile[100];

   int indexmatrix[1000][6];
   float weightmatrix[1000][1];
   
   lines=0;
//   output=fopen("total-centroid-weights-n20-v2.txt","a");
   for(snap=1; snap<=2; snap++){
      for(bead=1; bead<=32; bead++){
         sprintf(ifile,"%d.O%d.H%d.wGraph",bead,snap,snap);
         input=fopen(ifile,"r");
         sprintf(ofile,"%d.O%d.H%d.wGraphGeod",bead,snap,snap);
         output=fopen(ofile,"w");

      while(fscanf(input,"%d %d %d %d %d %d %f\n",&index1,&type1,&rindex1,&index2,&type2,&rindex2,&weight)==7){
        indexmatrix[lines][0]=index1;
	indexmatrix[lines][1]=type1;
	indexmatrix[lines][2]=rindex1;
	indexmatrix[lines][3]=index2;
	indexmatrix[lines][4]=type2;
	indexmatrix[lines][5]=rindex2;
	weightmatrix[lines][0]=weight; 
	lines++;
      }//end while-loop
      fclose(input);

      for(row=0; row<lines; row++){
         for(nrow=1; nrow<lines; nrow++){
	    if(indexmatrix[row][0]!=0 && indexmatrix[nrow][0]!=0){
	       if(indexmatrix[row][1]==indexmatrix[nrow][4] && indexmatrix[row][0]==indexmatrix[nrow][3]){
//if flag1 matches flag2 and index1 matches index2
	          if(indexmatrix[row][2]==indexmatrix[nrow][5] && weightmatrix[row][0]==weightmatrix[nrow][0]){
//if rindex1 matches rindex2 and weight
		     fprintf(output,"%d %d %d %0.3f\n",
		     snap,
		     indexmatrix[row][0],
//		     indexmatrix[row][1],
//		     indexmatrix[row][2],
		     indexmatrix[row][3],
//		     indexmatrix[row][4],
//		     indexmatrix[row][5],
		     weightmatrix[row][0]);
/*
		     fprintf(output,"%d %d %d %d %d %d %f\n",
		     indexmatrix[nrow][0],
		     indexmatrix[nrow][1],
		     indexmatrix[nrow][2],
		     indexmatrix[nrow][3],
		     indexmatrix[nrow][4],
		     indexmatrix[nrow][5],
		     weightmatrix[nrow][0]);
*/
		     indexmatrix[nrow][0]=0;
		     

		  }
	       }
	    }
	 }//end nrow-loop
      }//end row-loop

//Initialize matrix for next set of snapshots
      for(y=0; y<lines; y++){
	 for(z=0; z<6; z++){
	    indexmatrix[y][z]=0;
	    weightmatrix[y][z]=0.0;
	 }
      }
      lines=0;
      fclose(output);
      }//end bead-loop
   }//end snap-loop
}///end main-loop
