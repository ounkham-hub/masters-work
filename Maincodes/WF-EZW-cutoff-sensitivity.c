#include<stdio.h>
#include<stdlib.h>

//*************************************************************************
//Program created to calculate the coordination number for an O atom and 
//OH weight by employing a user-input cutoff.
//
//Input: R-avg-O#.H#.wGraphGeod
//Output: 
//      coln 1: O index
//      coln 2: Coordination #
//      coln 3: avg. weight
//      (optional) colns 4+: H indices and OH weights
//     
//
//Created by Lelee Ounkham
//Last modified September 13, 2018
//
//*************************************************************************

int main(){
   int snap,Oxy,Hyd,rep,sindex,Oindex,Hhold;
   int y,z,maxsnap,countOindex,clines,colns;
   int bcount,blines;
   float weight,whold;
   FILE *input,*output;
   char ifile[100];

   int CBmatrix[1000][20];
   float weightmatrix[1000][20];

   output=fopen("nocutoff-CN.txt","a");
   for(snap=1; snap<=61555; snap++){
//Intialize matrices and counters
      for(y=0; y<1000; y++){
	 for(z=0; z<20; z++){
	    CBmatrix[y][z]=0;
	    weightmatrix[y][z]=0.0;
	 }
      }
      clines=0;
      colns=2; //starts at to because of position in matrix
      countOindex=1;

      sprintf(ifile,"R-avg-O%d.H%d.wGraphGeod",snap,snap); //open current snap
      input=fopen(ifile,"r");
//Read in input
      while(fscanf(input,"%d %d %d %f\n",&sindex,&Oxy,&Hyd,&weight)==4){
	 if(countOindex==1){ //New line, store info on O atom
	     Oindex=Oxy;
	     Hhold=Hyd;
	     whold=weight;
         }
	 else if(countOindex>1){
	    if(Oindex==Oxy){ //if O index matches the previous one, store info
	       if(Hhold!=Hyd){
	          if(countOindex==2){
	             CBmatrix[clines][0]=snap;
	             CBmatrix[clines][1]=Oindex;
	             CBmatrix[clines][colns]=Hhold;
	             weightmatrix[clines][colns]=whold;
		     colns++;
		   }
		   CBmatrix[clines][0]=snap;
		   CBmatrix[clines][1]=Oxy;
		   CBmatrix[clines][colns]=Hyd;
		   weightmatrix[clines][colns]=weight;
		   Oindex=Oxy;
		   Hhold=Hyd;
		   whold=weight;
		   colns++;
	       }	
	       else if(Hhold==Hyd){ 
//For instances when the duplicates are the first two bonds of O, only count 1
	          CBmatrix[clines][0]=snap;
	          CBmatrix[clines][1]=Oindex;
	          CBmatrix[clines][colns]=Hhold;
	          weightmatrix[clines][colns]=whold;
		  colns++;
	       }
	    }
	    else if(Oindex!=Oxy){ //If O's don't match, then new O appeared - restart loop
//printf("%d %d\n",Oxy,Hyd);
	       countOindex=1;
	       colns=2;
	       clines++;

  	       Oindex=Oxy;
	       Hhold=Hyd;
	       whold=weight;
	    }
	 }
	 countOindex++;
      }//end while-loop
      fclose(input);
     
      bcount=0;
      for(y=0; y<clines+1; y++){
	 for(z=0; z<20; z++){
	    if(CBmatrix[y][z]!=0){
//	       printf("%d ",CBmatrix[y][z]);
	       bcount++;
	    }
	 }
fprintf(output,"%d %d %d\n",snap,CBmatrix[y][1],bcount-2);
	 blines=bcount;
//	 printf("\n");
	 for(z=2; z<blines; z++){
//	    printf("%0.3f ",weightmatrix[y][z]);
	}
	bcount=0;
//	printf("\n");
      }

   }//end snap-loop
   fclose(output);
}//end main-loop
