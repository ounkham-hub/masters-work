#include<stdio.h>
#include<stdlib.h>

//***************************************************************************
//Program created to calculate the # of snapshots after a PT and when 
//a new reaction is formed. This will yield information on how long
//often a proton is reacts. 
//
//Created by Lelee Ounkham
//Last modified on August 3,2018
//
//***************************************************************************


int main(){

  int Hindex,currentO,partnerO,snapshot,lines; 
  int row,crow,diff,x,y;
  int eventmatrix[13000][4];
  int comparematrix[13000][4];
  FILE *input; 


//Initialize matrix and counter
   for(x=0; x<13000; x++){
      for(y=0; y<4; y++){
	 eventmatrix[x][y]=0;
	 comparematrix[x][y]=0;
      }
   }
   
   lines=0;
   input=fopen("list-of-return-after-rxn.txt","r");
   while(fscanf(input,"%d %d %d %d\n",&Hindex,&currentO,&partnerO,&snapshot)==4){
      eventmatrix[lines][0]=Hindex;
      eventmatrix[lines][1]=currentO;
      eventmatrix[lines][2]=partnerO;
      eventmatrix[lines][3]=snapshot;

      comparematrix[lines][0]=Hindex;
      comparematrix[lines][1]=currentO;
      comparematrix[lines][2]=partnerO;
      comparematrix[lines][3]=snapshot;
      lines++;
   }//end while-loop
   fclose(input);
 
    for(row=0; row<lines; row++){
//printf("%d %d %d %d\n",eventmatrix[row][0],eventmatrix[row][1],eventmatrix[row][2],eventmatrix[row][3]);
     for(crow=1; crow<lines; crow++){
	if(eventmatrix[row][0]!=0 && comparematrix[crow][0]!=0){
	   if(eventmatrix[row][3]<comparematrix[crow][3]){ //If snapshot # don't match

	      if(eventmatrix[row][0]==comparematrix[crow][0] && eventmatrix[row][1]==comparematrix[crow][1] && eventmatrix[row][2]==comparematrix[crow][2]){
//If current and partner O index, and H matches
	          diff=comparematrix[crow][3]-eventmatrix[row][3];
/*
	          printf("%d %d %d %d %d %d\n",
	          eventmatrix[row][0],
                  eventmatrix[row][1],
                  eventmatrix[row][2],
                  comparematrix[crow][3],
                  eventmatrix[row][3],
	          diff);
*/

	          eventmatrix[row][0]=0;
	          eventmatrix[row][1]=0;
	          eventmatrix[row][2]=0;
	          eventmatrix[row][3]=0; 

	          comparematrix[crow][0]=0;
	          comparematrix[crow][1]=0;
	          comparematrix[crow][2]=0;
	          comparematrix[crow][3]=0; 
	       }
	    }     	   
	 }
      }//end crow-loop 
   }//end row-loop

//   for(row=0; row<lines; row++){
//      if(eventmatrix[row][0]!=0){
//         printf("%d %d %d %d\n",eventmatrix[row][0],eventmatrix[row][1],eventmatrix[row][2],eventmatrix[row][3]);
//      }
//   }//end row-loop
}//end main-loop
