#include<stdio.h>
#include<stdlib.h>

//**********************************************************
//Program written to identify and quantify the number of
//successful PT events.
//
//Note:
//Flags are either 2, Eigen, or 1, Zundel
//The 3 types of PT steps
// - 12 (Zundel to Eigen or ZZ intermediate)
// - 21 (Eigen to Zundel)
// - 3 (Zundel dissociation)
//
//Written by: Lelee Ounkham
//Last Modified: January 8, 2018
//
//**********************************************************

int main()
{
//Setting up variable parameters
   int type,intsnap,nextsnap,O1,O2,sharedH;
   int i,y,z,elines,plines;
   int erow,ecoln,prow,trow,tcoln;
   int PTeventlist[150000][8];
   int eigenlist[150000][8];
   FILE *input,*output,*output2;


//Initialize matrices and counters
   plines=0;
   elines=0;
   for(y=0; y<150000; y++){
      for(z=0; z<8; z++){
         PTeventlist[y][z]=0;
	 eigenlist[y][z]=0;
	 pweightlist[y][0]=0.0;
	 pweightlist[y][1]=0.0;
	 pweightlist[y][2]=0.0;
	 eweightlist[y][0]=0.0;
	 eweightlist[y][1]=0.0;
	 eweightlist[y][2]=0.0;
      }//end z-loop
   }//end y-loop

   input=fopen("total-seq-PT-steps.txt","r");
   output=fopen("total-succ-wPT-events.txt","w");
   output2=fopen("remainders.txt","w");
//Read in file and store entire list of PT steps into two matrices
   while(fscanf(input,"%d %d %d %d %d %d\n",&type,&intsnap,&nextsnap,&O1,&O2,&sharedH)==6){ 
      PTeventlist[plines][0]=type;
      PTeventlist[plines][1]=intsnap;
      PTeventlist[plines][2]=nextsnap;
      PTeventlist[plines][3]=O1;
      PTeventlist[plines][4]=O2;
      PTeventlist[plines][5]=sharedH;
      plines++;

      if(type==21){ //If E to Z transfer store into eigenlist
         eigenlist[elines][0]=type;
         eigenlist[elines][1]=intsnap;
         eigenlist[elines][2]=nextsnap;
         eigenlist[elines][3]=O1;
         eigenlist[elines][4]=O2;
         eigenlist[elines][5]=sharedH;
         elines++;
      }
   }//end while
   fclose(input);


   for(erow=0; erow<elines; erow++){
      for(prow=0; prow<plines; prow++){
        if(PTeventlist[prow][0]==3){ //only for Zundel diss.
	   if(eigenlist[erow][4]==PTeventlist[prow][4]){
	      fprintf(output,"%d %d %d %d %d %d %d\n",
		 eigenlist[erow][1],//initial E snap.
		 eigenlist[erow][2],//initial Z snap.
		 PTeventlist[prow][1],//final Z snap.
		 PTeventlist[prow][2],//final E snap.
		 eigenlist[erow][3], //initial E O index
		 eigenlist[erow][4], //intial Zundel O PARTNER - assume orig. E is the other O index
		 PTeventlist[prow][5]); //shared proton

		 PTeventlist[prow][0]=0; //Zero out type, so step is no recounted
		 eigenlist[erow][0]=0;
		 eigenlist[erow][1]=0;
		 eigenlist[erow][2]=0;
		 eigenlist[erow][3]=0;
		 eigenlist[erow][4]=0;
		 eigenlist[erow][5]=0;
	   }
	}
      }//end erow-loop
   }//end prow-loop

   for(z=0; z<elines; z++){
      if(eigenlist[z][0]!=0){
         fprintf(output2,"%d %d %d %d %d %d\n",
         eigenlist[z][0],
         eigenlist[z][1],
         eigenlist[z][2],
         eigenlist[z][3],
         eigenlist[z][4],
         eigenlist[z][5]);
       }//end z-loop
    }
    fclose(output);
    fclose(output2);  

}//end main
