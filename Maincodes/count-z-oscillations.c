#include<stdio.h>
#include<stdlib.h>

//***********************************************************************
//Program created to calculate the average lifetime of Zundels
//and the # of PTs during an oscillation period. 
//
//Input: duration-of-zundel.txt (output from find-zundel-pairs-v1.c)
//Output:
//
//Created by Lelee Ounkham
//Last Modified Marach 19, 2018
//
//***********************************************************************

int main(){
   
   int Hyd,O1,O2,startsnap,finalsnap;
   int y,z,slines,plines,zlines;
   int srow,zrow,trow;
   FILE *input,*output,*hinput;
   char hfile[100]; 

   int switchmatrix[23000][5];
   int tempmatrix[23000][5];
   int zpairmatrix[23000][7];

//Initialize matrices and counters
   for(y=0; y<23000; y++){
      for(z=0; z<5; z++){
	 zpairmatrix[y][z]=0;
	 zpairmatrix[y][5]=0;
	 zpairmatrix[y][6]=0;
	 switchmatrix[y][z]=0;
	 tempmatrix[y][z]=0;
      }
   }//end y-loop
   
   zlines=0;
   slines=0;
   plines=0; 

//Read in input file
   input=fopen("duration-of-zundel.txt","r");
   output=fopen("z-oscillations.txt","w");

   while(fscanf(input,"%d %d %d %d %d\n",&Hyd,&O1,&O2,&startsnap,&finalsnap)==5){
      switchmatrix[zlines][0]=Hyd;
      switchmatrix[zlines][1]=O1;
      switchmatrix[zlines][2]=O2;
      switchmatrix[zlines][3]=startsnap;
      switchmatrix[zlines][4]=finalsnap;

      zpairmatrix[zlines][0]=Hyd;
      zpairmatrix[zlines][1]=O1;
      zpairmatrix[zlines][2]=O2;
      zpairmatrix[zlines][3]=startsnap;
      zpairmatrix[zlines][4]=finalsnap;

      tempmatrix[zlines][0]=Hyd;
      tempmatrix[zlines][1]=O1;
      tempmatrix[zlines][2]=O2;
      tempmatrix[zlines][3]=startsnap;
      tempmatrix[zlines][4]=finalsnap;
      zlines++;
   }//end while-loop
   fclose(input);

//printf("%d\n",zlines);

//***********************************************************************
//Section I: Identify all unique Zundel pairs, such that the original
//O's are the first to be identified and PTs can be identified after.
//
//switchmatrix contains all the Zundel pairs that have formed after PT
//zpairmatrix contains all the Zundel pairs prior to PT
//
//***********************************************************************

   for(srow=0; srow<zlines; srow++){
      for(zrow=0; zrow<zlines; zrow++){
         if(switchmatrix[srow][0]!=0 && zpairmatrix[zrow][0]!=0){ //for nonzero elements
	    if(switchmatrix[srow][0]==zpairmatrix[zrow][0] && switchmatrix[srow][1]==zpairmatrix[zrow][1] && switchmatrix[srow][2]==zpairmatrix[zrow][2]){ //If H's and O's match in both arrays
               for(z=0; z<zlines; z++){
	          if(switchmatrix[srow][0]==zpairmatrix[z][0]){ //If H atoms matches
		     if(switchmatrix[srow][2]==zpairmatrix[z][1] && switchmatrix[srow][1]==zpairmatrix[z][2]){ 
//Remove all O2-O1 pairs which reflect a PT
		        zpairmatrix[z][0]=0;
		        zpairmatrix[z][1]=0;
		        zpairmatrix[z][2]=0;
		        zpairmatrix[z][3]=0;
		        zpairmatrix[z][4]=0;
		     }
	          }
	       }//end z-loop
	       switchmatrix[srow][0]=0;
	       switchmatrix[srow][1]=0;
	       switchmatrix[srow][2]=0;
	       switchmatrix[srow][3]=0;
               switchmatrix[srow][4]=0;
	    }
         }
      }//end zrow-loop
   }//end row-loop

/* Cross-check
   for(z=0; z<zlines; z++){
      if(switchmatrix[z][0]!=0){
         printf("switch: %d %d %d %d %d\n",
         switchmatrix[z][0],
         switchmatrix[z][1],
         switchmatrix[z][2],
         switchmatrix[z][3],
         switchmatrix[z][4]);
	 slines++;
      }
   }

   for(z=0; z<zlines; z++){
      if(zpairmatrix[z][0]!=0){
         printf("zpair: %d %d %d %d %d\n",
         zpairmatrix[z][0],
         zpairmatrix[z][1],
         zpairmatrix[z][2],
         zpairmatrix[z][3],
         zpairmatrix[z][4]);
	 plines++;
      }
   }
   printf("slines plines: %d %d\n",slines,plines);
*/
//***********************************************************************
//Section II: Count the number of PTs during between an Zundel pair
//by counting the number of times a Zundel pair appears in switchmatrix.
//***********************************************************************

    for(zrow=0; zrow<zlines; zrow++){
       for(srow=0; srow<zlines; srow++){
	  if(zpairmatrix[zrow][0]!=0 && switchmatrix[srow][0]!=0){
	     if(zpairmatrix[zrow][0]==switchmatrix[srow][0]){ //If shared H's matches
	        if(zpairmatrix[zrow][2]==switchmatrix[srow][1] && zpairmatrix[zrow][1]==switchmatrix[srow][2]){
//If O2 in zpairmatrix equals O1 in switchmatrix (i.e. O2-O38 to O38-O2)
	           for(z=0; z<zlines; z++){ //loop through to see how many times this pair appears
                      if(switchmatrix[z][0]!=0){
		         if(switchmatrix[srow][0]==switchmatrix[z][0] && switchmatrix[srow][1]==switchmatrix[z][1] && switchmatrix[srow][2]==switchmatrix[z][2]){ //if identical pair is found
		            zpairmatrix[zrow][5]=zpairmatrix[zrow][5]+1;
/*
if(switchmatrix[srow][0]==181){
printf("%d %d %d %d %d\n",
switchmatrix[z][0],
switchmatrix[z][1],
switchmatrix[z][2],
switchmatrix[z][3],
switchmatrix[z][4]);
}
*/		            switchmatrix[z][0]=0;
		         }
                       }
		    }//end z-loop
		 }
	      }
	   }
        }//end srow-loop 
     }//end zrow-loop

//***********************************************************************
//Section III: Identify the number of unsuccessful PTs by counting
//how many times the original Zundel reforms. 
//***********************************************************************
     for(trow=0; trow<zlines; trow++){
        for(zrow=0; zrow<zlines; zrow++){
	   if(zpairmatrix[zrow][0]!=0 && tempmatrix[trow][0]!=0){
	      if(tempmatrix[trow][0]==zpairmatrix[zrow][0] && tempmatrix[trow][1]==zpairmatrix[zrow][1] && tempmatrix[trow][2]==zpairmatrix[zrow][2]){ 
//If H index, orig. and partner O's matches
		 for(z=0; z<zlines; z++){
		   if(tempmatrix[trow][0]==tempmatrix[z][0] && tempmatrix[trow][1]==tempmatrix[z][1] && tempmatrix[trow][2]==tempmatrix[z][2]){
		      zpairmatrix[zrow][6]=zpairmatrix[zrow][6]+1; //Increase count
		      tempmatrix[z][0]=0; //zero out duplicate nonsuccessful events
		   } 
		 } //end z-loop
	      }
	   }
	}//end nrow-loop	
     }//end zrow-loop 

      for(z=0; z<zlines; z++){
         if(zpairmatrix[z][0]!=0 && zpairmatrix[z][6]!=0){
            fprintf(output,"%d %d %d %d %d %d\n",
            zpairmatrix[z][0], //Shared H
            zpairmatrix[z][1], //Orig. O index
            zpairmatrix[z][2], //Partner O index
            zpairmatrix[z][3], //Time of Zundel formation
	    zpairmatrix[z][5], //# PTs
	    zpairmatrix[z][6]); //# unsuccessful PTs
         }
      }
   fclose(output);
}//end main-loop
