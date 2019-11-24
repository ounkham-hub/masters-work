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
//Last Modified: October 5, 2017
//
//**********************************************************

int main()
{
//Setting up variable parameters
   int type,intsnap,nextsnap,O1,O2,sharedH;
   int intflag,nextflag;
   int i,y,z,elines,plines;
   int erow,ecoln,prow,trow,tcoln;
   int srow,slines;
   int PTeventlist[60000][8];
   int eigenlist[60000][8];
//   int maxsnap=12;
//  int tracklist[5000][maxsnap][2];
   int switchcount[5000][3];
   FILE *input,*output;


//Initialize matrices and counters
   slines=0;
   plines=0;
   elines=0;
   for(y=0; y<60000; y++){
      for(z=0; z<8; z++){
         PTeventlist[y][z]=0;
	 eigenlist[y][z]=0;
      }//end z-loop
   }//end y-loop

   for(y=0; y<5000; y++){
//      for(z=0; z<maxsnap; z++){
         switchcount[y][0]=0;
	 switchcount[y][1]=0;
	 switchcount[y][2]=0; 
//	 tracklist[y][z][0]=0;
//	 tracklist[y][z][1]=0;
//      }//end z-loop
   }//end y-loop
   
   input=fopen("labeled-PT-events-8m.txt","r");
   output=fopen("output.txt","w");

//Read in file and store entire list of PT steps into two matrices
   while(fscanf(input,"%d %d %d %d %d %d %d %d\n",&type,&intsnap,&nextsnap,&O1,&O2,&sharedH,&intflag,&nextflag)==8){ 
      PTeventlist[plines][0]=type;
      PTeventlist[plines][1]=intsnap;
      PTeventlist[plines][2]=nextsnap;
      PTeventlist[plines][3]=O1;
      PTeventlist[plines][4]=O2;
      PTeventlist[plines][5]=sharedH;
      PTeventlist[plines][6]=intflag;
      PTeventlist[plines][7]=nextflag;
      plines++;

      eigenlist[elines][0]=type;
      eigenlist[elines][1]=intsnap;
      eigenlist[elines][2]=nextsnap;
      eigenlist[elines][3]=O1;
      eigenlist[elines][4]=O2;
      eigenlist[elines][5]=sharedH;
      eigenlist[elines][6]=intflag;
      eigenlist[elines][7]=nextflag;
      elines++;
   }//end while
   fclose(input);

//Identify uniqes EZ PAIRS
   for(erow=0; erow<elines; erow++){
       for(ecoln=0; ecoln<8; ecoln++){
           if(eigenlist[erow][0]!=0){
              if(eigenlist[erow][6]==2){ //if flag is 2, or eigen
	         for(i=erow+1; i<elines; i++){ //search through the remaining list and remove duplicate eigens
                    if(eigenlist[erow][0]==eigenlist[i][0] && eigenlist[erow][3]==eigenlist[i][3] && eigenlist[erow][4]==eigenlist[i][4]){
 //matching step type, O1 and O2 - remove duplicate
		       eigenlist[i][0]=0;
		       eigenlist[i][1]=0;
		       eigenlist[i][2]=0;
		       eigenlist[i][3]=0;
		       eigenlist[i][4]=0;
		       eigenlist[i][5]=0;
		       eigenlist[i][6]=0;
		       eigenlist[i][7]=0;
		    }  
		    else if(eigenlist[erow][0]==eigenlist[i][0] && eigenlist[erow][3]==eigenlist[i][4] && eigenlist[erow][4]==eigenlist[i][3]){
//Remove if inverse pair is identified (PT occurred but it returned to it's partner)
		       eigenlist[i][0]=0;
		       eigenlist[i][1]=0;
		       eigenlist[i][2]=0;
		       eigenlist[i][3]=0;
		       eigenlist[i][4]=0;
		       eigenlist[i][5]=0;
		       eigenlist[i][6]=0;
		       eigenlist[i][7]=0;
		    }
		 }//end i-loop
              }//end flag==2
	      else if(eigenlist[erow][6]==1){//if intitial structure is Zundel - remove
		 eigenlist[erow][0]=0;
		 eigenlist[erow][1]=0;
		 eigenlist[erow][2]=0;
		 eigenlist[erow][3]=0;
		 eigenlist[erow][4]=0;
		 eigenlist[erow][5]=0;
		 eigenlist[erow][6]=0;
		 eigenlist[erow][7]=0;
	     }//end if flag==1
	  }//end if non-zero element
       }//end ecoln-loop
   }//end erow-loop
/*
   for(y=0; y<elines; y++){
      if(eigenlist[y][0]!=0){
         fprintf(output,"%d %d %d %d %d %d %d %d\n",
         eigenlist[y][0],
         eigenlist[y][1],
         eigenlist[y][2],
         eigenlist[y][3],
         eigenlist[y][4],
         eigenlist[y][5],
         eigenlist[y][6],
         eigenlist[y][7]);
      }
   }//end y-loop
*/

//Set up switchcount matrix to count the number of PTs w.r.t each O index

  for(erow=0; erow<elines; erow++){
     if(eigenlist[erow][0]!=0){
       switchcount[slines][0]=eigenlist[erow][3];//Corresponds to unique Eigen O index
       switchcount[slines][1]=eigenlist[erow][4];//Zundel partner 
       slines++;
     }//if nonzero element
  }//end e-row

//Now eigenlist needs to be compared to PTeventslist to determine the number successful PT events

   for(erow=0; erow<elines; erow++){
      if(eigenlist[erow][0]!=0){
         for(prow=0; prow<plines; prow++){
	    if(PTeventlist[prow][0]==3){ //Only focus on Z-dissociation
               if(eigenlist[erow][3]==PTeventlist[prow][3]){
	          if(eigenlist[erow][4]==PTeventlist[prow][4]){//If Z-diss to partner O, successful PT
		     for(srow=0; srow<slines; srow++){ //loop through switchcount matrix to find correct O
		        if(eigenlist[erow][3]==switchcount[srow][0] && eigenlist[erow][4]==switchcount[srow][1]){ //if O indices match
			switchcount[srow][2]++; //increment counter
		        }
		     }//end srow-loop 
	          } 
	       }//end if ZD step, 3
	    }//end if O's match
         }//end prow-loop
      }//end if eigenlist is nonzero
   }//end erow-loop 

//Print switchcount array

   for(srow=0; srow<slines; srow++){
     fprintf(output,"%d %d %d\n",switchcount[srow][0],switchcount[srow][1],switchcount[srow][2]);
   }//end srow-loop


/*
   for(erow=0; erow<elines; erow++){
      if(eigenlist[erow][0]!=0){
   	    tracklist[erow][0][0]=eigenlist[erow][3]; //copy O index
	    tracklist[erow][1][0]=switchcount[erow];// count # oscillations
	    tracklist[erow][2][0]=eigenlist[erow][1]; //initial snap
	    tracklist[erow][3][0]=eigenlist[erow][6]; //initial flag, type
	    tracklist[erow][4][0]=eigenlist[erow][2]; //next snap
	    tracklist[erow][5][0]=eigenlist[erow][7]; //final flag type
	    tracklist[erow][5][1]=eigenlist[erow][4]; //Zundel O partner
//Should correspond to EZ step
	 }
      }

   for(y=0; y<elines; y++){
      if(tracklist[y][0][0]!=0){
	printf("%d %d %d %d %d %d %d\n",
        tracklist[y][0][0],
	tracklist[y][1][0],
	tracklist[y][2][0],
	tracklist[y][2][1],
	tracklist[y][3][0],
	tracklist[y][4][0],
	tracklist[y][4][1]); //final flag, type
     }
//Should correspond to EZ step
   }//end y-loop
*/

//Compare eigenlist to PTeventlist to track step change and store into tracklist
/*
   for(trow=0; trow<elines; trow++){
      for(tcoln=5; tcoln<maxsnap; tcoln++){//continue tracking list
         if(tracklist[trow][tcoln][0]!=0){
            for(prow=0; prow<plines; prow++){
if(PTeventlist[prow][1]!=0){
	       if(tracklist[trow][0][0]==PTeventlist[prow][3]){// matching eigen indices
	          if(PTeventlist[prow][0]==12){// only Zundel intermediate step      
	          }
//Let's come back to these later
		  else if(PTeventlist[prow][0]==3){//only Zundel dissociation, 1 to 2 flag change
		     if(tracklist[trow][tcoln][1]==PTeventlist[prow][4]){ //If Zundel diss. occurs but O index differs than initial O
			printf("trigger: %d %d\n",PTeventlist[prow][4],PTeventlist[prow][1]);
			tracklist[trow][1][0]=switchcount[trow]+1; //count # of switches before actual PT event
			tracklist[trow][tcoln+1][0]=PTeventlist[prow][2];//corresponds to snap. with NEW E state,nextsnap
			tracklist[trow][tcoln+2][0]=PTeventlist[prow][7];//flag Eigen, has no partners
		     }
	          }
		  else if(PTeventlist[prow][0]==21){//only Eigen to Zundel, 2 to 1 flag change
		     if(tracklist[trow][tcoln-1][0]!=PTeventlist[prow][2]){//Unique event 
		        tracklist[trow][tcoln+1][0]=PTeventlist[prow][2];//corresponds to snap. with new Z state,nextsnap
			tracklist[trow][tcoln+1][1]=PTeventlist[prow][4];//Zundel partner
		     }
		  }
               }
}//if PTeventlist not zero
            }//end prow-loop 
	 }//if non-zero element
      }//end tcoln-loop
   }//end trow-loop

   for(y=0; y<elines; y++){
      for(z=0; z<maxsnap; z++){
         if(tracklist[y][z][0]!=0){
	    printf("%d ",tracklist[y][z][0]);
         }
      }//end z-loop
      printf("\n");
   }//end y-loop
*/
   fclose(output);

}//end main
