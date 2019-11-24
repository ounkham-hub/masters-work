#include<stdio.h>
#include<stdlib.h>

//*******************************************************************************
//Program created to identify how many replicas participated in a PT event
//given a +/- 5 time snapshot interval. 
//
//
//Input: #-total-PT-events.txt
//Output:
//*The file 0-total-PT-events.txt corresponds to the centroid data 
//
//Created by Lelee Ounkham
//Last modified: March 25, 2018
//*******************************************************************************

int main(){

  int bead,Hindex,currentO,partnerO,event; 
  int y,z,crow,row,clines,rlines;
  FILE *input,*output;
  char ifile[100];
 
  int centroidmatrix[20000][8];
  int replicamatrix[60000][4];

//Initialize matrices
   for(y=0; y<10000; y++){
      for(z=0; z<8; z++){
         centroidmatrix[y][z]=0;
      }//end z-loop
   }//end y-loop

   for(y=0; y<60000; y++){
      for(z=0; z<4; z++){
	 replicamatrix[y][z]=0;
      }
   } //end y-loop
   clines=0;
   rlines=0;
  
//Change variables accordingly
  output=fopen("all-Successful-PTs-20snap.txt","w");
  for(bead=0; bead<=32; bead++){
     sprintf(ifile,"%d-all-H-PT-events.txt",bead);
     input=fopen(ifile,"r");

      while(fscanf(input,"%d %d %d %d\n",&Hindex,&currentO,&partnerO,&event)==4){
         if(bead==0){
	    centroidmatrix[clines][0]=Hindex;
	    centroidmatrix[clines][1]=currentO;
	    centroidmatrix[clines][2]=partnerO;
	    centroidmatrix[clines][3]=event;
	    centroidmatrix[clines][4]=event-20; //minimum snap. range
	    centroidmatrix[clines][5]=event+20; //maximum snap. range
	    centroidmatrix[clines][6]=0; //counter
	    clines++;
         }//end if bead==0
         else if(bead!=0){
	    replicamatrix[rlines][0]=Hindex;
	    replicamatrix[rlines][1]=currentO;
	    replicamatrix[rlines][2]=partnerO;
	    replicamatrix[rlines][3]=event;
	    rlines++;
	 }
      }//end while-loop
      fclose(input);

//*******************************************************************************
//Section to compaire centroidmatrix to replicamatrix to determine which
//beads participated in PT events.
//
//*******************************************************************************
      for(crow=0; crow<clines; crow++){
         for(row=0; row<rlines; row++){
	    if(replicamatrix[row][0]!=0){
	       if(centroidmatrix[crow][0]==replicamatrix[row][0] && centroidmatrix[crow][1]==replicamatrix[row][1] && centroidmatrix[crow][2]==replicamatrix[row][2]){ //If H index and both O indices matches in centroid and replica matrix
	          if(replicamatrix[row][3] >= centroidmatrix[crow][4] && replicamatrix[row][3] <= centroidmatrix[crow][5]){
//If event in replica is within +/- snapshots of the PT in centroid
		     if(centroidmatrix[crow][7]==0){ //replica not counted for particular event
		        centroidmatrix[crow][6]=centroidmatrix[crow][6]+1; //update count
			centroidmatrix[crow][7]=1; //flag event

/*if(replicamatrix[row][0]==304 && replicamatrix[row][1]==56 && replicamatrix[row][2]==40){
printf("%d %d %d %d %d\n",
bead,
replicamatrix[row][0],
replicamatrix[row][1],
replicamatrix[row][2],
replicamatrix[row][3]);
}
*/		        replicamatrix[row][0]=0; //remove line from replicamatrix

		     }
		  }
	       }
	    }
	 }//end row-loop
      }//end crow-loop

//Initialize replicamatrix for next bead dataset to be read in
      for(y=0; y<rlines; y++){
         for(z=0; z<4; z++){
	    replicamatrix[y][z]=0;
         }
      } //end y-loop

      for(y=0; y<clines; y++){
         centroidmatrix[y][7]=0; //0 out flags for next set of beads
      }
      rlines=0;
   }//end bead-loop

//Print out PT events with associated # of particpating beads

   for(crow=0; crow<clines; crow++){
if(centroidmatrix[crow][6]==32){
      fprintf(output,"%d %d %d %d %d\n",
      centroidmatrix[crow][0], //Hindex
      centroidmatrix[crow][1], //current O index
      centroidmatrix[crow][2], //partner O index
      centroidmatrix[crow][3], //PT snapshot
      centroidmatrix[crow][6]);	//# of participating beads
}
   }//end crow-loop
   fclose(output);

}//end main-loop
