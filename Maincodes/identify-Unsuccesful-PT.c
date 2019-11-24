#include<stdio.h>
#include<stdlib.h>

//******************************************************************************
//**Initial program to use is track-specific-events-v2.exe to obtain 
//centroid data.
//
//Program created to track a specific proton in time output all
//associated structure changes especially 21 or Zundel to Eigen 
//flag types.
//
//An additional output will contain information on 
//successful PTs along with corresponding replica #'s. This will allow us to 
//determine which replicas have or have not undergone PT.
//
//To use this code, you need to know the H index and the two O atoms that the 
//proton is oscillating between. Time frames selected based on previous analysis.
//In the case of H106, O2 and O38 the time frame selected was 150 snapshots 
//
//**To identify an unsuccesful PT, we need to read in the input from the centroid
//data and see if the same event is observed in the replica data**
//
//Note the flags indicate the following:
//    - 0, water
//    - 1, eigen
//    - 2, zundel
//
//Input: H#-PT-events.txt
//Output: Count of # replicas that underwent PT in specified time window
//
//Created by Lelee Ounkham
//Last modified: February 7, 2018
//
//****************************************************************************

int main(){

//Parameters to read in input file
   int type,startsnap,nextsnap,Oindex1;
   int startH1,startH2,startH3,startZH;
   int nextH1,nextH2,nextH3,nextZH,Oindex2;
   int y,z,crow,coln,clines;
   int rrow,rcoln,rlines,newcount;
   int beadcount,bead,proton,EigenOIndex;
   int centroidlist[30000][17];
   int replicalist[30000][13];
   FILE *input,*beadinput,*output;
   char beadfile[100],ofile[100];


//Input proton of interest
      proton=106;
      EigenOIndex=2;

//Initialize matrices and counters 
      for(y=0; y<30000; y++){
         for(z=0; z<18; z++){
            centroidlist[y][z]=0;
         }//end z-loop
      }//end y-loop
      clines=0;

      input=fopen("centroid-UnsuccesfulPT-ZE-H106-O2-1-16600.txt","r");
//Store input file into an independent matrix,centroidlist
      while(fscanf(input,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
      &type,&startsnap,&nextsnap,&Oindex1,&startH1,&startH2,&startH3,
      &Oindex2,&nextH1,&nextH2,&nextH3,&startZH,&nextZH)==13){

         centroidlist[clines][0]=type;
         centroidlist[clines][1]=startsnap;
         centroidlist[clines][2]=nextsnap;
         centroidlist[clines][3]=Oindex1;
         centroidlist[clines][4]=startH1;
         centroidlist[clines][5]=startH2;
         centroidlist[clines][6]=startH3;
	 centroidlist[clines][7]=Oindex2;
         centroidlist[clines][8]=nextH1;
         centroidlist[clines][9]=nextH2;
         centroidlist[clines][10]=nextH3;;
         centroidlist[clines][11]=startZH;
         centroidlist[clines][12]=nextZH;

//************************************************************************
//Set up timewindow of interest. For example, if your PT of interest is
//380 then any replica that has a type 21 transfer within +-10 snaps.
//is counted as successful PT.
//
//************************************************************************
         centroidlist[clines][13]=startsnap-5; //minimum snap. time
	 centroidlist[clines][14]=startsnap+5; //maximum snap. time
         clines++;
      }//end while-loop
      fclose(input);

//for(z=0; z<clines; z++){
//printf("%d %d\n",centroidlist[z][13],centroidlist[z][14]);
//}

      beadcount=0;
      output=fopen("rep-Unsuccessful-PTs-H106-O2-1-16600.txt","w");
      for(bead=1; bead<=32; bead++){
         sprintf(beadfile,"H%d-rep%d-PTevents.txt",proton,bead);
         beadinput=fopen(beadfile,"r");
  
//Initialize matrices and counters 
         for(y=0; y<30000; y++){
            for(z=0; z<14; z++){
	       replicalist[y][z]=0;
            }//end z-loop
         }//end y-loop
         rlines=0;

//Store input file into an independent matrix,origlist
        while(fscanf(beadinput,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
         &type,&startsnap,&nextsnap,&Oindex1,&startH1,&startH2,&startH3,
         &Oindex2,&nextH1,&nextH2,&nextH3,&startZH,&nextZH)==13){
            if(type==21){ //if Zundel dissoctation
	       if(Oindex1==EigenOIndex){
                  replicalist[rlines][0]=type;
                  replicalist[rlines][1]=startsnap;
                  replicalist[rlines][2]=nextsnap;
                  replicalist[rlines][3]=Oindex1;
                  replicalist[rlines][4]=startH1;
                  replicalist[rlines][5]=startH2;
                  replicalist[rlines][6]=startH3;
	          replicalist[rlines][7]=Oindex2;
                  replicalist[rlines][8]=nextH1;
                  replicalist[rlines][9]=nextH2;
                  replicalist[rlines][10]=nextH3;;
                  replicalist[rlines][11]=startZH;
                  replicalist[rlines][12]=nextZH;
                  rlines++;
	       }
	    }
         }//end while-loop
         fclose(beadinput);

         for(crow=0; crow<clines; crow++){
	    for(coln=2; coln<13; coln++){
	       for(rrow=0; rrow<rlines; rrow++){
	          for(rcoln=2; rcoln<13; rcoln++){
		     if(replicalist[rrow][1] >= centroidlist[crow][13] && replicalist[rrow][1] <= centroidlist[crow][14]){
//If the PT in replica data falls with specified min. and max time window
//			if(replicalist[rrow][0]!=0){
	 	           if(replicalist[rrow][rcoln]==centroidlist[crow][coln]){ //if all the elements in centroid and replicalist match
/*			      printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			      bead,
			      replicalist[rrow][0],
			      replicalist[rrow][1],
			      replicalist[rrow][2],
			      replicalist[rrow][3],
			      replicalist[rrow][4],
			      replicalist[rrow][5],
			      replicalist[rrow][6],
			      replicalist[rrow][7],
			      replicalist[rrow][8],
			      replicalist[rrow][9],
			      replicalist[rrow][10],
			      replicalist[rrow][11],
			      replicalist[rrow][12]);
*/
//			      replicalist[rrow][0]=0;

			      if(centroidlist[crow][16]==0){ //if flag is turned off
			         if(bead==1){			
				    centroidlist[crow][16]=1; //One replica identified, flag on
				    centroidlist[crow][15]=bead; //Keeps count of how many replicas are founnd for PT in centroid
			         }	
			         else if(centroidlist[crow][15]!=bead-1){//Occurs when some replicas do not undergo an event
				    centroidlist[crow][16]=1;
				    newcount=centroidlist[crow][15]+1;
				    centroidlist[crow][15]=newcount; //add 1 to whatever current count is 
			         }
			         else if(centroidlist[crow][15]==bead-1){ //If beadcount is consistent, add 1
			            centroidlist[crow][16]=1;
				    centroidlist[crow][15]=bead;  
		 	         }	
			      }
			  }
//		       }			 
		    }
	         }//end rcoln-loop
	      }//end rrow-loop
	   } //end ccoln-loop
        } //end crow-loop

        for(z=0; z<clines; z++){
        centroidlist[z][16]=0;
       }

     }//end bead-loop
              
     for(crow=0; crow<clines; crow++){
        fprintf(output,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
        centroidlist[crow][0],
        centroidlist[crow][1],
        centroidlist[crow][2],
        centroidlist[crow][3],
        centroidlist[crow][4],
        centroidlist[crow][5],
        centroidlist[crow][6],
	centroidlist[crow][7],
        centroidlist[crow][8],
        centroidlist[crow][9],
        centroidlist[crow][10],
        centroidlist[crow][11],
        centroidlist[crow][12],
	centroidlist[crow][15]); //bead count
     }//end crow-loop
     fclose(output);
}//end main-loop

