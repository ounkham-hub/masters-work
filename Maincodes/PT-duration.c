#include<stdio.h>
#include<stdlib.h>

//******************************************************************************
//Program created to track a specific proton and two of neighboring
//O atoms by flagging each of the structures as Zundels,Eigens or Waters.
//Then calculate the duration to which the proton remains on the new 
//O atom or is "successful" transferred before returning to it's initial
//O or transforming into another structure.
//
//Slightly different for centoid data than PIMD
//*To use this program, you need to know the time frame that the proton
//and it's O's undergo structural oscillations. i.e. H106 undergoes
//structural changes with O2 and O38 from snaps. 1 to 16600. 
//
//Note the flags indicate the following:
//    - 0, water
//    - 1, eigen
//    - 2, zundel
//
//Input: #-event-list.txt
//Output: 
//
//Created by Lelee Ounkham
//Last modified: February 9, 2018
//
//****************************************************************************

int main(){

//Parameters to read in input file
   int type,startsnap,nextsnap,Oindex1;
   int startH1,startH2,startH3,startZH;
   int nextH1,nextH2,nextH3,nextZH,Oindex2;
   int y,z,row,row2,lines,bead,proton;
   int initialtimeframe,finaltimeframe;
   int targetO,snapdiff;
   int targetOlist[100000][13];
   FILE *input,*output;
   char ifile[100],ofile[100];

//CHANGE ACCORDINGLY
   proton=106; 
   targetO=5; //The O atom of new water received the proton during PT
   initialtimeframe=41230; //timeframe starts at initial to final
   finaltimeframe=61555;

   for(bead=1; bead<=1; bead++){
      sprintf(ifile,"%d-event-list.txt",bead);
      input=fopen(ifile,"r");
      output=fopen(ofile,"w");
  
//Initialize matrices and counters 
      for(y=0; y<100000; y++){
         for(z=0; z<14; z++){
            targetOlist[y][z]=0;
         }//end z-loop
      }//end y-loop
      lines=0;

//Store input file into an independent matrix, targetOlist
      while(fscanf(input,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
      &type,&startsnap,&nextsnap,&Oindex1,&startH1,&startH2,&startH3,
      &Oindex2,&nextH1,&nextH2,&nextH3,&startZH,&nextZH)==13){
         if(startsnap >= initialtimeframe && startsnap <= finaltimeframe){ //Only store information if it's within specified windowtime
            if(Oindex1==targetO){ //if the O from input is equals targetO (O of interest)
               targetOlist[lines][0]=type;
               targetOlist[lines][1]=startsnap;
               targetOlist[lines][2]=nextsnap;
               targetOlist[lines][3]=Oindex1;
               targetOlist[lines][4]=startH1;
               targetOlist[lines][5]=startH2;
               targetOlist[lines][6]=startH3;
	       targetOlist[lines][7]=Oindex2;
               targetOlist[lines][8]=nextH1;
               targetOlist[lines][9]=nextH2;
               targetOlist[lines][10]=nextH3;;
               targetOlist[lines][11]=startZH;
               targetOlist[lines][12]=nextZH;
               lines++;
	    }
	 }
      }//end while-loop
      fclose(input);

//*********************************************************************
//Section to identify successful PTs (when Z dissociates in E, ZE)
//and calculate the duration that the proton is transferred.
//*********************************************************************
     for(row=0; row<lines; row++){
        if(targetOlist[row][0]!=55){
           if(targetOlist[row][0]==21){ //if ZE transfer, indicative of PT 
	      if(targetOlist[row][0]!=targetOlist[row+1][0]){
		 snapdiff=(targetOlist[row+1][1]+1)-targetOlist[row][2]; //Accounts for snapdiffs that should be 1 but are 0
		
		 if(snapdiff>0){ //if snapdiff is non-negative
	             printf("%d %d %d %d %d %d %d\n",
		     bead,
		     targetOlist[row][3], //target O index
		     targetOlist[row][2], //snap when proton is on target O (nextsnap)
		     targetOlist[row+1][1], //final snap before proton goes to another state
		     snapdiff, //diff. in time
		     targetOlist[row][0], //initial type (should be Eigen or 1)
		     targetOlist[row+1][0]); //final type (whatever Eigen tranformed into)
		  }
		     targetOlist[row][0]=55; //marked as analyzed
	      } 
	   } 
        }
     }//end row-loop
   }//end bead-loop
}//end main-loop
