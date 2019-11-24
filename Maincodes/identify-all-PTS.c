#include<stdio.h>
#include<stdlib.h>

//**************************************************************
//Program created to identify the snapshot when a Zundel
//dissociates to its new counterpart water - PT. 
//
//Input: #-Zundel-switches.txt
//Output: total-PT-events.txt
//
//
//Created by Lelee Ounkham 
//Last Modified March 24, 2018
//
//**************************************************************

int main(){

   int Hyd,O1,O2,startsnap,finalsnap;
   int sharedH,origOindex,nextOindex,timeformed,finaltime;
   int row,nrow,y,z;
   int proton,Hmax,lines;
   FILE *input,*output;
   char ifile[100];

   int eventmatrix[70000][5];

//Change variable accordingly
   Hmax=314;
   lines=0;

//Read in input file
   for(proton=103; proton<=Hmax; proton++){

      sprintf(ifile,"%d-Zundel-switches.txt",proton);
      input=fopen(ifile,"r");

      output=fopen("total-PT-events.txt","a");

       while(fscanf(input,"%d %d %d %d %d\n",&Hyd,&O1,&O2,&startsnap,&finalsnap)==5){
	  eventmatrix[lines][0]=Hyd;
	  eventmatrix[lines][1]=O1;
	  eventmatrix[lines][2]=O2;
	  eventmatrix[lines][3]=startsnap;
	  eventmatrix[lines][4]=finalsnap;
	  lines++;
       }//end while-loop
       fclose(input);

//**************************************************************
//Section I: Store info for the first Zundel that appears to 
//that corresponding H atom. Determine whether the proton
//was transferred or returned to its original partner.
//
//**************************************************************

   for(row=0; row<lines-1; row++){
     if(row==0){ //store info on the first Zundel formed
	sharedH=eventmatrix[row][0];
        origOindex=eventmatrix[row][1];
	nextOindex=eventmatrix[row][2];
        timeformed=eventmatrix[row][3];
     }
     for(nrow=1; nrow<lines; nrow++){
        if(eventmatrix[row][0]!=0 && eventmatrix[nrow][0]!=0){
           if(origOindex==eventmatrix[nrow][1] && nextOindex==eventmatrix[nrow][2]){
//If orig. and next O index matches the next event - H returned to orig. O, no PT
	      finaltime=eventmatrix[nrow][4]; //update final snap
	      eventmatrix[nrow][0]=0; 
	    }
	    else if(origOindex==eventmatrix[nrow][2] && nextOindex==eventmatrix[nrow][1]){
//When PT occurs, orig. and next O index will be switched - same Zundel though
	       if(eventmatrix[nrow+1][1]!=nextOindex && eventmatrix[nrow+1][2]!=origOindex){
//If the next event is NOT the PT Zundel (i.e. O2-O38 to O38-O2)
//Print PT which will be the snapshot after it dissociated into an Eigen
	          fprintf(output,"%d %d %d %d\n",
	          sharedH,
	          eventmatrix[nrow][1],
	          eventmatrix[nrow][2],
	          finaltime+1); //Add 1 because the dissociation occurred in the next snapshot
	 
//if(eventmatrix[nrow][1]==51){
//printf("%d %d %d %d\n",eventmatrix[nrow][1],eventmatrix[nrow][2],eventmatrix[nrow][4],eventmatrix[nrow][5]);
//}
	          finaltime=eventmatrix[nrow][4]; //update final snap
	          eventmatrix[nrow][0]=0;
	       }
	    }
	    else if(origOindex!=eventmatrix[nrow][1]  && origOindex==eventmatrix[nrow][2]){
//When the orig. index no longer appears in the Zundel after PT
//Print out previous Zundel before updated variables.
	       sharedH=eventmatrix[nrow][0];
	       origOindex=eventmatrix[nrow][1]; 
	       nextOindex=eventmatrix[nrow][2];
	       timeformed=eventmatrix[nrow][3];
	       finaltime=eventmatrix[nrow][4];

	       eventmatrix[nrow][0]=0; 
	    }
            else if(origOindex==eventmatrix[nrow][1] && nextOindex!=eventmatrix[nrow][2]){
//When PT is not transferred and a new Zundel is formed between the orig. O and a new water - partner switch
//Print out previous Zundel before updated variables.
	       sharedH=eventmatrix[nrow][0];
	       origOindex=eventmatrix[nrow][1]; 
	       nextOindex=eventmatrix[nrow][2];
	       timeformed=eventmatrix[nrow][3];
	       finaltime=eventmatrix[nrow][4];

	       eventmatrix[nrow][0]=0; 
	   }
            else if(origOindex!=eventmatrix[nrow][1] && origOindex!=eventmatrix[nrow][2] && nextOindex==eventmatrix[nrow][1]){
//When PT is tranferred to new water and then forms a new Zundel pair - partner switch
//Print out previous Zundel before updated variables.
	       sharedH=eventmatrix[nrow][0];
	       origOindex=eventmatrix[nrow][1]; 
	       nextOindex=eventmatrix[nrow][2];
	       timeformed=eventmatrix[nrow][3];
	       finaltime=eventmatrix[nrow][4];

	       eventmatrix[nrow][0]=0; 
	   }
	}
     }//end nrow-loop
   }//end row-loop
  
//Intialize eventmatrix
      for(y=0; y<30000; y++){
         for(z=0; z<5; z++){
	    eventmatrix[y][z]=0;
	 }
      }
      lines=0;
   }//end proton-loop
   fclose(output);
}//end main-loop
