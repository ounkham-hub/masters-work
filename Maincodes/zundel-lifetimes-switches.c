#include<stdio.h>
#include<stdlib.h>

//**************************************************************
//Program created to find the initial formation of Zundels
//and when they completely dissociate or form a new Zundel.
//As a crosscheck, the # of oscillations between Zundels
//will also be counted. 
//
//Input: #-Zundel-switches.txt
//Output:
//
//
//Created by Lelee Ounkham 
//Last Modified March 21, 2018
//
//**************************************************************

int main(){

   int Hyd,O1,O2,startsnap,finalsnap;
   int sharedH,origOindex,nextOindex,timeformed,finaltime;
   int row,nrow,y,z;
   int proton,Hmax,lines,mlines,switches;
   FILE *input,*output,*output2,*PToutput;
   char ifile[100];

   int eventmatrix[30000][5];
   int mainmatrix[30000][6];

//Change variable accordingly
   Hmax=314;
   lines=0;
   switches=0;
   mlines=0;

//Read in input file
   for(proton=103; proton<=Hmax; proton++){
      sprintf(ifile,"%d-Zundel-switches.txt",proton);
      input=fopen(ifile,"r");


//Initialize matrix for the next proton
       for(y=0; y<lines; y++){
          for(z=0; z<5; z++){
	     eventmatrix[y][z]=0;
	  }
       }
       lines=0;
       switches=0;

//This is only for the example input file - change input as needed
       output=fopen("lifetimes-zundels.txt","w");
//       output2=fopen("partner-switch-list.txt","a");
       PToutput=fopen("all-H-PT-events.txt","a");

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

   for(row=0; row<lines; row++){
     if(row==0){ //store info on the first Zundel formed
	sharedH=eventmatrix[row][0];
        origOindex=eventmatrix[row][1];
	nextOindex=eventmatrix[row][2];
        timeformed=eventmatrix[row][3];
     }
     for(nrow=0; nrow<lines; nrow++){
        if(eventmatrix[row][0]!=0 && eventmatrix[nrow][0]!=0){
           if(origOindex==eventmatrix[nrow][1] && nextOindex==eventmatrix[nrow][2]){
//If orig. and next O index matches the next event - H returned to orig. O, no PT
	      finaltime=eventmatrix[nrow][4]; //update final snap
	      eventmatrix[nrow][0]=0; 

	    }
	    else if(origOindex==eventmatrix[nrow][2] && nextOindex==eventmatrix[nrow][1]){
//When PT occurs, orig. and next O index will be switched - same Zundel though
//Print out PT occurrence 
	       fprintf(PToutput,"%d %d %d %d\n",
	       eventmatrix[row][0],
	       eventmatrix[row][1],
	       eventmatrix[row][2],
	       finaltime+1);

	       finaltime=eventmatrix[nrow][4]; //update final snap
	       switches++; //increase PT count
	       eventmatrix[nrow][0]=0;
	    }
	    else if(origOindex!=eventmatrix[nrow][1]  && origOindex==eventmatrix[nrow][2] && nextOindex!=eventmatrix[nrow][1]){
//When the orig. index no longer appears in the Zundel after PT
	       mainmatrix[mlines][0]=sharedH;
	       mainmatrix[mlines][1]=origOindex;
	       mainmatrix[mlines][2]=nextOindex;
	       mainmatrix[mlines][3]=timeformed;
	       mainmatrix[mlines][4]=finaltime;
	       mainmatrix[mlines][5]=switches;
	       mlines++;

//Partner switch - print info.
/*               fprintf(output2,"%d %d %d %d %d %d %d\n",
	       eventmatrix[nrow][0],
	       origOindex,
	       nextOindex,
	       eventmatrix[nrow][1],
	       eventmatrix[nrow][2],
	       eventmatrix[nrow][3],
	       eventmatrix[nrow][4]); //snap when new Zundel formed
*/
//Print out previous Zundel before updated variables.
	       sharedH=eventmatrix[nrow][0];
	       origOindex=eventmatrix[nrow][1]; 
	       nextOindex=eventmatrix[nrow][2];
	       timeformed=eventmatrix[nrow][3];
	       finaltime=eventmatrix[nrow][4];
	       switches=0;

	       eventmatrix[nrow][0]=0; 
	    }
            else if(origOindex==eventmatrix[nrow][1] && nextOindex!=eventmatrix[nrow][2]){
//When PT is not transferred and a new Zundel is formed between the orig. O and a new water - partner switch
	       mainmatrix[mlines][0]=sharedH;
	       mainmatrix[mlines][1]=origOindex;
	       mainmatrix[mlines][2]=nextOindex;
	       mainmatrix[mlines][3]=timeformed;
	       mainmatrix[mlines][4]=finaltime;
	       mainmatrix[mlines][5]=switches;
	       mlines++;

//Partner switch - print info.
/*               fprintf(output2,"%d %d %d %d %d %d %d\n",
	       eventmatrix[nrow][0],
	       origOindex,
	       nextOindex,
	       eventmatrix[nrow][1],
	       eventmatrix[nrow][2],
	       eventmatrix[nrow][3],
	       eventmatrix[nrow][4]); //snap when new Zundel formed
*/

//Print out previous Zundel before updated variables.
	       sharedH=eventmatrix[nrow][0];
	       origOindex=eventmatrix[nrow][1]; 
	       nextOindex=eventmatrix[nrow][2];
	       timeformed=eventmatrix[nrow][3];
	       finaltime=eventmatrix[nrow][4];
	       switches=0;

	       eventmatrix[nrow][0]=0; 
	   }
            else if(origOindex!=eventmatrix[nrow][1] && origOindex!=eventmatrix[nrow][2] && nextOindex==eventmatrix[nrow][1]){
//When PT is tranferred to new water and then forms a new Zundel pair - partner switch
	       mainmatrix[mlines][0]=sharedH;
	       mainmatrix[mlines][1]=origOindex;
	       mainmatrix[mlines][2]=nextOindex;
	       mainmatrix[mlines][3]=timeformed;
	       mainmatrix[mlines][4]=finaltime;
	       mainmatrix[mlines][5]=switches;
	       mlines++;

//Partner switch - print info.
/*               fprintf(output2,"%d %d %d %d %d %d %d\n",
	       eventmatrix[nrow][0],
	       origOindex,
	       nextOindex,
	       eventmatrix[nrow][1],
	       eventmatrix[nrow][2],
	       eventmatrix[nrow][3],
	       eventmatrix[nrow][4]); //snap when new Zundel formed
*/

//Print out previous Zundel before updated variables.
	       sharedH=eventmatrix[nrow][0];
	       origOindex=eventmatrix[nrow][1]; 
	       nextOindex=eventmatrix[nrow][2];
	       timeformed=eventmatrix[nrow][3];
	       finaltime=eventmatrix[nrow][4];
	       switches=0;

	       eventmatrix[nrow][0]=0; 
	   }
	}
     }//end nrow-loop
   }//end row-loop

//Print out last Zundel on list
    mainmatrix[mlines][0]=sharedH;
    mainmatrix[mlines][1]=origOindex;
    mainmatrix[mlines][2]=nextOindex;
    mainmatrix[mlines][3]=timeformed;
    mainmatrix[mlines][4]=finaltime;
    mainmatrix[mlines][5]=switches;
    mlines++;

   for(z=0; z<mlines; z++){
      if(mainmatrix[z][0]!=0){
         fprintf(output,"%d %d %d %d %d %d\n",
         mainmatrix[z][0],
         mainmatrix[z][1],
         mainmatrix[z][2],
         mainmatrix[z][3],
         mainmatrix[z][4],
         mainmatrix[z][5]);
      }
   }//end z-loop

   fclose(output);
//   fclose(output2);
   fclose(PToutput);
   }//end proton-loop

}//end main-loop
