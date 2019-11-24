#include<stdio.h>
#include<stdlib.h>

//*********************************************************************************
//Program created to apply select O-H*-O cutoff to Zundels to correct for
//bonds that are transiently broken or when a Zundel disappears and later returns
//and is counted new entity with a long lifetime. 
//
//Note: Even after PT, if the same O's are sharing the same proton then it
//the pair remains unchanged. O1 and O2 does NOT correspond to which O the
//proton is currently bound to.
//
//Created by Lelee Ounkham
//Last modified August 15, 2018
//
//*********************************************************************************

int main(){

   int proton,Hmax,Hyd,O1,O2,startsnap,finalsnap;
   int sharedH,origOindex,nextOindex,timeformed,finaltime;
   int lines,x,y,z,row,nrow,diff,cutoff,mlines,mrow; 
   FILE *input,*output;
   char ifile[100];

   int eventmatrix[1000][5];
   int mainmatrix[1000][5];

//Change variable accordingly
   Hmax=314;
   cutoff=100; //snapshots

//Initialize matrices and variables
   for(x=0; x<1000; x++){
      for(y=0; y<5; y++){
	 eventmatrix[x][y]=0;
      }
   }
   for(x=0; x<1000; x++){
      for(y=0; y<5; y++){
	 mainmatrix[x][y]=0;
      }
   }

   lines=0;
   mlines=0;

 //Read input file and open output file

   output=fopen("complete-list-Zundels-using-100snap-cutoff.txt","a"); 
   for(proton=103; proton<=Hmax; proton++){
//      sprintf(ifile,"%d-example.txt",proton);
      sprintf(ifile,"%d-Zundel-switches.txt",proton);
      input=fopen(ifile,"r");

//*********************************************************************************
//Section I: Store info for the first Zundel that appears to that corresponding
//H atom. Determine whether the proton was transferred or returned to its original 
//partner.
//*********************************************************************************

       while(fscanf(input,"%d %d %d %d %d\n",&Hyd,&O1,&O2,&startsnap,&finalsnap)==5){
          eventmatrix[lines][0]=Hyd;
          eventmatrix[lines][1]=O1;
          eventmatrix[lines][2]=O2;
          eventmatrix[lines][3]=startsnap;
          eventmatrix[lines][4]=finalsnap;
          lines++;
       }//end while-loop
       fclose(input);

       for(row=0; row<lines; row++){
	  if(row==0){ //Store info of first Zundel into respective variables
             sharedH=eventmatrix[row][0];
             origOindex=eventmatrix[row][1];
             nextOindex=eventmatrix[row][2];
             timeformed=eventmatrix[row][3];
             finaltime=eventmatrix[row][4];
	  }
          for(nrow=1; nrow<lines; nrow++){
             if(eventmatrix[row][0]!=0 && eventmatrix[nrow][0]!=0){
                if(origOindex==eventmatrix[nrow][1] && nextOindex==eventmatrix[nrow][2]){       
//if orig. O and next O matches the next event - H returned to orig O, no PT
//printf("nrow: %d %d %d %d %d\n",sharedH,origOindex,nextOindex,timeformed,finaltime); 
		    if(origOindex==eventmatrix[nrow][1] && nextOindex==eventmatrix[nrow][2]){
//printf("Return: %d %d %d %d %d\n",sharedH,origOindex,nextOindex,eventmatrix[nrow][3],eventmatrix[nrow][4]); 

//If the next event corresponds to the SAME H index and current O, examine time difference
		       diff=eventmatrix[nrow][3]-finaltime;

//printf("Diff 1: %d %d %d %d %d %d\n",sharedH,origOindex,nextOindex,eventmatrix[nrow][3],finaltime,diff); 
		       if(diff<=cutoff){
                         finaltime=eventmatrix[nrow][4]; //if diff is less than cutoff, bond NOT broken
//printf("bond: %d %d %d %d\n",nrow,eventmatrix[nrow][3],eventmatrix[nrow][4],diff);
 
                         mainmatrix[mlines-1][4]=eventmatrix[nrow][4];
                         eventmatrix[nrow][0]=0; 
		       }
		       else if(diff>cutoff){ //if diff is greater than cutoff, bond is break - new Zundel
//printf("broken: %d %d %d %d\n",nrow,eventmatrix[nrow][3],eventmatrix[nrow][4],diff);
 			  if(nrow==1){ //for the second line
                             mainmatrix[mlines][0]=sharedH;
                             mainmatrix[mlines][1]=origOindex;
                             mainmatrix[mlines][2]=nextOindex;
                             mainmatrix[mlines][3]=timeformed;
                             mainmatrix[mlines][4]=finaltime;
                             mlines++;
			  }
                          mainmatrix[mlines][0]=sharedH;
                          mainmatrix[mlines][1]=origOindex;
                          mainmatrix[mlines][2]=nextOindex;
                          mainmatrix[mlines][3]=eventmatrix[nrow][3];
                          mainmatrix[mlines][4]=eventmatrix[nrow][4];
                          mlines++;
			  
//printf("Main 1: %d %d %d %d %d\n",sharedH,origOindex,nextOindex,eventmatrix[nrow][3],eventmatrix[nrow][4]); 

//Copy new Zundel before updated variables.
                          sharedH=eventmatrix[nrow][0];
                          origOindex=eventmatrix[nrow][1];
                          nextOindex=eventmatrix[nrow][2];
                          timeformed=eventmatrix[nrow][3];
                          finaltime=eventmatrix[nrow][4];

                          eventmatrix[nrow][0]=0;
		       }
		    }  
	         }
		 else if(origOindex==eventmatrix[nrow][2] && nextOindex==eventmatrix[nrow][1]){
//If PT occurs, it will be bonded to nextO instead of O1 - same pair though
//If the next event corresponds to the SAME H index and current O, examine time difference
		    diff=eventmatrix[nrow][3]-finaltime;

//printf("Diff 1: %d %d %d %d %d %d\n",sharedH,origOindex,nextOindex,eventmatrix[nrow][3],finaltime,diff); 
		    if(diff<=cutoff){
                       finaltime=eventmatrix[nrow][4]; //if diff is less than cutoff, bond NOT broken
//printf("bond: %d %d %d %d\n",nrow,eventmatrix[nrow][3],eventmatrix[nrow][4],diff);
 
                       mainmatrix[mlines-1][4]=eventmatrix[nrow][4];
                       eventmatrix[nrow][0]=0; 
		    }
		    else if(diff>cutoff){ //if diff is greater than cutoff, bond is break - new Zundel
//printf("broken: %d %d %d %d\n",nrow,eventmatrix[nrow][3],eventmatrix[nrow][4],diff);
 		       if(nrow==1){ //for the second line
                          mainmatrix[mlines][0]=sharedH;
                          mainmatrix[mlines][1]=origOindex;
                          mainmatrix[mlines][2]=nextOindex;
                          mainmatrix[mlines][3]=timeformed;
                          mainmatrix[mlines][4]=finaltime;
                          mlines++;

                          eventmatrix[nrow][0]=0;
		       }
                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=eventmatrix[nrow][3];
                       mainmatrix[mlines][4]=eventmatrix[nrow][4];
                       mlines++;
			  
//Copy new Zundel before updated variables.
                       sharedH=eventmatrix[nrow][0];
                       origOindex=eventmatrix[nrow][2];
                       nextOindex=eventmatrix[nrow][1];
                       timeformed=eventmatrix[nrow][3];
                       finaltime=eventmatrix[nrow][4];

                       eventmatrix[nrow][0]=0;
		    }	    
		 }
		 else if(origOindex!=eventmatrix[nrow][1] && origOindex!=eventmatrix[nrow][2] && nextOindex==eventmatrix[nrow][1]){
//Partner switch case 1: when origO switches partner to form new Zundel with nextO
                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=timeformed;
                       mainmatrix[mlines][4]=finaltime;
                       mlines++;

//Copy new Zundel before updated variables.
                       sharedH=eventmatrix[nrow][0];
                       origOindex=eventmatrix[nrow][1];
                       nextOindex=eventmatrix[nrow][2];
                       timeformed=eventmatrix[nrow][3];
                       finaltime=eventmatrix[nrow][4];

                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=timeformed;
                       mainmatrix[mlines][4]=finaltime;
                       mlines++;

//printf("trigger 1: %d %d %d %d %d\n",sharedH,origOindex,nextOindex,timeformed,finaltime);

                       eventmatrix[nrow][0]=0;

		 }
		 else if(origOindex!=eventmatrix[nrow][1] && origOindex!=eventmatrix[nrow][2] && nextOindex==eventmatrix[nrow][2]){
//Partner switch case 2: when origO switches partner to form new Zundel with nextO

                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=timeformed;
                       mainmatrix[mlines][4]=finaltime;
                       mlines++;
			  
//Copy new Zundel before updated variables.
                       sharedH=eventmatrix[nrow][0];
                       origOindex=eventmatrix[nrow][1];
                       nextOindex=eventmatrix[nrow][2];
                       timeformed=eventmatrix[nrow][3];
                       finaltime=eventmatrix[nrow][4];

                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=timeformed;
                       mainmatrix[mlines][4]=finaltime;
                       mlines++;

//printf("trigger 2: %d %d %d %d %d\n",sharedH,origOindex,nextOindex,timeformed,finaltime);

                       eventmatrix[nrow][0]=0;

		 }
		 else if(nextOindex!=eventmatrix[nrow][1] && nextOindex!=eventmatrix[nrow][2] && origOindex==eventmatrix[nrow][1]){
//Partner switch case 3: when newO switches partner to form new Zundel with origO
                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=timeformed;
                       mainmatrix[mlines][4]=finaltime;
                       mlines++;
			  
//Copy new Zundel before updated variables.
                       sharedH=eventmatrix[nrow][0];
                       origOindex=eventmatrix[nrow][1];
                       nextOindex=eventmatrix[nrow][2];
                       timeformed=eventmatrix[nrow][3];
                       finaltime=eventmatrix[nrow][4];

                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=timeformed;
                       mainmatrix[mlines][4]=finaltime;
                       mlines++;

//printf("trigger 3: %d %d %d %d %d\n",sharedH,origOindex,nextOindex,timeformed,finaltime);
                       eventmatrix[nrow][0]=0;

		 }
		 else if(nextOindex!=eventmatrix[nrow][1] && nextOindex!=eventmatrix[nrow][2] && origOindex==eventmatrix[nrow][2]){
//Partner switch case 4: when newO switches partner to form new Zundel with origO
                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=timeformed;
                       mainmatrix[mlines][4]=finaltime;
                       mlines++;
			  
//Copy new Zundel before updated variables.
                       sharedH=eventmatrix[nrow][0];
                       origOindex=eventmatrix[nrow][1];
                       nextOindex=eventmatrix[nrow][2];
                       timeformed=eventmatrix[nrow][3];
                       finaltime=eventmatrix[nrow][4];

                       mainmatrix[mlines][0]=sharedH;
                       mainmatrix[mlines][1]=origOindex;
                       mainmatrix[mlines][2]=nextOindex;
                       mainmatrix[mlines][3]=timeformed;
                       mainmatrix[mlines][4]=finaltime;
                       mlines++;

//printf("trigger 4: %d %d %d %d %d\n",sharedH,origOindex,nextOindex,timeformed,finaltime);

                       eventmatrix[nrow][0]=0;

		 }
	      }
	   }//end nrow-loop
	}//end row-loop

            for(z=0; z<mlines; z++){
               if(mainmatrix[z][0]!=0){
                  fprintf(output,"%d %d %d %d %d\n",
                  mainmatrix[z][0],
                  mainmatrix[z][1],
                  mainmatrix[z][2],
                  mainmatrix[z][3],
                  mainmatrix[z][4]);
               }
            }//end z-loop

//Initialize eventmatrix and mainmatrix for the next H atom
	   for(x=0; x<1000; x++){
              for(y=0; y<5; y++){
	         eventmatrix[x][y]=0;
              } 
           }
	   lines=0;

           for(x=0; x<1000; x++){
              for(y=0; y<5; y++){
	          mainmatrix[x][y]=0;
              }
           }
	   mlines=0;
   }//end proton-loop 
   fclose(output);

}//end main-loop
