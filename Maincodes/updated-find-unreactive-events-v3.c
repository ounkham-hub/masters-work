#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//******************************************************************************
//Program created to identify H atoms, associated with water molecuels, that
//remain in a stable state for a duration in time. This requires all 32
//beads contained with the ring polymer to remain unstable too. 
//
//This is an extension of a previous program that did not fully account for
//potential fluctuations in connectivity.
//
//Created by Lelee Ounkham
//Last modified October 24, 2018
//
//******************************************************************************

int main(){

   int Hindex,Oindex,flag,initialtime,finaltime;
   int Hatom,mlines,repno,coln,var;
   int x,y,z,clines,lines,bead,diff;
   int crow,row,mrow,time;
   FILE *input,*houtput,*newinput,*output;
   char ifile[100],hfile[100],nfile[100];

   int centroidmatrix[195000][6][33];
   int mainlist[1000000][6][33];
   int currentbead[10000][5];

   time=20; //snapshot window size

/*
   for(x=0; x<1000000; x++){
      for(y=0; y<7; y++){
         for(z=0; z<33; z++){
            centroidmatrix[x][y][z]=0;
         }
      }
   }
   clines=0;

      for(bead=1; bead<=32; bead++){
         sprintf(ifile,"%d-undynamic-events-diff.txt",bead);
         input=fopen(ifile,"r");

//Read in input file
         while(fscanf(input,"%d %d %d %d %d %d\n",&Hindex,&Oindex,&flag,&initialtime,&finaltime,&diff)==6){
            if(diff>1){
               centroidmatrix[clines][0][0]=Hindex;
               centroidmatrix[clines][1][0]=Oindex;
               centroidmatrix[clines][2][0]=initialtime;
               centroidmatrix[clines][3][0]=finaltime;
               centroidmatrix[clines][4][0]=diff;
               centroidmatrix[clines][5][0]=0;
               clines++;
//printf("%d %d %d %d %d\n",Hindex,Oindex,initialtime,finaltime,diff); 
             }
          }//end while-loop
          fclose(input);

  for(crow=0; crow<clines; crow++){
      printf("%d %d %d %d %d %d\n",
      centroidmatrix[crow][0][0],
      centroidmatrix[crow][1][0],
      centroidmatrix[crow][2][0],
      centroidmatrix[crow][3][0],
      centroidmatrix[crow][4][0],
      centroidmatrix[crow][5][0]);      
   }//end crow-loop


//    printf("rlines tlines: %d %d\n",rlines,tlines);
//    printf("%d %d\n",bead,rlines);

//*******************************************************************************
//Section to separate connections based on H index and output into a independent
//files organized by OH initial times.
//*******************************************************************************

      for(Hatom=212; Hatom<=314; Hatom++){
         sprintf(hfile,"%d-H%d-persistences.txt",bead,Hatom);
         houtput=fopen(hfile,"a");
 
         for(crow=0; crow<clines; crow++){
            if(centroidmatrix[crow][0][0]==Hatom){ //if H index matches store connection into new file
	       fprintf(houtput,"%d %d %d %d %d\n",
	       centroidmatrix[crow][0][0], //Hindex
	       centroidmatrix[crow][1][0], //Oindex
	       centroidmatrix[crow][2][0], //initial snapshot
	       centroidmatrix[crow][3][0], //final snapshot
	       centroidmatrix[crow][4][0]); //diff
	   }//end row-loop
        }//end crow-loop
     }//end Hatom-loop

   for(x=0; x<195000; x++){
      for(y=0; y<7; y++){
         for(z=0; z<33; z++){
            centroidmatrix[x][y][z]=0;
         }
      }
   }
   clines=0;
  }//end bead-loop 
  fclose(houtput);
*/
//*******************************************************************************
//Section to compare connections in each replica and determine the windowsize
//when all 32 beads agree. For every Hatom exam persistences in each of the
//replicas.
//*******************************************************************************

    output=fopen("All-32-agreed-connections-20snap-v3.txt","a");

     for(x=0; x<1000000; x++){
        for(y=0; y<7; y++){
           for(z=0; z<33; z++){
             mainlist[x][y][z]=0;
           }
        }
     }
     mlines=0;

     for(x=0; x<10000; x++){
        for(y=0; y<5; y++){
            currentbead[x][y]=0;
        }
     }
     lines=0;
//Let's make element mainlist[mlines][0][1] the counter element

     for(Hatom=212; Hatom<=314; Hatom++){
        for(bead=1; bead<=32; bead++){
           sprintf(nfile,"%d-H%d-persistences.txt",bead,Hatom);
           newinput=fopen(nfile,"r");

	   if(bead==1){ //for bead #1, store info
	      while(fscanf(newinput,"%d %d %d %d %d\n",&Hindex,&Oindex,&initialtime,&finaltime,&diff)==5){
                 mainlist[mlines][0][0]=Hindex;
                 mainlist[mlines][1][0]=Oindex;
                 mainlist[mlines][2][0]=initialtime;
                 mainlist[mlines][3][0]=finaltime;
                 mainlist[mlines][4][0]=diff;
		 mainlist[mlines][5][0]=0;

		 mainlist[mlines][0][1]=1; //start count at 1

		 mainlist[mlines][2][1]=initialtime;
		 mainlist[mlines][3][1]=finaltime;
		 mainlist[mlines][4][1]=diff;
                 mlines++;
//printf("%d %d %d %d %d\n",Hindex,Oindex,initialtime,finaltime,diff);

              }//end while-loop
              fclose(newinput);

//printf("1) %d %d %d\n",bead,Hatom,mlines);

	   }
	   else if(bead>1){
	      while(fscanf(newinput,"%d %d %d %d %d\n",&Hindex,&Oindex,&initialtime,&finaltime,&diff)==5){
                 currentbead[lines][0]=Hindex;
                 currentbead[lines][1]=Oindex;
                 currentbead[lines][2]=initialtime;
                 currentbead[lines][3]=finaltime;
                 currentbead[lines][4]=diff;
                 lines++;
//printf("trigger: %d %d %d %d %d\n",Hindex,Oindex,initialtime,finaltime,diff); 


	      }//end while-loop
	      fclose(newinput);
/*
        for(mrow=0; mrow<mlines; mrow++){
	   printf("%d %d %d %d %d\n",
	   mainlist[mrow][0][0],
	   mainlist[mrow][1][0],
           mainlist[mrow][2][0],
           mainlist[mrow][3][0],
           mainlist[mrow][4][0]);
	}//end mrow-loop 
*/
//*******************************************************************************
//Section to compare current replica to next replica and store into matrix
//*******************************************************************************
//printf("%d %d %d\n",bead,lines,mlines);


	      if(lines >= mlines){ //to ensure all lines are being compared in matrix properly
	         for(crow=0; crow<lines; crow++){
		    for(mrow=0; mrow<mlines; mrow++){
			if(mainlist[mrow][5][0]==0 && currentbead[crow][0]!=0){// only examine nonzero elements
			   if(mainlist[mrow][0][0]==currentbead[crow][0] && mainlist[mrow][1][0]==currentbead[crow][1]){
//If O and H indices matches in both replicas
			      if(currentbead[crow][2] >= mainlist[mrow][2][0]-time && currentbead[crow][2] <= mainlist[mrow][2][0]+time && currentbead[crow][3] >= mainlist[mrow][3][0]-time && currentbead[crow][3] <= mainlist[mrow][3][0]+time){
//Initial and final snapshot must be within bounds of current bead replica, store data into currentbead matrix
			         repno=bead;
			         mainlist[mrow][2][repno]=currentbead[crow][2];
			         mainlist[mrow][3][repno]=currentbead[crow][3];
				 mainlist[mrow][4][repno]=currentbead[crow][4];
				 mainlist[mrow][5][0]=1; //turn switch on

 				 mainlist[mrow][0][1]=mainlist[mrow][0][1]+1; //increment counter
				 
				 currentbead[crow][0]=0; //zero out element
			      }
			   }
		        }//end mrow-loop
	             }//end crow-loop
		  }
	       }//end if lines > clines
	       else if(mlines >= lines){
	         for(mrow=0; mrow<mlines; mrow++){
		    for(crow=0; crow<lines; crow++){
			if(mainlist[mrow][5][0]==0 && currentbead[crow][0]!=0){// only examine nonzero elements
			   if(mainlist[mrow][0][0]==currentbead[crow][0] && mainlist[mrow][1][0]==currentbead[crow][1]){
//If O and H indices matches in both replicas
			      if(currentbead[crow][2] >= mainlist[mrow][2][0]-time && currentbead[crow][2] <= mainlist[mrow][2][0]+time && currentbead[crow][3] >= mainlist[mrow][3][0]-time && currentbead[crow][3] <= mainlist[mrow][3][0]+time){
//Initial and final snapshot must be within bounds of current bead replica, store data into currentbead matrix
			         repno=bead;
			         mainlist[mrow][2][repno]=currentbead[crow][2];
			         mainlist[mrow][3][repno]=currentbead[crow][3];
				 mainlist[mrow][4][repno]=currentbead[crow][4];
				 mainlist[mrow][5][0]=1; //turn switch on

 				 mainlist[mrow][0][1]=mainlist[mrow][0][1]+1; //increment counter

				 currentbead[crow][0]=0; //zero out element
			      }
			   }
		        }
	             }//end mrow-loop
		  }//end crow-loop
	       }//else if clines > lines
//*******************************************************************************
//Section to add new connections for next iteration 
//*******************************************************************************

	   for(crow=0; crow<lines; crow++){
	      if(currentbead[crow][0]!=0){ //for unmarked connections
		 repno=bead;
		 mainlist[mlines][0][0]=currentbead[crow][0];
		 mainlist[mlines][1][0]=currentbead[crow][1];
		 mainlist[mlines][2][0]=currentbead[crow][2];
		 mainlist[mlines][3][0]=currentbead[crow][3];
		 mainlist[mlines][4][0]=currentbead[crow][4];
                 mainlist[mlines][2][repno]=currentbead[crow][2];
		 mainlist[mlines][3][repno]=currentbead[crow][3];
		 mainlist[mlines][4][repno]=currentbead[crow][4];
		 mlines++;
		 currentbead[crow][0]=0; //zero out element
	      }
	   }

	   }//end elseif-bead>1 loop

//*******************************************************************************
//Section to initialize currentbead matrix
//*******************************************************************************
              for(x=0; x<10000; x++){
                 for(y=0; y<7; y++){
                     currentbead[x][y]=0;
                 }
             }
             lines=0;


	     for(x=0; x<mlines; x++){
                mainlist[x][5][0]=0;
             }

//*******************************************************************************
//Section to print out mainmatrix containing list of persistences and how many
//beads agree using a +/- 20 windowtime
//*******************************************************************************

//         if(Hatom==314){
	    for(mrow=0; mrow<mlines; mrow++){
	       var=1;
	    for(coln=1; coln<33; coln++){
	       var=var*mainlist[mrow][2][coln];
	    }
	    if(var!=0 && var!=1){
printf("%d\n",mainlist[mrow][0][1]); //print counter, should always be 32
		     
	             fprintf(output,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	             mainlist[mrow][0][0],
	             mainlist[mrow][1][0],
	             mainlist[mrow][2][0],
	             mainlist[mrow][3][0],
	             mainlist[mrow][4][0],
	             mainlist[mrow][2][1],
	             mainlist[mrow][2][2],
                     mainlist[mrow][2][3],
	             mainlist[mrow][2][4],
	             mainlist[mrow][2][5],
	             mainlist[mrow][2][6],
	             mainlist[mrow][2][7],
                     mainlist[mrow][2][8],
	             mainlist[mrow][2][9],
	             mainlist[mrow][2][10],
	             mainlist[mrow][2][11],
	             mainlist[mrow][2][12],
                     mainlist[mrow][2][13],
	             mainlist[mrow][2][14],
	             mainlist[mrow][2][15],
	             mainlist[mrow][2][16],
	             mainlist[mrow][2][17],
                     mainlist[mrow][2][18],
	             mainlist[mrow][2][19],
	             mainlist[mrow][2][20],
	             mainlist[mrow][2][21],
	             mainlist[mrow][2][22],
	             mainlist[mrow][2][23],
                     mainlist[mrow][2][24],
	             mainlist[mrow][2][25],
	             mainlist[mrow][2][26],
	             mainlist[mrow][2][27],
	             mainlist[mrow][2][28],
                     mainlist[mrow][2][29],
	             mainlist[mrow][2][30],
	             mainlist[mrow][2][31],
	             mainlist[mrow][2][32]);
	         }
	      }//end mrow-loop
           }   
//        }//end bead-loop

         for(x=0; x<mlines; x++){
           for(y=0; y<7; y++){
              for(z=0; z<33; z++){
                 mainlist[x][y][z]=0;
              }
           }
        }
        mlines=0;


    }//end Hatom-loop
    fclose(output);

}//end main-loop

