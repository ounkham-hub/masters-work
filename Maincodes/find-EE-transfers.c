#include<stdio.h>
#include<stdlib.h>

//*********************************************************************
//Program created to track a specific proton in time by comparing the 
//change in structure type in current snap to the next.
//
//Modified to output the O-H+ in Zundels and E-E transfers
//
//Note the flags indicate the following:
//    - 0, water
//    - 1, eigen
//    - 2, zundel
//
//Input: R-HCl.input.solB#.xyz.solC#.xyz.GraphGeod
//Output: event-list.txt
//
//Created by Lelee Ounkham
//Last modified: February 6, 2018
//
//*********************************************************************

int main(){

//Parameters to read O-H GraphGeods
   int Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5;
   int snap,maxsnap,countOindex,y,z;
   int Oindex,Hhold1,Hhold2;
   int clines,crow,coln;
   int nlines,nrow,ncoln;
   float dist,angle,disthold1,disthold2;
   FILE *input,*nextinput,*output,*output2,*output3;
   char ifile[100],nfile[100];

//Parameters for matrices
   int currentmatrix[300][7];
   float currentdist[300][5];
   int nextmatrix[300][7];
   float nextdist[300][5];
 
//Intialize matrices
         for(y=0; y<300; y++){
            for(z=0; z<8; z++){
               currentmatrix[y][z]=0;
               nextmatrix[y][z]=0;
            }//end y-loop
         }//end z-loop

         for(y=0; y<300; y++){
            for(z=0; z<6; z++){
               currentdist[y][z]=0.0;
               nextdist[y][z]=0.0;
            }//end y-loop
         }//end z-loop
         clines=0;
         countOindex=1;
  
   maxsnap=61555; //CHANGE ACCORDINGLY 
//   output=fopen("1-event-list.txt","w");
   output2=fopen("1-EE-transfers.txt","w");
   output3=fopen("cent-Zundel-dist.txt","w");
   for(snap=1; snap<=maxsnap; snap++){
      if(snap==maxsnap){ //END OF PROGRAM
	 exit(0);
      }
      else{

     if(snap==1){
         sprintf(ifile,"R-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
         input=fopen(ifile,"r");

//*******************************************************
//Section I(a): Identify which O's are a part of a 
//water OR Eigen structure for current snapshot.
//
//*******************************************************
         while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
              if(countOindex==1){
                 Oindex=Oxy;
                 Hhold1=Hyd;
                 disthold1=dist;
              }
              else if(countOindex==2){
                 if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;
		    disthold2=dist;

    	            currentmatrix[clines][0]=snap;
		    currentmatrix[clines][1]=Oxy;
		    currentmatrix[clines][2]=Hhold1;
		    currentmatrix[clines][3]=Hhold2;
		    currentdist[clines][2]=disthold1;
		    currentdist[clines][3]=disthold2;
		    currentmatrix[clines][5]=0; //flagged as water
		    
		    clines++;

                 }
                 else{ //Should never trigger unless OH- exist       
	            countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
		    disthold1=dist;
                 }
              }
              else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
                 if(Oindex==Oxy){
                    currentmatrix[clines][0]=snap;
                    currentmatrix[clines][1]=Oxy;
                    currentmatrix[clines][2]=Hhold1;
                    currentmatrix[clines][3]=Hhold2;
                    currentmatrix[clines][4]=Hyd;
		    currentdist[clines][2]=disthold1;
		    currentdist[clines][3]=disthold2;
		    currentdist[clines][4]=dist;
		    currentmatrix[clines][5]=1; //flagged as eigen or H3O+

//Remove previous data on O atom, when H3O+ was no identified
                    currentmatrix[clines-1][0]=0;  
                    currentmatrix[clines-1][1]=0;
                    currentmatrix[clines-1][2]=0;
                    currentmatrix[clines-1][3]=0;
                    currentmatrix[clines-1][4]=0;
		    currentdist[clines-1][2]=0.0;
		    currentdist[clines-1][3]=0.0;
		    currentdist[clines-1][4]=0.0;
		    clines++;		  

                    countOindex=0; //restarts Oindex to read the next line;
                 }
                 else{ //Water molecule identified and countOindex==2

                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
		    disthold1=dist;
                 }
              }
              countOindex++;
         }//end while-loop
         fclose(input);

/* 	for(z=0; z<clines; z++){
	   if(currentmatrix[z][0]!=0){
	      printf("%d %d %d %d %d %d\n",
	      currentmatrix[z][0],
	      currentmatrix[z][1],
	      currentmatrix[z][2],
	      currentmatrix[z][3],
	      currentmatrix[z][4],
	      currentmatrix[z][5]);
	   }
	}//end z-loop
*/

//*******************************************************
//Section I(b): Identify which Eigens are a part of a 
//Zundel - when two Eigens share a proton.
//
//*******************************************************
      for(crow=0; crow<clines; crow++){
         for(coln=2; coln<5; coln++){
	    for(nrow=1; nrow<clines; nrow++){
	       for(ncoln=2; ncoln<5; ncoln++){
	          if(currentmatrix[crow][coln]!=0){ //if matrix element is nonzero
                     if(currentmatrix[crow][coln]==currentmatrix[nrow][ncoln]){ //If any H's atoms match
		        if(currentmatrix[crow][1]!=currentmatrix[nrow][1]){ //and if O' atoms differ - Zundel identified
			   if(currentmatrix[crow][5]!=2 && currentmatrix[nrow][5]!=2){
		              fprintf(output3,"%d %d %d %d %f %f\n",
	      	              currentmatrix[crow][0], //snapshot #
		              currentmatrix[crow][1], //initial O
	      	              currentmatrix[nrow][1], //partner O
	      	              currentmatrix[crow][coln], //shared proton
	      	              currentdist[crow][coln], //dist. between initial O and H+
	      	              currentdist[nrow][ncoln]); //dist. between final O and H+
			   }
//Change to flag 2
		           currentmatrix[crow][5]=2;
		           currentmatrix[nrow][5]=2;
			   currentmatrix[crow][6]=currentmatrix[crow][coln]; //store shared proton
			   currentmatrix[nrow][6]=currentmatrix[nrow][ncoln]; //store shared proton
                        }
		     }
		  }
               }
	    }//end coln-loop 
	 }//end nrow-loop
      }//end crow-loop  
/* 
 	for(z=0; z<clines; z++){
	   if(currentmatrix[z][0]!=0){
	      printf("%d %d %d %d %d %d\n",
	      currentmatrix[z][0],
	      currentmatrix[z][1],
	      currentmatrix[z][2],
	      currentmatrix[z][3],
	      currentmatrix[z][4],
	      currentmatrix[z][5]);
	   }
	}//end z-loop
*/   
     }//end if snap==1 loop

//*******************************************************
//Section II(a): Identify which O's are a part of a 
//water OR Eigen structure for NEXT snapshot.
//
//*******************************************************
           sprintf(nfile,"R-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap+1,snap+1); //open next snap
           nextinput=fopen(nfile,"r");
           countOindex==1;
           nlines=0;

           while(fscanf(nextinput,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
              if(countOindex==1){
                 Oindex=Oxy;
                 Hhold1=Hyd;
		 disthold1=dist;
              }
              else if(countOindex==2){
                 if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;
		    disthold2=dist;

    	            nextmatrix[nlines][0]=snap+1;
		    nextmatrix[nlines][1]=Oxy;
		    nextmatrix[nlines][2]=Hhold1;
		    nextmatrix[nlines][3]=Hhold2;
	 	    nextdist[nlines][2]=disthold1;
		    nextdist[nlines][3]=disthold2;
		    nextmatrix[nlines][5]=0; //flagged as water
		    nlines++;

                 }
                 else{ //Should never trigger unless OH- exist       
	            countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
		    disthold1=dist;
                 }
              }
              else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
                 if(Oindex==Oxy){
                    nextmatrix[nlines][0]=snap+1;
                    nextmatrix[nlines][1]=Oxy;
                    nextmatrix[nlines][2]=Hhold1;
                    nextmatrix[nlines][3]=Hhold2;
                    nextmatrix[nlines][4]=Hyd;
		    nextdist[nlines][2]=disthold1;
		    nextdist[nlines][3]=disthold2;
		    nextdist[nlines][4]=dist;
		    nextmatrix[nlines][5]=1; //flagged as eigen or H3O+

//Remove previous data on O atom, when H3O+ was no identified
                    nextmatrix[nlines-1][0]=0;  
                    nextmatrix[nlines-1][1]=0;
                    nextmatrix[nlines-1][2]=0;
                    nextmatrix[nlines-1][3]=0;
                    nextmatrix[nlines-1][4]=0;
		    nextdist[nlines-1][2]=0.0;
		    nextdist[nlines-1][3]=0.0;
		    nextdist[nlines-1][4]=0.0;
		    nlines++;		  

                    countOindex=0; //restarts Oindex to read the next line;
                 }
                 else{ //Water molecule identified and countOindex==2

                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
		    disthold1=dist;
                 }
              }
              countOindex++;
         }//end while-loop
         fclose(nextinput);
/*
  	for(z=0; z<nlines; z++){
	   if(currentmatrix[z][0]!=0){
	      printf("%d %d %d %d %d %d\n",
	      nextmatrix[z][0],
	      nextmatrix[z][1],
	      nextmatrix[z][2],
	      nextmatrix[z][3],
	      nextmatrix[z][4],
	      nextmatrix[z][5]);
	   }
	}//end z-loop
*/


//*******************************************************
//Section II(b): Identify which Eigens are a part of a 
//Zundel - when two Eigens share a proton.
//
//*******************************************************
      for(crow=0; crow<nlines; crow++){
         for(coln=2; coln<5; coln++){
	    for(nrow=1; nrow<nlines; nrow++){
	       for(ncoln=2; ncoln<5; ncoln++){
	          if(nextmatrix[crow][coln]!=0){ //if matrix element is nonzero
                     if(nextmatrix[crow][coln]==nextmatrix[nrow][ncoln]){ //If any H's atoms match
		        if(nextmatrix[crow][1]!=nextmatrix[nrow][1]){ //and if O' atoms differ - Zundel identified
			   if(nextmatrix[crow][5]!=2 && nextmatrix[nrow][5]!=2){ //if not already flagged
		  	      fprintf(output3,"%d %d %d %d %f %f\n",
			         nextmatrix[crow][0], //snapshot #
				 nextmatrix[crow][1], //initial O
				 nextmatrix[nrow][1], //partner O
				 nextmatrix[crow][coln], //shared H
				 nextdist[crow][coln], //dist. between initial O and H+
				 nextdist[nrow][ncoln]); //dist. between final O and H+
			   }
//Change flag type to 2
		           nextmatrix[crow][5]=2;
		           nextmatrix[nrow][5]=2;
			   nextmatrix[crow][6]=nextmatrix[crow][coln];
			   nextmatrix[nrow][6]=nextmatrix[nrow][ncoln]; 
                        }
		     }
		  }
               }
	    }//end coln-loop 
	 }//end nrow-loop
      }//end crow-loop  
/* 
 	for(z=0; z<clines; z++){
	   if(currentmatrix[z][0]!=0){
	      printf("%d %d %d %d %d %d\n",
	      currentmatrix[z][0],
	      currentmatrix[z][1],
	      currentmatrix[z][2],
	      currentmatrix[z][3],
	      currentmatrix[z][4],
	      currentmatrix[z][5]);
	   }
	}//end z-loop
*/

//*******************************************************
//Section III: Track when flag or structure
//changes occur in time. 
//
//*******************************************************

        for(crow=0; crow<clines; crow++){
	   for(nrow=0; nrow<nlines; nrow++){
	      if(currentmatrix[crow][0]!=0){ //if nonzero matrix elements
		 if(currentmatrix[crow][1]==nextmatrix[nrow][1]){ //if O's in currentmatrix match nextmatrix
		    if(currentmatrix[crow][5]!=nextmatrix[nrow][5]){ //if flags are different - CHANGE IN STRUCTURE
/*		        fprintf(output,"%d%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			currentmatrix[crow][5], //current structure
			nextmatrix[nrow][5], //structure in next snapshot
			currentmatrix[crow][0],// current snap
			nextmatrix[nrow][0], //next snap
			currentmatrix[crow][1], //O index
		        currentmatrix[crow][2], //CBH1 
			currentmatrix[crow][3], //CBH2
		        currentmatrix[crow][4], //CBH3
		        nextmatrix[nrow][1], //O index
			nextmatrix[nrow][2], //CBH1
			nextmatrix[nrow][3], //CBH2
			nextmatrix[nrow][4], //CBH4
			currentmatrix[crow][6], //shared proton for Zundels
			nextmatrix[nrow][6]); //shared proton for Zundels
*/		    }
		 }
	      }
	   }//end nrow-loop
	}//end crow-loop

//********************************************************
//Section IV: Identifying Eigen to Eigen PTs
//
//********************************************************
      for(crow=0; crow<nlines; crow++){
         for(coln=2; coln<5; coln++){
	    for(nrow=1; nrow<nlines; nrow++){
	       for(ncoln=2; ncoln<5; ncoln++){
	          if(nextmatrix[crow][coln]!=0){ //if matrix element is nonzero
		     if(currentmatrix[crow][coln]==nextmatrix[nrow][ncoln]){ //if H atoms matches
			if(currentmatrix[crow][1]!=nextmatrix[nrow][1]){ //If O atoms do NOT match
			   if(currentmatrix[crow][5]==1 && nextmatrix[nrow][5]==1){ //if Eigen-Eigen flag identified
			      fprintf(output2,"%d%d %d %d %d %d %d %f %f\n",
			      currentmatrix[crow][5], //initial flag
			      nextmatrix[nrow][5], //final flag
			      currentmatrix[crow][0], //initial snap
			      nextmatrix[nrow][0], //final snap
			      currentmatrix[crow][1], //initial O
			      nextmatrix[nrow][1], //final O
			      currentmatrix[crow][coln], //shared proton
			      currentdist[crow][coln], //dist. between initial O and H+
			      nextdist[nrow][ncoln]); //dist. between final O and H+
			   }
			}
		     }
		  }
	       }//end ncoln-loop
            }//end nrow-loop
         }//end coln-loop
      }//end crow-loop  
 
//*******************************************************
//Section V: Copy nextmatrix into currentmatrix 
//for proceeding set of comparisons.
//
//*******************************************************
//Intialize currentmatrix
         for(y=0; y<clines; y++){
	    for(z=0; z<8; z++){
	       currentmatrix[y][z]=0;
	    }//end z-loop
	 }//end y-loop

         for(y=0; y<300; y++){
            for(z=0; z<6; z++){
               currentdist[y][z]=0.0;
            }//end y-loop
         }//end z-loop

         clines=nlines;
         for(nrow=0; nrow<nlines; nrow++){
	    for(ncoln=0; ncoln<8; ncoln++){
	       currentmatrix[nrow][ncoln]=nextmatrix[nrow][ncoln];
	    }
	 }

	 for(nrow=0; nrow<nlines; nrow++){
	    for(ncoln=2; ncoln<5; ncoln++){
	       currentdist[nrow][ncoln]=nextdist[nrow][ncoln];
	    }
	 }

//Intialize nexetmatrix
         for(y=0; y<clines; y++){
	    for(z=0; z<8; z++){
	       nextmatrix[y][z]=0;
	    }//end z-loop
	 }//end y-loop

         for(y=0; y<300; y++){
            for(z=0; z<6; z++){
               nextdist[y][z]=0.0;
            }//end y-loop
         }//end z-loop

         nlines=0;
      }//end else-condition
   }// end of snap-loop 
   fclose(output);
   fclose(output2);
   fclose(output3);
}//end main-loop


