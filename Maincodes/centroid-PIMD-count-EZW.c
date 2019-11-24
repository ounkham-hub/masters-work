#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//*************************************************************************
//Program created to count the number of Eigens,Zundels, and waters
//per snapshot and across all snapshots. 
//
//There are sections of the program that account for the following
//artificial Zundel-like structures.
//
//1) Eigen-water, not counted as Zundel contributes single proton
//2) Double Zundel, occurs when two H's of an H3O+ form Zundels with other
//waters, still only contributes 1 proton
//
//Created by Lelee Ounkham
//Last modified on August 25, 2018
//
//*************************************************************************


int main(){

//Parameters to read O-H GraphGeod

   int Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5;
   int snap,maxsnap,countOindex,x,y,z;
   int Oindex,Hhold1,Hhold2,clines;
   int crow,coln,nrow,ncoln,findproton,waterO;
   int counteigen,countzundel,countwater,countproton,numofprotons;
   float dist,angle,numofzundels;
   FILE *input,*output;
   char ifile[100];

//Parameters for matrices
   int currentmatrix[300][7];

//Initialize matrices
   for(y=0; y<300; y++){
     for(z=0; z<8; z++){
	currentmatrix[y][z]=0;
     }
   }
   clines=0;
   counteigen=0;
   countzundel=0;
   countwater=0;
   countOindex=1;
   countproton=0;

//Adjust maxsnap accordingly
   maxsnap=61555;
  
   output=fopen("PIMD-complete-count-of-EZW.txt","a"); 
   for(snap=1; snap<=maxsnap; snap++){
         sprintf(ifile,"R-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
         input=fopen(ifile,"r");

      while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
         if(countOindex==1){
             Oindex=Oxy;
             Hhold1=Hyd;
          }
          else if(countOindex==2){
             if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;

                    currentmatrix[clines][0]=snap;
                    currentmatrix[clines][1]=Oxy;
                    currentmatrix[clines][2]=Hhold1;
                    currentmatrix[clines][3]=Hhold2;
                    currentmatrix[clines][5]=55; //flagged as water
                    clines++;
                 }
                 else{ //Should never trigger unless OH- exist       
                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
                 }
              }
              else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
                 if(Oindex==Oxy){
                    currentmatrix[clines][0]=snap;
                    currentmatrix[clines][1]=Oxy;
                    currentmatrix[clines][2]=Hhold1;
                    currentmatrix[clines][3]=Hhold2;
                    currentmatrix[clines][4]=Hyd;
                    currentmatrix[clines][5]=1; //flagged as eigen or H3O+

//Remove previous data on O atom, when H3O+ was no identified
                    currentmatrix[clines-1][0]=0;
                    currentmatrix[clines-1][1]=0;
                    currentmatrix[clines-1][2]=0;
                    currentmatrix[clines-1][3]=0;
                    currentmatrix[clines-1][4]=0;
                    clines++;

                    countOindex=0; //restarts Oindex to read the next line;
                 }
                 else{ //Water molecule identified and countOindex==2
                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
                 }
              }
              countOindex++;
         }//end while-loop
         fclose(input);
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

//*************************************************************************
//Section I(b): Identify which Eigens are a part of a Zundel - when two 
//Eigens share a proton.
//
//*************************************************************************

         for(crow=0; crow<clines; crow++){
            for(coln=2; coln<5; coln++){
               for(nrow=1; nrow<clines; nrow++){
                  for(ncoln=2; ncoln<5; ncoln++){
                     if(currentmatrix[crow][coln]!=0){ //if matrix element is nonzero
                        if(currentmatrix[crow][coln]==currentmatrix[nrow][ncoln]){ //If any H's atoms match
                           if(currentmatrix[crow][1]!=currentmatrix[nrow][1]){ //and if O' atoms differ - Zundel identified
//Change flag type to 2
                              currentmatrix[crow][5]=2;
                              currentmatrix[nrow][5]=2;
                              currentmatrix[crow][6]=currentmatrix[crow][coln]; //store shared proton
                              currentmatrix[nrow][6]=currentmatrix[nrow][ncoln]; //store shared proton

/*                              printf("%d %d %d %d %d %d\n",
                              currentmatrix[crow][0],
                              currentmatrix[crow][1],
                              currentmatrix[crow][2],
                              currentmatrix[crow][3],
                              currentmatrix[crow][4],
                              currentmatrix[crow][5]);
*/                          }
                        }
                     }
                  }
               }//end coln-loop 
            }//end nrow-loop
         }//end crow-loop


//*************************************************************************
//Section I(c): Loop through all Zundels and remove any that artifacts
//where there's only a water and an H3O+.
//
//*************************************************************************

        for(crow=0; crow<clines; crow++){
	   for(coln=2; coln<5; coln++){
	      if(currentmatrix[crow][5]==2){ //if Zundel is identified
		 if(currentmatrix[crow][coln]==0){ //if any of the H indices is empty - NOT a Zundel
/*                    printf("%d %d %d %d %d %d %d\n",
                    currentmatrix[crow][0],
                    currentmatrix[crow][1],
                    currentmatrix[crow][2],
                    currentmatrix[crow][3],
                    currentmatrix[crow][4],
                    currentmatrix[crow][5],
		    currentmatrix[crow][6]);
*/
		    currentmatrix[crow][5]=55; //these are actually waters
		    findproton=currentmatrix[crow][6];//store proton index
		    waterO=currentmatrix[crow][1]; //store O index

		    for(y=0; y<clines; y++){
		       if(currentmatrix[y][5]==2){//if Zundel identified
			  if(currentmatrix[y][6]==findproton){//if h index found
			     if(currentmatrix[y][1]!=waterO){ //if O indices are NOT the same, partner found
/*	                        printf("%d %d %d %d %d %d %d\n",
                                currentmatrix[y][0],
                                currentmatrix[y][1],
                                currentmatrix[y][2],
                                currentmatrix[y][3],
                                currentmatrix[y][4],
                                currentmatrix[y][5],
		                currentmatrix[y][6]);
*/
				currentmatrix[y][5]=1; //change to Eigen
			     }
			  }
		       }//end y-loop
		    }//end y-loop			        	
		 }
	      }
	   }//end coln-loop
	}//end crow-loop

	for(crow=0; crow<clines; crow++){
	   for(coln=2; coln<5; coln++){
	      if(currentmatrix[crow][5]==1){
	         if(currentmatrix[crow][coln]==0){ //if any of the H indeices is empty - not an Eigen but a water
		    currentmatrix[crow][5]=55;
		 }
	      }
	   }
	}

//*************************************************************************
//Section II: Count the number of unique protons (H3O+), how many
//of them are involved in an Eigen or Zundel.
//
//*************************************************************************

	for(crow=0; crow<clines; crow++){
	   if(currentmatrix[crow][0]!=0){
	      if(currentmatrix[crow][5]==1){
	         counteigen++;
	         countproton++;
	      }
	      else if(currentmatrix[crow][5]==2){
	         countzundel++;
	      }
	      else if(currentmatrix[crow][5]==55){
	         countwater++;
	      }
	   }
	}

//To account for double Zundels - it should count as 1 proton and 1 water, not 2 protons
	if(countzundel==1){
	   fprintf(output,"%d %d %d %d %d\n",
   	   snap,
	   countproton+(countzundel/2),
	   counteigen,
	   (countzundel/2),
	   countwater+1);
	}
	else if(countzundel==3){
	   fprintf(output,"%d %d %d %d %d\n",
   	   snap,
	   countproton+(countzundel/2),
	   counteigen,
	   (countzundel/2),
	   countwater+1);
	}
	else if(countzundel==5){
	   fprintf(output,"%d %d %d %d %d\n",
   	   snap,
	   countproton+(countzundel/2),
	   counteigen,
	   (countzundel/2),
	   countwater+1);
	}
	else if(countzundel==7){
	   fprintf(output,"%d %d %d %d %d\n",
   	   snap,
	   countproton+(countzundel/2),
	   counteigen,
	   (countzundel/2),
	   countwater+1);
	}
	else if(countzundel!=1 && countzundel!=3 && countzundel!=5 && countzundel!=7){
	   fprintf(output,"%d %d %d %d %d\n",
   	   snap,
	   countproton+(countzundel/2),
	   counteigen,
	   (countzundel/2),
	   countwater);
	}

/*
	   fprintf(output,"%d %d %d %d %d\n",
   	   snap,
	   countproton+(countzundel/2),
	   counteigen,
	   (countzundel/2),
	   countwater);
*/
/*
	for(crow=0; crow<clines; crow++){
if(snap==8161){
	   if(currentmatrix[crow][0]!=0){
              printf("%d %d %d %d %d %d\n",
              currentmatrix[crow][0],
              currentmatrix[crow][1],
              currentmatrix[crow][2],
              currentmatrix[crow][3],
              currentmatrix[crow][4],
              currentmatrix[crow][5]);
}
	   }
	}
*/

//*************************************************************************
//Section III: Initialize counters and matrices for next snapshot
//*************************************************************************

         for(y=0; y<clines; y++){
            for(z=0; z<8; z++){
               currentmatrix[y][z]=0;
            }//end z-loop
         }//end y-loop

	 clines=0;
	 countOindex=1;

	 countwater=0;
	 counteigen=0;
	 countzundel=0;
	 countproton=0;
     }//end snap-loop
     fclose(output);
} //end main-loop
