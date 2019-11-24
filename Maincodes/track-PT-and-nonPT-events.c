#include<stdio.h>
#include<stdlib.h>

//*********************************************************************
//Program created to track a specific proton in time by comparing the 
//change in structure type in current snap to the next.
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
//Last modified: January 22, 2018
//
//*********************************************************************

int main(){

//Parameters to read O-H GraphGeods
   int Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5;
   int snap,maxsnap,countOindex,x,y,z;
   int Oindex,Hhold1,Hhold2;
   int clines,crow,coln;
   int nlines,nrow,ncoln,wcoln;
   int wlines,wrow,xlines,xrow,xcoln;
   float dist,angle;
   FILE *input,*nextinput,*output,*output2;
   char ifile[100],nfile[100];

//Parameters for matrices
   int currentmatrix[300][7];
   int nextmatrix[300][7];
   int watermatrix[300][5][4];
   int nextwatermatrix[300][5][4];
 
//Intialize matrices
         for(y=0; y<300; y++){
            for(z=0; z<8; z++){
               currentmatrix[y][z]=0;
               nextmatrix[y][z]=0;
            }//end z-loop
         }//end y-loop

	 for(x=0; x<300; x++){
	    for(y=0; y<5; y++){
	       for(z=0; z<4; z++){
	          watermatrix[x][y][z]=0;
	       }//end z-loop
	    }//end y-loop
	 }//end z-loop
         clines=0;
	 wlines=0;
         countOindex=1;
  
   maxsnap=10; //CHANGE ACCORDINGLY 

   output2=fopen("undynamic-events-3.txt","a");
//   output=fopen("1-event-list.txt","w");
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
//Change flag type to 2
		         currentmatrix[crow][5]=2;
		         currentmatrix[nrow][5]=2;
			 currentmatrix[crow][6]=currentmatrix[crow][coln]; //store shared proton
			 currentmatrix[nrow][6]=currentmatrix[nrow][ncoln]; //store shared proton

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
               }
	    }//end coln-loop 
	 }//end nrow-loop
      }//end crow-loop
 
//*******************************************************
//Section I(c): Section to store only waters into 
//an array, startwaterlist.
//
//*******************************************************

      for(crow=0; crow<clines; crow++){
	 if(currentmatrix[crow][5]==55){ //If water is identified in currentmatrix
	    watermatrix[wlines][0][0]=currentmatrix[crow][0]; //snap
	    watermatrix[wlines][1][0]=currentmatrix[crow][1]; //Oxy
	    watermatrix[wlines][2][0]=currentmatrix[crow][2]; //H1
	    watermatrix[wlines][3][0]=currentmatrix[crow][3]; //H2
	    wlines++;
	 }
       }//end srow-loop

/*
       for(z=0; z<wlines; z++){
	  if(watermatrix[z][0]!=0){
	     printf("%d %d %d %d %d\n",
	     watermatrix[z][0][0],
	     watermatrix[z][1][0],
	     watermatrix[z][2][0],
	     watermatrix[z][3][0],
	     watermatrix[z][4][0]);
	  }
        }
   
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
	   xlines=0;

           while(fscanf(nextinput,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
              if(countOindex==1){
                 Oindex=Oxy;
                 Hhold1=Hyd;
              }
              else if(countOindex==2){
                 if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;

    	            nextmatrix[nlines][0]=snap+1;
		    nextmatrix[nlines][1]=Oxy;
		    nextmatrix[nlines][2]=Hhold1;
		    nextmatrix[nlines][3]=Hhold2;
		    nextmatrix[nlines][5]=55; //flagged as water
		    nlines++;

                 }
                 else{ //Should never trigger unless OH- exist       
	            countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
                 }
              }
              else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
                 if(Oindex==Oxy){
                    nextmatrix[nlines][0]=snap+1;
                    nextmatrix[nlines][1]=Oxy;
                    nextmatrix[nlines][2]=Hhold1;
                    nextmatrix[nlines][3]=Hhold2;
                    nextmatrix[nlines][4]=Hyd;
		    nextmatrix[nlines][5]=1; //flagged as eigen or H3O+

//Remove previous data on O atom, when H3O+ was no identified
                    nextmatrix[nlines-1][0]=0;  
                    nextmatrix[nlines-1][1]=0;
                    nextmatrix[nlines-1][2]=0;
                    nextmatrix[nlines-1][3]=0;
                    nextmatrix[nlines-1][4]=0;
		    nlines++;		  

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
         fclose(nextinput);
/*
  	for(z=0; z<nlines; z++){
	   if(nextmatrix[z][0]!=0){
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
//Change flag type to 2
		         nextmatrix[crow][5]=2;
		         nextmatrix[nrow][5]=2;
			 nextmatrix[crow][6]=nextmatrix[crow][coln];
			 nextmatrix[nrow][6]=nextmatrix[nrow][ncoln]; 

/*		         printf("%d %d %d %d %d %d\n",
	      	         nextmatrix[crow][0],
	      	         nextmatrix[crow][1],
	      	         nextmatrix[crow][2],
	      	         nextmatrix[crow][3],
	      	         nextmatrix[crow][4],
	      	         nextmatrix[crow][5]);
*/                       }
		     }
		  }
               }
	    }//end coln-loop 
	 }//end nrow-loop
      }//end crow-loop

 
//*******************************************************
//Section II(c): Section to store only waters into 
//an array, nextwaterlist.
//
//*******************************************************

       for(nrow=0; nrow<nlines; nrow++){
	 if(nextmatrix[nrow][5]==55){ //If water is identified in nextwatermatrix
	    nextwatermatrix[xlines][0][0]=nextmatrix[nrow][0]; //snap
	    nextwatermatrix[xlines][1][0]=nextmatrix[nrow][1]; //Oxy
	    nextwatermatrix[xlines][2][0]=nextmatrix[nrow][2]; //H1
	    nextwatermatrix[xlines][3][0]=nextmatrix[nrow][3]; //H2
	    xlines++;
	 }
       }//end nrow-loop



/*
       for(z=0; z<xlines; z++){
	  if(nextwatermatrix[z][0]!=0){
	     printf("%d %d %d %d\n",
	     nextwatermatrix[z][0][0],
	     nextwatermatrix[z][1][0],
	     nextwatermatrix[z][2][0],
	     nextwatermatrix[z][3][0]);
	  }
        }
*/
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
		        fprintf(output,"%d%d %d %d %d %d %d %d %d %d %d %d %d\n",
			currentmatrix[crow][5], //current structure
			nextmatrix[nrow][5], //structure in next snapshot
			currentmatrix[crow][0],// current snap
			nextmatrix[nrow][0], //next snap
			currentmatrix[crow][1], //O index
		        currentmatrix[crow][2], //CBH1 
			currentmatrix[crow][3], //CBH2
		        currentmatrix[crow][4], //CBH3
			nextmatrix[nrow][2], //CBH1
			nextmatrix[nrow][3], //CBH2
			nextmatrix[nrow][4], //CBH3
			currentmatrix[crow][6], //shared proton for Zundels
			nextmatrix[nrow][6]); //shared proton for Zundels
		    }
		 }
	      }
	   }//end nrow-loop
	}//end crow-loop

//*******************************************************
//Section IV: Identifying waters that do NOT participate
//in any PT events.  
//
//*******************************************************
/*      for(z=0; z<xlines; z++){
	  if(nextwatermatrix[z][0][0]!=0){
	     printf("next: %d %d %d %d\n",
	     nextwatermatrix[z][0][0],
	     nextwatermatrix[z][1][0],
	     nextwatermatrix[z][2][0],
	     nextwatermatrix[z][3][0]);
	  }
        }


       for(z=0; z<wlines; z++){
	  if(watermatrix[z][0][0]!=0){
	     printf("water: %d %d %d %d %d\n",
	     watermatrix[z][0][0],
	     watermatrix[z][1][0],
	     watermatrix[z][2][0],
	     watermatrix[z][3][0],
	     watermatrix[z][4][0]);
	  }
        }
*/
	for(wrow=0; wrow<wlines; wrow++){
           for(wcoln=2; wcoln<4; wcoln++){
   	      for(xrow=0; xrow<xlines; xrow++){
	         for(xcoln=2; xcoln<4; xcoln++){
	            if(watermatrix[wrow][0][0]!=0){
	               if(watermatrix[wrow][1][0]==nextwatermatrix[xrow][1][0]){ //if O index matches
			  if(watermatrix[wrow][wcoln][0]==nextwatermatrix[xrow][xcoln][0]){ //If any of the H atoms matches

		 	     if(watermatrix[wrow][wcoln][wcoln]==0 && nextwatermatrix[xrow][xcoln][xcoln]==0){ //if unflagged
				watermatrix[wrow][wcoln][wcoln]=1;
				nextwatermatrix[xrow][xcoln][xcoln]=1;
			     }
	                  }
		       }
	            }
		 }//end xcoln-loop 
	      }//end xrow-loop
           }//end wcoln-loop
	}//end wrow-loop
/*
      for(z=0; z<xlines; z++){
	  if(nextwatermatrix[z][0][0]!=0){
	     printf("next: %d %d %d %d %d %d\n",
	     nextwatermatrix[z][0][0],
	     nextwatermatrix[z][1][0],
	     nextwatermatrix[z][2][0],
	     nextwatermatrix[z][3][0],
	     nextwatermatrix[z][2][2],
	     nextwatermatrix[z][3][3]);
	  }
        }

       for(z=0; z<wlines; z++){
	  if(watermatrix[z][0][0]!=0){
	     printf("water: %d %d %d %d %d %d %d\n",
	     watermatrix[z][0][0],
	     watermatrix[z][1][0],
	     watermatrix[z][2][0],
	     watermatrix[z][3][0],
	     watermatrix[z][4][0],
	     watermatrix[z][2][2],
	     watermatrix[z][3][3]);
	  }
        }

*/

	for(wrow=0; wrow<wlines; wrow++){
	   for(wcoln=2; wcoln<4; wcoln++){
	      for(xrow=0; xrow<xlines; xrow++){
		 for(xcoln=2; xcoln<4; xcoln++){   
		    if(watermatrix[wrow][1][0]==nextwatermatrix[xrow][1][0]){ //if O index matches
		       if(watermatrix[wrow][wcoln][2]==1 && nextwatermatrix[xrow][xcoln][2]==1){ //If H is flagged
			  if(watermatrix[wrow][wcoln+1][3]==1 && nextwatermatrix[xrow][xcoln+1][3]==1){ //if BOTH H's are flagged
    		 	     if(nextwatermatrix[xrow][0][0]!=watermatrix[wrow][0][0]){ //if water remains unchanged in next snap
		                 watermatrix[wrow][4][0]=nextwatermatrix[xrow][0][0]; //update final time in watermatrix
			     }
			  }
       		       }
		       else if(watermatrix[wrow][wcoln][3]==1 && nextwatermatrix[xrow][xcoln][3]==1){ //If H is flagged
			  if(watermatrix[wrow][wcoln-1][2]==1 && nextwatermatrix[xrow][xcoln-1][2]==1){ //if BOTH H's are flagged
    		 	     if(nextwatermatrix[xrow][0][0]!=watermatrix[wrow][0][0]){ //if water remains unchanged in next snap
		                 watermatrix[wrow][4][0]=nextwatermatrix[xrow][0][0]; //update final time in watermatrix
			     }
			  }
       		       }
		    }
	         }//end xcoln-loop
	      }//end xrow-loop
	   }//end wcoln-loop
	}//end wrow-loop


/*
        if(snap==maxsnap-1){
           for(wrow=0; wrow<wlines; wrow++){
              if(watermatrix[wrow][0][0]!=0){
                 printf("%d %d %d %d %d\n",
                 watermatrix[wrow][0][0], //snap 
                 watermatrix[wrow][4][0], //final snap
                 watermatrix[wrow][1][0], //O index
                 watermatrix[wrow][2][0], //H1
                 watermatrix[wrow][3][0]); //H2
              }
           }//end wrow-loop
        }
*/
   
//*******************************************************
//Section V: Copy nextmatrix into currentmatrix and
//nextwatermatrix into startwatermatrix for 
//proceeding set of comparisons.
//
//*******************************************************
//Intialize currentmatrix
         for(y=0; y<clines; y++){
	    for(z=0; z<8; z++){
	       currentmatrix[y][z]=0;
	    }//end z-loop
	 }//end y-loop
       
         clines=nlines;
         for(nrow=0; nrow<nlines; nrow++){
	    for(ncoln=0; ncoln<8; ncoln++){
	       currentmatrix[nrow][ncoln]=nextmatrix[nrow][ncoln];
	    }
	 }

//Intialize nextmatrix
         for(y=0; y<clines; y++){
	    for(z=0; z<8; z++){
	       nextmatrix[y][z]=0;
	    }//end z-loop
	 }//end y-loop

	 for(y=0; y<xlines; y++){
	    for(z=0; z<5; z++){
	       nextwatermatrix[y][z][0]=0;
	       nextwatermatrix[y][z][0]=0;
	       nextwatermatrix[y][z][0]=0;
	       nextwatermatrix[y][z][0]=0;
	       nextwatermatrix[y][2][2]=0;
	       nextwatermatrix[y][3][3]=0;
	    }
	 }

	 for(y=0; y<300; y++){
	    watermatrix[y][2][2]=0; //initialize flags
	    watermatrix[y][3][3]=0;
	 }
	 xlines=0;
         nlines=0;
      }//end else-condition
   }// end of snap-loop 
  // fclose(output);
   fclose(output2);
}//end main-loop


