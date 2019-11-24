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
//**Include calculating the duration of the waters
//
//Input: R-HCl.input.solB#.xyz.solC#.xyz.GraphGeod
//Output: event-list.txt
//
//Created by Lelee Ounkham
//Last modified: April 20, 2018
//
//*********************************************************************

int main(){

//Parameters to read O-H GraphGeods
   int Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5;
   int snap,maxsnap,countOindex,x,y,z;
   int Oindex,Hhold1,Hhold2;
   int clines,crow,coln,diff;
   int nlines,nrow,ncoln;
   int mrow;
   int numH,index;
   float dist,angle;
   FILE *input,*nextinput,*output;
   char ifile[100],nfile[100];

   numH=213; //add 1 since matrix starts at 1 and not 0
//Parameters for matrices
   int currentmatrix[300][7];
   int nextmatrix[300][7];
   int mainmatrix[numH][5];
 
//Intialize matrices
         for(y=0; y<300; y++){
            for(z=0; z<8; z++){
               currentmatrix[y][z]=0;
               nextmatrix[y][z]=0;
            }//end z-loop
         }//end y-loop
         clines=0;
         countOindex=1;
	 index=102;

//Number mainmatrix

   for(y=1; y<numH; y++){
      index++;
      mainmatrix[y][0]=index; 
//      printf("%d %d %d\n",mainmatrix[y][0],mainmatrix[y][1],mainmatrix[y][2]);
      for(z=1; z<5; z++){	
         mainmatrix[y][z]=0;
      }
   } 
  
   maxsnap=61555; //CHANGE ACCORDINGLY 

   output=fopen("undynamic-events-diff.txt","a");
   for(snap=1; snap<maxsnap; snap++){
      if(snap==maxsnap){ //END OF PROGRAM
//Print mainmatrix, may contain O atoms that have not participated in PT
/*         for(y=1; y<numH; y++){
            printf("%d %d %d %d %d\n",
	    mainmatrix[y][0],
	    mainmatrix[y][1],
	    mainmatrix[y][2],
	    mainmatrix[y][3],
	    mainmatrix[y][4]);
	 }//end y-loop
*/
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

/*		         printf("%d %d %d %d %d %d\n",
	      	         currentmatrix[crow][0],
	      	         currentmatrix[crow][1],
	      	         currentmatrix[crow][2],
	      	         currentmatrix[crow][3],
	      	         currentmatrix[crow][4],
	      	         currentmatrix[crow][5]);
*/                      }
		     }
		  }
               }
	    }//end coln-loop 
	 }//end nrow-loop
      }//end crow-loop

//*******************************************************
//Section I(c): Section identify which O's and H's are 
//associated with water,Zundel, or Eigen.
//
//*******************************************************

      for(mrow=1; mrow<numH; mrow++){
         for(crow=0; crow<clines; crow++){
            for(coln=2; coln<5; coln++){
	       if(mainmatrix[mrow][0]==currentmatrix[crow][coln]){ //if H index matches
	          mainmatrix[mrow][1]=currentmatrix[crow][1]; //O index
	          mainmatrix[mrow][2]=currentmatrix[crow][5]; //flag
	          mainmatrix[mrow][3]=currentmatrix[crow][0]; //initial snap
	          mainmatrix[mrow][4]=currentmatrix[crow][0]; //final snap
/*
                  printf("%d %d %d %d %d\n",
	          mainmatrix[mrow][0],
	          mainmatrix[mrow][1],
	          mainmatrix[mrow][2],
	          mainmatrix[mrow][3],
	          mainmatrix[mrow][4]);
*/	       }
	    }
	 }
      }
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
	    for(nrow=0; nrow<nlines; nrow++){
	       for(ncoln=2; ncoln<5; ncoln++){
	          if(nextmatrix[crow][coln]!=0){ //if matrix element is nonzero
                     if(nextmatrix[crow][coln]==nextmatrix[nrow][ncoln]){ //If any H's atoms match
		        if(nextmatrix[crow][1]!=nextmatrix[nrow][1]){ //and if O' atoms differ - Zundel identified
//Change flag type to 2
		         nextmatrix[crow][5]=2;
		         nextmatrix[nrow][5]=2;
			 nextmatrix[crow][6]=nextmatrix[crow][coln];
			 nextmatrix[nrow][6]=nextmatrix[nrow][ncoln]; 
/*
		         printf("%d %d %d %d %d %d\n",
	      	         nextmatrix[crow][0],
	      	         nextmatrix[crow][1],
	      	         nextmatrix[crow][2],
	      	         nextmatrix[crow][3],
	      	         nextmatrix[crow][4],
	      	         nextmatrix[crow][5]);
*/                     }
		     }
		  }
               }
	    }//end coln-loop 
	 }//end nrow-loop
      }//end crow-loop

//*******************************************************
//Section III: Identifying waters that do NOT participate
//in any PT events.  
//
//*******************************************************

   for(mrow=1; mrow<numH; mrow++){
      for(nrow=0; nrow<nlines; nrow++){
         for(ncoln=2; ncoln<5; ncoln++){
            if(mainmatrix[mrow][0]==nextmatrix[nrow][ncoln]){ //If H indices matches
               if(mainmatrix[mrow][2]==nextmatrix[nrow][5]){ //if flags matches
	          if(mainmatrix[mrow][1]==nextmatrix[nrow][1]){ //If O indices matches
	             mainmatrix[mrow][4]=mainmatrix[mrow][4]+1; //update final snapshot 
                  }
               }
	       else if(mainmatrix[mrow][2]==55 && nextmatrix[nrow][5]!=55){ //change in structure!
//Calculate the lifetimes

	          diff=mainmatrix[mrow][4]-mainmatrix[mrow][3];
//Print initial and final snap. water appeared
	          fprintf(output,"%d %d %d %d %d %d\n",
		  mainmatrix[mrow][0],
		  mainmatrix[mrow][1],
		  mainmatrix[mrow][2],
		  mainmatrix[mrow][3],
		  mainmatrix[mrow][4],
		  diff+1);
//Update structure and time info for O index
		  mainmatrix[mrow][2]=nextmatrix[nrow][5]; //update flag
		  mainmatrix[mrow][1]=nextmatrix[nrow][1]; //update O index
		  mainmatrix[mrow][3]=nextmatrix[nrow][0]; //initial snap
		  mainmatrix[mrow][4]=nextmatrix[nrow][0]; //restart final snap 
	       }
	       else if(mainmatrix[mrow][2]!=55 && nextmatrix[nrow][5]==55){ //If zundel or eigen becomes water
		  mainmatrix[mrow][2]=nextmatrix[nrow][5]; //update flag
		  mainmatrix[mrow][1]=nextmatrix[nrow][1]; ////update O index
		  mainmatrix[mrow][3]=nextmatrix[nrow][0]; //update initial snap
		  mainmatrix[mrow][4]=nextmatrix[nrow][0]; //restart final snap
	       }
	       else if(mainmatrix[mrow][4]!=55 && nextmatrix[nrow][5]!=55){ //If Eigen or Zundel transforms to Zundel or Eigen
                  mainmatrix[mrow][2]=nextmatrix[nrow][5]; //update flag
		  mainmatrix[mrow][1]=nextmatrix[nrow][1]; //update O index
		  mainmatrix[mrow][3]=nextmatrix[nrow][0]; //update initial snap
		  mainmatrix[mrow][4]=nextmatrix[nrow][0]; //restart final snap
	       }
	    }
	 }//end ncoln-loop
      }//end nrow-loop
   }//end mrow-loop 
/*
         for(y=1; y<numH; y++){
            printf("%d %d %d %d %d\n",
	    mainmatrix[y][0],
	    mainmatrix[y][1],
	    mainmatrix[y][2],
	    mainmatrix[y][3],
	    mainmatrix[y][4]);
	 }//end y-loop
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
         nlines=0;

      }//end else-condition
   }// end of snap-loop 
   fclose(output);
//Print mainmatrix, may contain O atoms that have not participated in PT
/*         for(y=1; y<numH; y++){
            printf("%d %d %d %d %d\n",
	    mainmatrix[y][0],
	    mainmatrix[y][1],
	    mainmatrix[y][2],
	    mainmatrix[y][3],
	    mainmatrix[y][4]);
	 }//end y-loop
*/
}//end main-loop
