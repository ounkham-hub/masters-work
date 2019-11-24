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
   int Oxy,Hyd;
   int snap,maxsnap,countOindex,x,y,z;
   int Oindex,Hhold1,Hhold2;
   int clines,crow,coln,diff;
   int nlines,nrow,ncoln;
   int mrow;
   int numH,index;
   float dist,weight,whold1,whold2;
   FILE *input,*nextinput,*output,*output2,*output3;
   char ifile[100],nfile[100],ofile2[100],ofile3[100];

   numH=213; //add 1 since matrix starts at 1 and not 0
//Parameters for matrices
   int currentmatrix[300][7];
   float currentweight[300][5];
   int nextmatrix[300][7];
   float nextweight[300][7];
   int mainmatrix[numH][5];
   float mainweight[numH][5];

 
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

   output=fopen("weighted-unreactive-zundels.txt","a");
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
         sprintf(ifile,"HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap,snap); //open current snap
//         sprintf(ofile2,"OCO-list-%d.txt",snap); //open current snap
         input=fopen(ifile,"r");
//         output2=fopen(ofile2,"w");
//*******************************************************
//Section I(a): Identify which O's are a part of a 
//water OR Eigen structure for current snapshot.
//
//*******************************************************
         while(fscanf(input,"%d %d %f %f\n",&Oxy,&Hyd,&dist,&weight)==4){
              if(countOindex==1){
                 Oindex=Oxy;
                 Hhold1=Hyd;
		 whold1=weight;
              }
              else if(countOindex==2){
                 if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;
		    whold2=weight;

    	            currentmatrix[clines][0]=snap;
		    currentmatrix[clines][1]=Oxy;
		    currentmatrix[clines][2]=Hhold1;
		    currentmatrix[clines][3]=Hhold2;
		    currentmatrix[clines][5]=55; //flagged as water

		    currentweight[clines][2]=whold1;
		    currentweight[clines][3]=whold2;
		    clines++;
                 }
                 else{ //Should never trigger unless OH- exist       
	            countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
		    whold1=weight;
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

		    currentweight[clines][2]=whold1;
		    currentweight[clines][3]=whold2;
		    currentweight[clines][4]=weight;

//Remove previous data on O atom, when H3O+ was no identified
                    currentmatrix[clines-1][0]=0;  
                    currentmatrix[clines-1][1]=0;
                    currentmatrix[clines-1][2]=0;
                    currentmatrix[clines-1][3]=0;
                    currentmatrix[clines-1][4]=0;

		    currentweight[clines-1][2]=0.0;
		    currentweight[clines-1][3]=0.0;
		    currentweight[clines-1][4]=0.0;
		    clines++;		  

                    countOindex=0; //restarts Oindex to read the next line;
                 }
                 else{ //Water molecule identified and countOindex==2
                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
		    whold1=weight;
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

/*
	for(z=0; z<clines; z++){
	   if(currentmatrix[z][0]!=0){
	      printf("%d %d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	      currentmatrix[z][0],
	      currentmatrix[z][1],
	      currentmatrix[z][2],
	      currentmatrix[z][3],
	      currentmatrix[z][4],
	      currentmatrix[z][5],
	      currentmatrix[z][6],
	      currentweight[z][2],
	      currentweight[z][3],
	      currentweight[z][4]);
	   }
	}//end z-loop
*/
//      fclose(output2);

//*******************************************************
//Section I(c): Section to identify waters, Eigens, 
//Zundels applying predetermined weight O-H cutoff
//=> 0.700 for O-H in Zundels and Eigens and 1.000
//for waters.
//
//*******************************************************

     for(crow=0; crow<clines; crow++){
	for(coln=2; coln<5; coln++){
	   if(currentmatrix[crow][1]!=0){
	      if(currentmatrix[crow][5]==55){ //for all waters
	         if(currentweight[crow][coln]>0.000 && currentweight[crow][coln]<1.000){ //if O-H bond has a weight less than 1.000, remove
	     	    currentmatrix[crow][0]=0;
		    currentmatrix[crow][1]=0;
		    currentmatrix[crow][2]=0;
		    currentmatrix[crow][3]=0;
		    currentmatrix[crow][4]=0;
		    currentmatrix[crow][5]=0;

		    currentmatrix[crow][2]=0.0;
		    currentmatrix[crow][3]=0.0;
		    currentmatrix[crow][4]=0.0; 
	         }  
	      }
	      else if(currentmatrix[crow][5]==1){ //for all eigens
		 if(currentweight[crow][coln]>0.700){
//Leave as Eigen
		 }
		 else if(currentweight[crow][coln]<0.700){ //Structure is more water-like, remove
		    currentmatrix[crow][0]=0;
		    currentmatrix[crow][1]=0;
		    currentmatrix[crow][2]=0;
		    currentmatrix[crow][3]=0;
		    currentmatrix[crow][4]=0;
		    currentmatrix[crow][5]=0;

		    currentmatrix[crow][2]=0.0;
		    currentmatrix[crow][3]=0.0;
		    currentmatrix[crow][4]=0.0; 
		 }
	      }
	      else if(currentmatrix[crow][5]==2){ //for all zundels
		 if(currentweight[crow][coln]>0.700){
//Leave as Zundel
		 }
		 else if(currentweight[crow][coln]<=0.700){ //change to Eigen
		    currentmatrix[crow][5]=1;
//Loop through and change partner flag to Eigen too
		    for(y=0; y<clines; y++){
	 	       for(z=2; z<5; z++){
			  if(currentmatrix[y][5]==2){ //if Zundel
			     if(currentmatrix[crow][coln]==currentmatrix[y][z]){ //If H index matches
			        if(currentmatrix[crow][1]!=currentmatrix[y][1]){ //BUT O indices are different - Zundel partner identified
			           currentmatrix[y][5]=1; //change flag type to Eigen
			        }
			     }
		          }//end z-loop	
		       }
		    }//end y-loop
		 }
	      }
	   }
	}//end coln-loop
     }//end crow-loop

//Remove all Eigens that are more water-like and have an O-H bond <0.700
     for(crow=0; crow<clines; crow++){
	for(coln=2; coln<5; coln++){
	   if(currentmatrix[crow][1]!=0){
	      if(currentmatrix[crow][5]==1){ //for all waters
	         if(currentweight[crow][coln]<0.700){ //if O-H bond has a weight less than 1.000, remove
	     	    currentmatrix[crow][0]=0;
		    currentmatrix[crow][1]=0;
		    currentmatrix[crow][2]=0;
		    currentmatrix[crow][3]=0;
		    currentmatrix[crow][4]=0;
		    currentmatrix[crow][5]=0;

		    currentmatrix[crow][2]=0.0;
		    currentmatrix[crow][3]=0.0;
		    currentmatrix[crow][4]=0.0; 
	         }  
	      }
	   }
	}//end coln-loop
     }//end crow-loop
/*
	for(z=0; z<clines; z++){
	   if(currentmatrix[z][0]!=0){
	      printf("%d %d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	      currentmatrix[z][0],
	      currentmatrix[z][1],
	      currentmatrix[z][2],
	      currentmatrix[z][3],
	      currentmatrix[z][4],
	      currentmatrix[z][5],
	      currentmatrix[z][6],
	      currentweight[z][2],
	      currentweight[z][3],
	      currentweight[z][4]);
	   }
	}//end z-loop
*/
//*******************************************************
//Section I(d): Section identify which O's and H's are 
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
           sprintf(nfile,"HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap+1,snap+1); //open next snap
//           sprintf(ofile3,"OCO-list-%d.txt",snap+1); //open current snap
           nextinput=fopen(nfile,"r");
//           output3=fopen(ofile3,"w");

           countOindex==1;
           nlines=0;

           while(fscanf(nextinput,"%d %d %f %f\n",&Oxy,&Hyd,&dist,&weight)==4){
              if(countOindex==1){
                 Oindex=Oxy;
                 Hhold1=Hyd;	
		 whold1=weight;
              }
              else if(countOindex==2){
                 if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;
		    whold2=weight;

    	            nextmatrix[nlines][0]=snap+1;
		    nextmatrix[nlines][1]=Oxy;
		    nextmatrix[nlines][2]=Hhold1;
		    nextmatrix[nlines][3]=Hhold2;
		    nextmatrix[nlines][5]=55; //flagged as water

		    nextweight[nlines][2]=whold1;
		    nextweight[nlines][3]=whold2;
		    nlines++;

                 }
                 else{ //Should never trigger unless OH- exist       
	            countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
		    whold1=weight;
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

		    nextweight[nlines][2]=whold1;
		    nextweight[nlines][3]=whold2;
		    nextweight[nlines][4]=weight;

//Remove previous data on O atom, when H3O+ was no identified
                    nextmatrix[nlines-1][0]=0;  
                    nextmatrix[nlines-1][1]=0;
                    nextmatrix[nlines-1][2]=0;
                    nextmatrix[nlines-1][3]=0;
                    nextmatrix[nlines-1][4]=0;

		    nextweight[nlines-1][2]=0.0;
		    nextweight[nlines-1][3]=0.0;
		    nextweight[nlines-1][4]=0.0;
		    nlines++;		  

                    countOindex=0; //restarts Oindex to read the next line;
                 }
                 else{ //Water molecule identified and countOindex==2

                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
		    whold1=weight;
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

/*
  	for(z=0; z<nlines; z++){
	   if(nextmatrix[z][0]!=0){
	      printf("%d %d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	      nextmatrix[z][0],
	      nextmatrix[z][1],
	      nextmatrix[z][2],
	      nextmatrix[z][3],
	      nextmatrix[z][4],
	      nextmatrix[z][5],
	      nextmatrix[z][6],
	      nextweight[z][2],
	      nextweight[z][3],	
	      nextweight[z][4]);
	   }
	}//end z-loop
*/
//	fclose(output3);


//*******************************************************
//Section II(c): Section to identify waters, Eigens, 
//Zundels applying predetermined weight O-H cutoff
//=> 0.700 for O-H in Zundels and Eigens and 1.000
//for waters.
//
//*******************************************************

     for(nrow=0; nrow<nlines; nrow++){
	for(ncoln=2; ncoln<5; ncoln++){
	   if(nextmatrix[nrow][1]!=0){
	      if(nextmatrix[nrow][5]==55){ //for all waters
	         if(nextweight[nrow][ncoln]>0.000 && nextweight[nrow][ncoln]<1.000){ //if O-H bond has a weight less than 1.000, remove
	     	    nextmatrix[nrow][0]=0;
		    nextmatrix[nrow][1]=0;
		    nextmatrix[nrow][2]=0;
		    nextmatrix[nrow][3]=0;
		    nextmatrix[nrow][4]=0;
		    nextmatrix[nrow][5]=0;

		    nextmatrix[nrow][2]=0.0;
		    nextmatrix[nrow][3]=0.0;
		    nextmatrix[nrow][4]=0.0; 
	         }  
	      }
	      else if(nextmatrix[nrow][5]==1){ //for all eigens
		 if(nextweight[nrow][ncoln]>0.700){
//Leave as Eigen
		 }
		 else if(nextweight[nrow][ncoln]<0.700){ //Structure is more water-like, remove
		    nextmatrix[nrow][0]=0;
		    nextmatrix[nrow][1]=0;
		    nextmatrix[nrow][2]=0;
		    nextmatrix[nrow][3]=0;
		    nextmatrix[nrow][4]=0;
		    nextmatrix[nrow][5]=0;

		    nextmatrix[nrow][2]=0.0;
		    nextmatrix[nrow][3]=0.0;
		    nextmatrix[nrow][4]=0.0; 
		 }
	      }
	      else if(nextmatrix[nrow][5]==2){ //for all zundels
		 if(nextweight[nrow][ncoln]>0.700){
//Leave as Zundel
		 }
		 else if(nextweight[nrow][ncoln]<=0.700){ //change to Eigen
		    nextmatrix[nrow][5]=1;
//Loop through and change partner flag to Eigen too
		    for(y=0; y<nlines; y++){
	 	       for(z=2; z<5; z++){
			  if(nextmatrix[y][5]==2){ //if Zundel
			     if(nextmatrix[nrow][ncoln]==nextmatrix[y][z]){ //If H index matches
			        if(nextmatrix[nrow][1]!=nextmatrix[y][1]){ //BUT O indices are different - Zundel partner identified
			           nextmatrix[y][5]=1; //change flag type to Eigen
			        }
			     }
		          }//end z-loop	
		       }
		    }//end y-loop
		 }
	      }
	   }
	}//end ncoln-loop
     }//end nrow-loop

//Remove all Eigens that are more water-like and have an O-H bond <0.700
     for(nrow=0; nrow<nlines; nrow++){
	for(ncoln=2; ncoln<5; ncoln++){
	   if(nextmatrix[nrow][1]!=0){
	      if(nextmatrix[nrow][5]==1){ //for all waters
	         if(nextweight[nrow][ncoln]<0.700){ //if O-H bond has a weight less than 1.000, remove
	     	    nextmatrix[nrow][0]=0;
		    nextmatrix[nrow][1]=0;
		    nextmatrix[nrow][2]=0;
		    nextmatrix[nrow][3]=0;
		    nextmatrix[nrow][4]=0;
		    nextmatrix[nrow][5]=0;

		    nextmatrix[nrow][2]=0.0;
		    nextmatrix[nrow][3]=0.0;
		    nextmatrix[nrow][4]=0.0; 
	         }  
	      }
	   }
	}//end coln-loop
     }//end crow-loop
/*
	for(z=0; z<nlines; z++){
	   if(nextmatrix[z][0]!=0){
	      printf("%d %d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	      nextmatrix[z][0],
	      nextmatrix[z][1],
	      nextmatrix[z][2],
	      nextmatrix[z][3],
	      nextmatrix[z][4],
	      nextmatrix[z][5],
	      nextmatrix[z][6],
	      nextweight[z][2],
	      nextweight[z][3],
	      nextweight[z][4]);
	   }
	}//end z-loop
*/
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
	       else if(mainmatrix[mrow][2]==2 && nextmatrix[nrow][5]!=2){ //change in structure!
//Calculate the lifetimes

	          diff=mainmatrix[mrow][4]-mainmatrix[mrow][3];
if(diff>1){
//Print initial and final snap. water appeared
	          fprintf(output,"%d %d %d %d %d %d\n",
		  mainmatrix[mrow][0],
		  mainmatrix[mrow][1],
		  mainmatrix[mrow][2],
		  mainmatrix[mrow][3],
		  mainmatrix[mrow][4],
		  diff+1);
}
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
      
	 for(y=0; y<clines; y++){
	    for(z=0; z<5; z++){
	       currentweight[y][z]=0.0;
	    }
	 } 

         clines=nlines;
         for(nrow=0; nrow<nlines; nrow++){
	    for(ncoln=0; ncoln<8; ncoln++){
	       currentmatrix[nrow][ncoln]=nextmatrix[nrow][ncoln];
	    }
	 }

         for(nrow=0; nrow<nlines; nrow++){
	    for(ncoln=0; ncoln<5; ncoln++){
	       currentweight[nrow][ncoln]=nextweight[nrow][ncoln];
	    }
	 }

//Intialize nextmatrix
         for(y=0; y<clines; y++){
	    for(z=0; z<8; z++){
	       nextmatrix[y][z]=0;
	    }//end z-loop
	 }//end y-loop
      
	 for(y=0; y<clines; y++){
	    for(z=0; z<5; z++){
	       nextweight[y][z]=0.0;
	    }
	 } 

         nlines=0;

      }//end else-condition
   }// end of snap-loop 
   fclose(output);
}//end main-loop
