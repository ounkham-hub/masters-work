#include<stdio.h>
#include<stdlib.h>

//*****************************************************************************
//Program created to identify all waters,eigens, and zundels based on a weight 
//cutoff then output associated weight in each snapshot.
//
//Flags for each structure:
//Eigens 1, O-H > 0.900
//Zundels 2 O-H > 0.700
//Waters 55 O-H = 1.000
//
//Input: HCl.covalent.O#.xyz.H#.xyz.wGraphGeod
//Output: c-EZW-O#.H#.wGraphGeod
//
//Created by Lelee Ounkham
//Last Modified June 14, 2018
//
//******************************************************************************

int main(){

  int Oxy,Hyd;
  int snap,maxsnap,Oindex,Hhold1,Hhold2;
  int y,z,clines,countOindex;
  int crow,coln,nrow,ncoln;
  float dist,weight,whold1,whold2;
  FILE *input,*output;
  char ifile[100],ofile[100];

  int currentmatrix[300][10];
  float weightmatrix[300][7];

//Change variables accordingly 
   maxsnap=61555;
   clines=0; 
   countOindex=1;

//Intialize matrices
   for(y=0; y<300; y++){
      for(z=0; z<10; z++){
	 currentmatrix[y][z]=0;
      }
   }

   for(snap=1; snap<=maxsnap; snap++){
      sprintf(ifile,"HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap,snap); //open current snap
      sprintf(ofile,"0.700E-0.600Z-EZW.O%d.H%d.wGraphGeod",snap,snap); 
      input=fopen(ifile,"r");
      output=fopen(ofile,"w");

//********************************************************************************
//Section I(a): Identify O's that are a part of a water or Eigen structure and 
//assign flags.
//********************************************************************************
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
	
		     weightmatrix[clines][2]=whold1; //the element index 2 and 3 matches currentmatrix
		     weightmatrix[clines][3]=whold2; //O1-H1 weight should match
                     clines++;
                  }
                  else{ //Should never trigger unless OH- exist       
                      countOindex=1;
                      Oindex=Oxy;
                      Hhold1=Hyd;
		      whold1=weight;
                  }
               }//end else-if
               else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
                  if(Oindex==Oxy){
                     currentmatrix[clines][0]=snap;
                     currentmatrix[clines][1]=Oxy;
                     currentmatrix[clines][2]=Hhold1;
                     currentmatrix[clines][3]=Hhold2;
                     currentmatrix[clines][4]=Hyd;
                     currentmatrix[clines][5]=1; //flagged as eigen or H3O+

		     weightmatrix[clines][2]=whold1;
		     weightmatrix[clines][3]=whold2;
		     weightmatrix[clines][4]=weight;
//Remove previous data on O atom, when H3O+ was no identified
                     currentmatrix[clines-1][0]=0;
                     currentmatrix[clines-1][1]=0;
                     currentmatrix[clines-1][2]=0;
                     currentmatrix[clines-1][3]=0;
                     currentmatrix[clines-1][4]=0;

		     weightmatrix[clines-1][2]=0.0;
		     weightmatrix[clines-1][3]=0.0;
		     weightmatrix[clines-1][4]=0.0;
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
	      if(currentmatrix[z][5]==55){ //if water
                 printf("%d %d %d %d %f %f\n",
                 currentmatrix[z][0], //snap
                 currentmatrix[z][5], //flag
                 currentmatrix[z][1], //O index
                 currentmatrix[z][2], //H1
                 currentmatrix[z][3], //H2
	         weightmatrix[z][2], //O-H1 wgt.
	         weightmatrix[z][3]); //O-H2 wgt.
	      }
	      else{
                 printf("%d %d %d %d %d %d %f %f %f\n",
                 currentmatrix[z][0], //snap
                 currentmatrix[z][5], //flag
                 currentmatrix[z][1], //O index
                 currentmatrix[z][2], //H1
                 currentmatrix[z][3], //H2
                 currentmatrix[z][4], //H3
	         weightmatrix[z][2], //O-H1 wgt.
	         weightmatrix[z][3], //O-H2 wgt.
	         weightmatrix[z][4]); //O-H3 wgt.
	      }
           }
        }//end z-loop
*/
//********************************************************************************
//Section I(b): Identify which Eigesn are a part of a Zundel - when two Eigens
//share a proton. In the the first snap, the first O is assumed to be the 
//original one.
//********************************************************************************
              for(crow=0; crow<clines; crow++){
                 for(coln=2; coln<5; coln++){
                    for(nrow=1; nrow<clines; nrow++){
                       for(ncoln=2; ncoln<5; ncoln++){
                          if(currentmatrix[crow][coln]!=0){ //if matrix element is nonzero
                             if(currentmatrix[crow][coln]==currentmatrix[nrow][ncoln]){ //If any H's atoms match
                                if(currentmatrix[crow][1]!=currentmatrix[nrow][1]){ //and if O' atoms differ - Zundel identified
			           if(currentmatrix[crow][9]!=1){ //If new Zundel is identified
                                      currentmatrix[crow][5]=2; 
                                      currentmatrix[nrow][5]=2;
                                      currentmatrix[crow][6]=currentmatrix[crow][1]; //assume this is the original O atom
                                      currentmatrix[crow][7]=currentmatrix[nrow][1]; //store O partner 
                                      currentmatrix[nrow][6]=currentmatrix[crow][1]; //store original O atom in partner 
                                      currentmatrix[nrow][7]=currentmatrix[nrow][1];
                                      currentmatrix[crow][8]=currentmatrix[crow][coln]; //store shared proton
                                      currentmatrix[nrow][8]=currentmatrix[nrow][ncoln]; //store shared proton

				      weightmatrix[crow][5]=weightmatrix[crow][coln]; //orig O and shared H weight
				      weightmatrix[crow][6]=weightmatrix[nrow][ncoln]; //partner o and shared H weight
				      weightmatrix[nrow][5]=weightmatrix[crow][coln];
				      weightmatrix[nrow][6]=weightmatrix[nrow][ncoln];	
				      currentmatrix[nrow][9]=1; //flag partner O info

//Print O-H* weights in Zundels
/*                 		      printf("%d %d %d %d %d %0.3f %0.3f\n",
                 		      currentmatrix[crow][0], //snap
                 		      currentmatrix[crow][5], //flag
                 		      currentmatrix[crow][8], //shared proton
                 		      currentmatrix[crow][6], //orig O
                 		      currentmatrix[crow][7], //partner O
	         		      weightmatrix[crow][5], //O-H1 wgt.
	         		      weightmatrix[crow][6]); //O-H2 wgt.
*/			           }
                                }
                             }
                          }
                       }
                    }//end coln-loop 
                 }//end nrow-loop
              }//end crow-loop

//********************************************************************************
//Section II: Apply weight cutoffs to stuctures and remove structures that
//do not satisfy conditions.
//********************************************************************************

       for(crow=0; crow<clines; crow++){
	  for(coln=2; coln<5; coln++){
	     if(currentmatrix[crow][0]!=0){
	        if(currentmatrix[crow][5]==55){ //Only water molecules
	           if(weightmatrix[crow][coln]<1.000 && weightmatrix[crow][coln]!=0.000){ //If any O-H bond is less than 1.000, remove
		      currentmatrix[crow][0]=0;
		      currentmatrix[crow][1]=0;
		      currentmatrix[crow][2]=0;
		      currentmatrix[crow][3]=0;
		      currentmatrix[crow][4]=0;
		      currentmatrix[crow][5]=0;

		      weightmatrix[crow][2]=0.0;
		      weightmatrix[crow][3]=0.0;
		      weightmatrix[crow][4]=0.0;
		   }
	        }
                else if(currentmatrix[crow][5]==1){
	           if(weightmatrix[crow][coln]<=0.700 && weightmatrix[crow][coln]!=0.000){ // # If any O-H bond is less than 0.900, remove
		      currentmatrix[crow][0]=0;
		      currentmatrix[crow][1]=0;
		      currentmatrix[crow][2]=0;
		      currentmatrix[crow][3]=0;
		      currentmatrix[crow][4]=0;
		      currentmatrix[crow][5]=0;

		      weightmatrix[crow][2]=0.0;
		      weightmatrix[crow][3]=0.0;
		      weightmatrix[crow][4]=0.0;
		   }
	        }
	        else if(currentmatrix[crow][5]==2){
		   if(weightmatrix[crow][coln]<=0.600 && weightmatrix[crow][coln]!=0.000){ // @ If any O-H bond is less than 0.700. remove
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
		   else if(weightmatrix[crow][coln]==0.000){ //If any O-H bond is 0.000 remove - not a real zundel
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

//********************************************************************************
//Section III: Loop through all Eigens and ensure they meet the cutoff 
//requirement.
//********************************************************************************

     for(crow=0; crow<clines; crow++){
        for(coln=2; coln<5; coln++){
           if(currentmatrix[crow][1]!=0){
              if(currentmatrix[crow][5]==1){ //for all eigens
                 if(weightmatrix[crow][coln]<0.700){ // # if O-H bond has a weight less than 1.000, remove
                    currentmatrix[crow][0]=0;
                    currentmatrix[crow][1]=0;
                    currentmatrix[crow][2]=0;
                    currentmatrix[crow][3]=0;
                    currentmatrix[crow][4]=0;
                    currentmatrix[crow][5]=0;

                    weightmatrix[crow][2]=0.0;
                    weightmatrix[crow][3]=0.0;
                    weightmatrix[crow][4]=0.0;
                 }
              }
           }
        }//end coln-loop
     }//end crow-loop

//********************************************************************************
//Section IV: Print out flagged GraphGeod files
//********************************************************************************

        for(z=0; z<clines; z++){
           if(currentmatrix[z][0]!=0){
              fprintf(output,"%d %d %d %d %d %0.3f %0.3f %0.3f\n",
              currentmatrix[z][5], //flag
              currentmatrix[z][1], //O index
              currentmatrix[z][2], //H1
              currentmatrix[z][3], //H2
              currentmatrix[z][4], //H3
	      weightmatrix[z][2], //O-H1 wgt.
	      weightmatrix[z][3], //O-H2 wgt.
	      weightmatrix[z][4]); //O-H3 wgt.
           }
	}//end z-loop
	fclose(output);

//********************************************************************************
//Copy info from nextmatrix into currentmatrix for proceeding comparisons.
//********************************************************************************
      for(y=0; y<300; y++){
         for(z=0; z<10; z++){
	    currentmatrix[y][z]=0;
         }
      }
      clines=0; 
      countOindex=1;
 
      for(y=0; y<300; y++){
         for(z=0; z<7; z++){	
	    weightmatrix[y][z]=0.0;
	 }
      }
   }//end snap
}//end main-loop
