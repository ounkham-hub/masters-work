#include<stdio.h>
#include<stdlib.h>

//*****************************************************************************
//Program created to identify all waters,eigens, and zundels and associated
//weights in each snapshot.
//
//Flags for each structure:
//Eigens 1
//Zundels 2
//Waters 55
//
//Input: HCl.covalent.O#.xyz.H#.xyz.wGraphGeod
//Output: EZW-#.wGraphGeod
//
//Created by Lelee Ounkham
//Last Modified May 23, 2018
//
//******************************************************************************

int main(){

  int Oxy,Hyd,index,crow,coln,nrow,ncoln;
  int snap,Oindex,Hhold1,Hhold2;
  int y,z,clines,countOindex;
  FILE *input,*eoutput,*zoutput;

  int currentmatrix[4800000][7];

//Change variables accordingly 
   clines=0; 
   countOindex=1;

//Intialize matrices
   for(y=0; y<4800000; y++){
      for(z=0; z<7; z++){
	 currentmatrix[y][z]=0;
      }
   }

   input=fopen("compiled-1-20000-unreactive-centroid.txt","r");
   eoutput=fopen("complete-list-of-UNR-eigen.txt","a");
   zoutput=fopen("complete-list-of-UNR-zundel.txt","a");
   spoutput=fopen("complete-list-of-UNR-SP.txt","a");

//********************************************************************************
//Section I(a): Identify O's that are a part of a water or Eigen structure and 
//assign flags.
//********************************************************************************
            while(fscanf(input,"%d %d %d\n",&index,&Oxy,&Hyd)==3){
               if(countOindex==1){
		  snap=index;
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
		      snap=index;      
                      countOindex=1;
                      Oindex=Oxy;
                      Hhold1=Hyd;
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
		     snap=index;
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
	      if(currentmatrix[z][5]==55){ //if water
                 printf("%d %d %d %d %0.3f %0.3f\n",
                 currentmatrix[z][0], //snap
                 currentmatrix[z][5], //flag
                 currentmatrix[z][1], //O index
                 currentmatrix[z][2], //H1
                 currentmatrix[z][3], //H2
	         weightmatrix[z][2], //O-H1 wgt.
	         weightmatrix[z][3]); //O-H2 wgt.
	      }
	      if(currentmatrix[z][5]!=55){
                 printf("%d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
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
if(currentmatrix[crow][0]==currentmatrix[nrow][0]){ //if snapshot index matches
                             if(currentmatrix[crow][coln]==currentmatrix[nrow][ncoln]){ //If any H's atoms match
                                if(currentmatrix[crow][1]!=currentmatrix[nrow][1]){ //and if O' atoms differ - Zundel identified
			           if(currentmatrix[crow][6]!=1){ //If new Zundel is identified
                                      currentmatrix[crow][5]=2; 
                                      currentmatrix[nrow][5]=2;

				      if(coln==2){// for H1

					 fprintf(spoutout,"%d 2 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][coln], //shared proton);
					 currentmatrix[crow][1]);

					 fprintf(zoutput,"%d 4 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][3], //shared proton);
					 currentmatrix[crow][1]);

					 fprintf(zoutput,"%d 4 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][4], //shared proton);
					 currentmatrix[crow][1]);

				         currentmatrix[crow][6]=1; //flag partner O info
				      }
				      else if(coln==3){ //for H2

					 fprintf(spoutput,"%d 2 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][coln], //shared proton);
					 currentmatrix[crow][1]);

					 fprintf(zoutput,"%d 4 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][2], //shared proton);
					 currentmatrix[crow][1]);

					 fprintf(zoutput,"%d 4 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][4], //shared proton);
					 currentmatrix[crow][1]);

				         currentmatrix[crow][6]=1; //flag partner O info

				      }
				      else if(coln==4){ //for H3
					 fprintf(spoutput,"%d 2 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][coln], //shared proton);
					 currentmatrix[crow][1]);

					 fprintf(zoutput,"%d 4 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][2], //shared proton);
					 currentmatrix[crow][1]);

					 fprintf(zoutput,"%d 4 %d %d\n",
					 currentmatrix[crow][0],//snap
					 currentmatrix[crow][3], //shared proton);
					 currentmatrix[crow][1]);

				         currentmatrix[crow][6]=1; //flag partner O info

				      }

			           }
                                }
                             }
                          }
                       }
                    }//end coln-loop 
                 }//end nrow-loop
              }//end crow-loop
	   }

        for(z=0; z<clines; z++){
           if(currentmatrix[z][0]!=0){
	      if(currentmatrix[z][5]==1){
                 fprintf(eoutput,"%d %d %d %d\n",
                 currentmatrix[z][0], //snap
                 currentmatrix[z][5], //flag
                 currentmatrix[z][1], //O index
                 currentmatrix[z][2]); //H1

                 fprintf(eoutput,"%d %d %d %d\n",
                 currentmatrix[z][0], //snap
                 currentmatrix[z][5], //flag
                 currentmatrix[z][1], //O index
                 currentmatrix[z][3]); //H1

                 fprintf(eoutput,"%d %d %d %d\n",
                 currentmatrix[z][0], //snap
                 currentmatrix[z][5], //flag
                 currentmatrix[z][1], //O index
                 currentmatrix[z][4]); //H1

	      }
           }
        }//end z-loop
	fclose(eoutput);
	fclose(zoutput);
	fclose(spoutput);

}//end main-loop
