#include<stdio.h>
#include<stdlib.h>

//*****************************************************************************
//Program created to identify all Zundels, the time formed until transformed.
//
//Flags for each structure:
//Eigens 1
//Zundels 2
//Waters 55
//
//Input: R-HCl.input.solB#.xyz.solC#.xyz.GraphGeod
//Output: duration-of-zundel.txt
//
//Created by Lelee Ounkham
//Last Modified March 18, 2018
//
//******************************************************************************

int main(){

  int Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5;
  int snap,maxsnap,Oindex,Hhold1,Hhold2;
  int y,z,clines,countOindex,zlines,nlines;
  int crow,coln,nrow,ncoln,zrow;
  float dist,angle;
  FILE *input,*nextinput,*output;
  char ifile[100],nfile[100];

  int currentmatrix[300][10];
  int nextmatrix[300][10];
  int zundelmatrix[300000][6];

//Change variables accordingly 
   maxsnap=61555;
   clines=0; 
   countOindex=1;
   zlines=0;

//Intialize matrices
   for(y=0; y<300; y++){
      for(z=0; z<10; z++){
	 currentmatrix[y][z]=0;
      }
   }

   output=fopen("duration-of-zundel.txt","w");
   for(snap=1; snap<=maxsnap; snap++){
      if(snap==maxsnap){
         for(z=0; z<zlines; z++){
	    fprintf(output,"%d %d %d %d %d\n",
            zundelmatrix[z][0],
	    zundelmatrix[z][1],
	    zundelmatrix[z][2],
	    zundelmatrix[z][3],
	    zundelmatrix[z][4]);
         }
	 fclose(output);
	 exit(0);
      }
      else{
         if(snap==1){
	    sprintf(ifile,"R-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
            input=fopen(ifile,"r");

//********************************************************************************
//Section I(a): Identify O's that are a part of a water or Eigen structure and 
//assign flags.
//********************************************************************************
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
                                      currentmatrix[crow][7]=currentmatrix[nrow][1];//store O partner 
                                      currentmatrix[nrow][6]=currentmatrix[crow][1]; //store original O atom in partner 
                                      currentmatrix[nrow][7]=currentmatrix[nrow][1];
                                      currentmatrix[crow][8]=currentmatrix[crow][coln]; //store shared proton
                                      currentmatrix[nrow][8]=currentmatrix[nrow][ncoln]; //store shared proton
				      currentmatrix[nrow][9]=1; //flag partner O info
			           }
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
              printf("%d %d %d %d %d %d %d %d %d %d\n",
              currentmatrix[z][0],
              currentmatrix[z][1],
              currentmatrix[z][2],
              currentmatrix[z][3],
              currentmatrix[z][4],
              currentmatrix[z][5],
	      currentmatrix[z][6],
	      currentmatrix[z][7],
	      currentmatrix[z][8],
	      currentmatrix[z][9]);
           }
        }//end z-loop
*/

//********************************************************************************
//Add Zundels to zundelmatrix - only original O's 
//********************************************************************************
           for(crow=0; crow<clines; crow++){
	      if(currentmatrix[crow][9]==1){ //if O in apart of Zundel
	         if(currentmatrix[crow][5]==2){ //Only store original O's
	            zundelmatrix[zlines][0]=currentmatrix[crow][8]; //shared proton
                    zundelmatrix[zlines][1]=currentmatrix[crow][6]; //original O atom
	            zundelmatrix[zlines][2]=currentmatrix[crow][7]; //partner O atom
	            zundelmatrix[zlines][3]=currentmatrix[crow][0]; //initial snap
	            zundelmatrix[zlines][4]=currentmatrix[crow][0]; //final snap
	            zlines++;
	         }
	      }
           }//end crow-loop
/*
      for(z=0; z<zlines; z++){
	 printf("%d %d %d %d %d\n",
         zundelmatrix[z][0],
	 zundelmatrix[z][1],
	 zundelmatrix[z][2],
	 zundelmatrix[z][3],
	 zundelmatrix[z][4]);
      }
*/
        }//end snap==1
//********************************************************************************
//Section II: Identify which O's are a part of a water or Eigen structure for 
//the NEXT snapshot.
//********************************************************************************
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

//Remove previous data on O atom, when H3O+ was not identifed
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
            printf("%d %d %d %d %d %d %d %d\n",
            nextmatrix[z][0],
            nextmatrix[z][1],
            nextmatrix[z][2],
            nextmatrix[z][3],
            nextmatrix[z][4],
            nextmatrix[z][5],
            nextmatrix[z][6],
            nextmatrix[z][7],
	    nextmatrix[z][8],
	    nextmatrix[z][9]);
         }
       }//end z-loop
*/
//********************************************************************************
//Section III: Add Zundels to zundelmatrix for NEXT snapshot - only original O's. 
//Need to compare nextmatrix to currentmatrix to seem if new Zundels are formed 
//and keep track of which O index is the corresponds to the original Eigen of
//a Zundel. 
//********************************************************************************
       for(crow=0; crow<nlines; crow++){
          for(coln=2; coln<5; coln++){
             for(nrow=1; nrow<nlines; nrow++){
                for(ncoln=2; ncoln<5; ncoln++){
                   if(nextmatrix[crow][coln]!=0){ //if matrix element is nonzero
                      if(nextmatrix[crow][coln]==nextmatrix[nrow][ncoln]){ //If any H's atoms match
                         if(nextmatrix[crow][1]!=nextmatrix[nrow][1]){ //and if O' atoms differ - Zundel identified
                            nextmatrix[crow][5]=2; 
                            nextmatrix[nrow][5]=2;
                            nextmatrix[crow][6]=nextmatrix[crow][1]; //assume this is the original O atom
                            nextmatrix[crow][7]=nextmatrix[nrow][1];//store O partner 
                            nextmatrix[nrow][6]=nextmatrix[crow][1]; //store original O atom in partner 
                            nextmatrix[nrow][7]=nextmatrix[nrow][1];
                            nextmatrix[crow][8]=nextmatrix[crow][coln]; //store shared proton
                            nextmatrix[nrow][8]=nextmatrix[nrow][ncoln]; //store shared proton
                          }
                       }
                    }
                 }
              }//end coln-loop 
           }//end nrow-loop
        }//end crow-loop
/*
      for(z=0; z<nlines; z++){
         if(nextmatrix[z][0]!=0){
            printf("%d %d %d %d %d %d %d %d %d %d\n",
            nextmatrix[z][0],
            nextmatrix[z][1],
            nextmatrix[z][2],
            nextmatrix[z][3],
            nextmatrix[z][4],
            nextmatrix[z][5],
            nextmatrix[z][6],
            nextmatrix[z][7],
	    nextmatrix[z][8],
	    nextmatrix[z][9]);
         }
       }//end z-loop

*/
//********************************************************************************
//Determine which O index corresponds to the original Eigen in newly formed Zundel
//********************************************************************************
       for(crow=0; crow<clines; crow++){
          for(nrow=0; nrow<nlines; nrow++){
	     if(currentmatrix[crow][1]==nextmatrix[nrow][1]){ //if O indices matches
	        if(currentmatrix[crow][6]==0 && currentmatrix[crow][7]==0){ //if new Zundel formed
	           if(currentmatrix[crow][5]==1 && nextmatrix[nrow][5]==2){
//If O index of an Eigen transforms to Zundel - assume original O
		      if(nextmatrix[nrow][6]!=nextmatrix[nrow][1]){ //must be different partner
//printf("trigger 1: %d %d\n",nextmatrix[nrow][1],nextmatrix[nrow][6],snap+1);
	                 zundelmatrix[zlines][0]=nextmatrix[nrow][8]; //shared proton
                         zundelmatrix[zlines][1]=nextmatrix[nrow][1]; //original O atom
	                 zundelmatrix[zlines][2]=nextmatrix[nrow][6]; //partner O atom
	                 zundelmatrix[zlines][3]=nextmatrix[nrow][0]; //initial snap
	                 zundelmatrix[zlines][4]=nextmatrix[nrow][0]; //final snap
	                 zlines++;
		      }
		      else if(nextmatrix[nrow][7]!=nextmatrix[nrow][1]){
//printf("trigger 2: %d %d\n",nextmatrix[nrow][1],nextmatrix[nrow][7],snap+1);
	                 zundelmatrix[zlines][0]=nextmatrix[nrow][8]; //shared proton
                         zundelmatrix[zlines][1]=nextmatrix[nrow][1]; //original O atom
	                 zundelmatrix[zlines][2]=nextmatrix[nrow][7]; //partner O atom
	                 zundelmatrix[zlines][3]=nextmatrix[nrow][0]; //initial snap
	                 zundelmatrix[zlines][4]=nextmatrix[nrow][0]; //final snap
	                 zlines++;
	              }
		   }
	        }
//********************************************************************************
//For pre-existing Zundels - update final snap in zundelmatrix
//********************************************************************************
		 else if(currentmatrix[crow][6]!=0 && currentmatrix[crow][7]!=0){
		    if(currentmatrix[crow][5]==2 && nextmatrix[nrow][5]==2){ //if O index remains a Zundel in next snap
		       for(zrow=0; zrow<zlines; zrow++){
			  if(nextmatrix[nrow][1]==zundelmatrix[zrow][1] && nextmatrix[nrow][8]==zundelmatrix[zrow][0]){
//If O index matches and shared proton
			     if(nextmatrix[nrow][0]==zundelmatrix[zrow][4]+1){
			        zundelmatrix[zrow][4]=nextmatrix[nrow][0]; //update snap.
			     }		        
			  }
		       }//end zrow-loop
		    }
		 }
	      }
	      else if(currentmatrix[crow][1]!=nextmatrix[nrow][1]){ //If O indices DON'T match
		 if(currentmatrix[crow][6]!=0 && currentmatrix[crow][7]!=0){
		    if(currentmatrix[crow][5]==2 && nextmatrix[nrow][5]==2){
		       for(zrow=0; zrow<zlines; zrow++){
			  
		       }//end zrow-loop
		    }   
		 }
	      }
           }//end nrow-loop
        }//end crow-loop

/*
       for(z=0; z<zlines; z++){
	  printf("%d %d %d %d %d\n",
          zundelmatrix[z][0],
	  zundelmatrix[z][1],
	  zundelmatrix[z][2],
	  zundelmatrix[z][3],
	  zundelmatrix[z][4]);
      }
*/


//********************************************************************************
//Copy info from nextmatrix into currentmatrix for proceeding comparisons.
//********************************************************************************
      for(y=0; y<300; y++){
         for(z=0; z<10; z++){
	    currentmatrix[y][z]=0;
         }
      }

      for(y=0; y<nlines; y++){
         for(z=0; z<10; z++){
            currentmatrix[y][z]=nextmatrix[y][z];
         }
      }
      clines=nlines;

      for(y=0; y<300; y++){
	 for(z=0; z<10; z++){
	    nextmatrix[y][z]=0;
	 }//end z-loop
      }

      }//end else-loop
   }//end snap
}//end main-loop
