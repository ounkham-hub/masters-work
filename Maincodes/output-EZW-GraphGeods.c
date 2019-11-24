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

  int Oxy,Hyd;
  int snap,maxsnap,Oindex,Hhold1,Hhold2;
  int y,z,clines,countOindex;
  int crow,coln,nrow,ncoln;
  float dist,weight,whold1,whold2;
  FILE *input,*output,*woutput,*eoutput,*zoutput;
  char ifile[100],ofile[100],efile[100],zfile[100],wfile[100];

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
//      sprintf(ofile,"EZW.O%d.H%d.wGraphGeod",snap,snap);
      sprintf(efile,"w-Eigen-%d.wGraphGeod",snap); 
      sprintf(wfile,"w-Water-%d.wGraphGeod",snap);  
      sprintf(zfile,"w-Zundel-%d.wGraphGeod",snap);
    
      input=fopen(ifile,"r");
//      output=fopen(ofile,"w");
      eoutput=fopen(efile,"w");
      woutput=fopen(wfile,"w");
      zoutput=fopen(zfile,"w");

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
/*
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
*/

        for(z=0; z<clines; z++){
           if(currentmatrix[z][0]!=0 && currentmatrix[z][5]==55){
              fprintf(woutput,"%d %d %d %d %0.3f %0.3f\n",
              currentmatrix[z][5], //flag
              currentmatrix[z][1], //O index
              currentmatrix[z][2], //H1
              currentmatrix[z][3], //H2
	      weightmatrix[z][2], //O-H1 wgt.
	      weightmatrix[z][3]); //O-H2 wgt.
	   }
	}
	fclose(woutput);

        for(z=0; z<clines; z++){
           if(currentmatrix[z][0]!=0 && currentmatrix[z][5]==1){
              fprintf(eoutput,"%d %d %d %d %d %0.3f %0.3f %0.3f\n",
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
	fclose(eoutput);

        for(z=0; z<clines; z++){
           if(currentmatrix[z][0]!=0 && currentmatrix[z][5]==2){
              fprintf(zoutput,"%d %d %d %d %d %0.3f %0.3f %0.3f\n",
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
	fclose(zoutput);

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
