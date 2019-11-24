#include<stdio.h>
#include<stdlib.h>

//***********************************************************************
//Program created to examine list of unreactive H atoms in the PIMD 
//dataset. IF eigen is associated with an]y other stucture other than
//water then a reaction occurred and line is output.
//
//Created by Lelee Ounkham
//Last Modified: May 13, 2018
//
//***********************************************************************

int main(){ 
   int snap,snapindex,H,lines,hlines;
   int Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5;
   int x,y,z,storeindex,slist,index,maxsnap;
   int Oindex,countOindex,Hhold1,Hhold2;
   int clines,srow,erow,ecoln,hrow;
   int crow,coln,nrow,ncoln;
   float dist,angle; 
   FILE *input,*hinput,*unrinput,*output;
   char ifile[100];

   int snaplist[600][3];
   int Hmatrix[20];
   int tempmatrix[1][3];
   int EZWmatrix[300][6];

//Initialize matrices and counters
    lines=0;
    hlines=0;
    clines=0;
    countOindex=0;
    maxsnap=1000;

     output=fopen("list-reactive-events.txt","a");
//Store H index and selected snapshots 
     unrinput=fopen("extracted-unreactive-snaps.txt","r");
        while(fscanf(unrinput,"%d\n",&snapindex)==1){
           snaplist[lines][0]=snapindex;
           snaplist[lines][1]=snapindex-20;
           snaplist[lines][2]=snapindex+20;
           lines++;
        }//end while-loop
        fclose(unrinput);
  
//      for(x=0; x<lines; x++){
//        printf("%d %d %d\n",snaplist[x][0],snaplist[x][1],snaplist[x][2]);
//      }

     hinput=fopen("H-index.txt","r");
        while(fscanf(hinput,"%d\n",&H)==1){
           Hmatrix[hlines]=H;
           hlines++;
        }//end while-loop
        fclose(hinput);

//     for(x=0; x<hlines; x++){
//       printf("%d\n",Hmatrix[x]);
//     }

//Store GraphGeod into a matrix and flag H type (i.e. water=55, zundel=2, eigen=1)
    for(srow=0; srow<lines; srow++){
       for(snap=1; snap<maxsnap; snap++){
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

                    EZWmatrix[clines][0]=snap;
                    EZWmatrix[clines][1]=Oxy;
                    EZWmatrix[clines][2]=Hhold1;
                    EZWmatrix[clines][3]=Hhold2;
                    EZWmatrix[clines][5]=55; //flagged as water
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
                   EZWmatrix[clines][0]=snap;
                   EZWmatrix[clines][1]=Oxy;
                   EZWmatrix[clines][2]=Hhold1;
                   EZWmatrix[clines][3]=Hhold2;
                   EZWmatrix[clines][4]=Hyd;
                   EZWmatrix[clines][5]=1; //flagged as eigen or H3O+

                   EZWmatrix[clines-1][0]=0;
                   EZWmatrix[clines-1][1]=0;
                   EZWmatrix[clines-1][2]=0;
                   EZWmatrix[clines-1][3]=0;
                   EZWmatrix[clines-1][4]=0;
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

//***********************************************************************
//Section to identify structures involved in Zundels
//***********************************************************************
      for(crow=0; crow<clines; crow++){
         for(coln=2; coln<5; coln++){
            for(nrow=1; nrow<clines; nrow++){
               for(ncoln=2; ncoln<5; ncoln++){
                  if(EZWmatrix[crow][coln]!=0){ //if matrix element is nonzero
                     if(EZWmatrix[crow][coln]==EZWmatrix[nrow][ncoln]){ //If any H's atoms match
                        if(EZWmatrix[crow][1]!=EZWmatrix[nrow][1]){ //and if O' atoms differ - Zundel identified
                           EZWmatrix[crow][5]=2;
                           EZWmatrix[nrow][5]=2;
//                           EZWmatrix[crow][6]=currentmatrix[crow][coln]; //store shared proton
//                           EZWmatrix[nrow][6]=currentmatrix[nrow][ncoln]; //store shared proton
			}
                     }
                  }
               }
            }//end coln-loop 
         }//end nrow-loop
      }//end crow-loop

/*
         for(z=0; z<clines; z++){
	    if(EZWmatrix[z][0]!=0){
	       printf("%d %d %d %d %d %d\n",
	       EZWmatrix[z][0],
	       EZWmatrix[z][1],
	       EZWmatrix[z][2],
	       EZWmatrix[z][3],
	       EZWmatrix[z][4],
	       EZWmatrix[z][5]);
	    }
	 }
*/

//***********************************************************************
//Section to identify structure of H atom at specified snapshot.
//***********************************************************************

             if(snap >= snaplist[srow][1] && snap <= snaplist[srow][2]){ //If snapshot index is within specified +/- 20 window
	        for(erow=0; erow<clines; erow++){
		   for(ecoln=2; ecoln<4; ecoln++){ 
	              for(hrow=0; hrow<hlines; hrow++){
		         if(EZWmatrix[erow][ecoln]==Hmatrix[hrow]){
//	       printf("%d %d %d %d %d %d\n",
//	       EZWmatrix[erow][0],
//	       EZWmatrix[erow][1],
//	       EZWmatrix[erow][2],
//	       EZWmatrix[erow][3],
//	       EZWmatrix[erow][4],
//	       EZWmatrix[erow][5]);
			    if(EZWmatrix[erow][5]!=55){ //Only print when H atom is NOT a water
                               fprintf(output,"%d %d %d %d\n",
			       snap,
			       EZWmatrix[erow][1],
			       EZWmatrix[erow][ecoln],
			       EZWmatrix[erow][5]);
			    }
			 }
		      }//end hrow-loop
	           }//end ecoln-loop
		}//end erow-loop
	     }
	     else{ //Is NOT within window - do NOT analyze and initialize matrices and counters
                for(z=0; z<clines; z++){
	           for(y=0; y<6; y++){
	              EZWmatrix[z][y]=0;
	           }//end y-loop
                }//end z-loop
	        clines=0;
		countOindex=0;
	     }

//Initialize GraphGeod matrix
          for(z=0; z<clines; z++){
	     for(y=0; y<6; y++){
	        EZWmatrix[z][y]=0;
	     }//end y-loop
          }//end z-loop
	  clines=0;
	  countOindex=0;
      }//end snap-loop
   }//end srow-loop
   fclose(output);
}// end main-loop
