#include<stdio.h>
#include<stdlib.h>

//*********************************************************************************
//Program created to remove all predetermiend CBs and associated weights involved
//in proton transfer event. Essentially, an updated list of weighted networks
//will be outputted and will contain the weights of true, nonreactive Eigens and
//Zundels.
//
//Input: PTsnaps-using-20persist.txt
//       H-HCl.covalent.O#.xyz.H#.xyz.wGraphGeod
//
//Created by Lelee Ounkham
//Last modified on March 2, 2019
//
//*********************************************************************************

int main(){

   int H,O1,O2,initransfer,fintransfer,y,z,plines;
   int Oxy,Hyd,wlines,maxsnap,snap,prow,wrow,pcoln;
   int PBC1,PBC2,PBC3,PBC4,PBC5;
   float dist,angle;
   FILE *pinput,*winput,*output,*output2,*input;
   char ifile[100];

   int persistmatrix[2000][5];
   int indexmatrix[4800000][3];

//Initialize counters and matrices

  for(y=0; y<2000; y++){
     for(z=0; z<5; z++){
	persistmatrix[y][z]=0;
     }
  }
  plines=0;

  for(y=0; y<4800000; y++){
     for(z=0; z<3; z++){
	indexmatrix[y][z]=0;
     }
  }
  wlines=0;
  

//Section to store information from input files to appropriate matrices

     output=fopen("compiled-1-20000-unreactive-centroid.txt","w");
//       output=fopen("compiled-40001-61555-centroid.txt","a");


     winput=fopen("compiled-1-20000-centroid.txt","r");
     while(fscanf(winput,"%d %d %d\n",&snap,&Oxy,&Hyd)==3){
	  indexmatrix[wlines][0]=snap;
          indexmatrix[wlines][1]=Hyd;
          indexmatrix[wlines][2]=Oxy;
          wlines++;
       }//end while-loop
       fclose(winput);     

     pinput=fopen("exchange-events-1-20000.txt","r");
     while(fscanf(pinput,"%d %d %d %d %d\n",&H,&O1,&O2,&initransfer,&fintransfer)==5){
         persistmatrix[plines][0]=H;
         persistmatrix[plines][1]=O1;
         persistmatrix[plines][2]=O2;
         persistmatrix[plines][3]=initransfer;
         persistmatrix[plines][4]=fintransfer;
         plines++;
     }//end while-loop
     fclose(pinput);

/*
   for(y=0; y<plines; y++){
      printf("%d %d %d %d %d\n",
      persistmatrix[y][0],
      persistmatrix[y][1],
      persistmatrix[y][2],
      persistmatrix[y][3],
      persistmatrix[y][4]);
   } 
*/

// *This section was used to compile the weighted networks into a 3 separate files:
//  a.  1 - 20,000 snapshots
//  b.  20,001 - 40,000 snapshots
//  c.  40,001 - 61,555 snapshots
/*
   for(snap=40001; snap<=61555; snap++){
      sprintf(ifile,"R-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
      input=fopen(ifile,"r");

       while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
          indexmatrix[wlines][0]=Hyd;
          indexmatrix[wlines][1]=Oxy;
          wlines++;
       }//end while-loop
       fclose(input);

       for(y=0; y<wlines; y++){
	  fprintf(output,"%d %d %d\n",
	  snap,
	  indexmatrix[y][1],
	  indexmatrix[y][0]);
       }

       for(y=0; y<wlines; y++){
	  indexmatrix[y][0]=0;
	  indexmatrix[y][1]=0;
       }
       wlines=0;
   }//end snap-loop
   fclose(output);
*/
//*********************************************************************************
//Section I: Remove weights associated with reactive O and H indices stored in
//persistmatrix. First identify which snapshots are involved in proton transfer 
//from initransfer and fintransfer snapshots.
//
//*********************************************************************************

     for(prow=0; prow<plines; prow++){
	 for(wrow=0; wrow<wlines; wrow++){
            if(indexmatrix[wrow][0] >= persistmatrix[prow][3] && indexmatrix[wrow][0] <= persistmatrix[prow][4]){
	        if(persistmatrix[prow][1]==indexmatrix[wrow][2]){ 
//If current O index matches, remove weight from indexmatrix and weight matrix
/*
		   printf("O1: %d %d %d %d %d\n",
		   persistmatrix[prow][3],
		   persistmatrix[prow][4],
		   indexmatrix[wrow][0],
		   indexmatrix[wrow][1],
		   indexmatrix[wrow][2]);
*/
		   indexmatrix[wrow][0]=0;
		   indexmatrix[wrow][1]=0;
		   indexmatrix[wrow][2]=0;

		}
		else if(persistmatrix[prow][2]==indexmatrix[wrow][2]){
//If previous O index matches, remove weight from indexmatrix and weight matrix
/*		   printf("O2: %d %d %d %d %d\n",
		   persistmatrix[prow][3],
		   persistmatrix[prow][4],
		   indexmatrix[wrow][0],
		   indexmatrix[wrow][1],
		   indexmatrix[wrow][2]);
*/
		   indexmatrix[wrow][0]=0;
		   indexmatrix[wrow][1]=0;
		   indexmatrix[wrow][2]=0;
		}
	     }//end wrow-loop
	  }//end if snap is within windowtime
       }//end prow-loop
//Section to print out updated unreactive set of weights
          for(y=0; y<wlines; y++){
	     if(indexmatrix[y][0]!=0){
               fprintf(output,"%d %d %d\n",
	       indexmatrix[y][0],
               indexmatrix[y][2],
	       indexmatrix[y][1]);
	    }
         }
        fclose(output);

}//end main-loop
