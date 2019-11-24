#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//***********************************************************************
//Program created to read file containing all the list of H atoms
//that remain unreactive for more than 10k snapshots and output
//the weights associated with that particular bond.
//
//Input: weighted-unreactive-waters.txt 
//       weighted-unreactive-eigens-0.9cutoff.txt
//       EZW-O#.H#.wGraphGeod
//Output: unr-weight-using-1.0cutoff.txt
//
//Created by Lelee Ounkham
//Last modified June 13, 2018
//
//***********************************************************************

int main(){

  int Hindex,Oindex,flag,initialsnap,finalsnap;
  int H1,H2,H3,maxsnap,elines,y,z,snap,diff;
  int ulines,urow,Oxy,erow,ecoln;
  float weight1,weight2,weight3,tenpercent;
  FILE *input,*unrinput,*output,*zoutput;
  char ifile[100]; 

  int unrmatrix[100000][6]; 
  int ezwmatrix[103][5];
  float weightmatrix[103][5];

//Intialize counters 
   ulines=0;
   elines=0;
   maxsnap=61555;

//***********************************************************************
//Read in list of unreactive H atoms, store into matrix, and apply 
//cutoffs.
//***********************************************************************

   unrinput=fopen("eigen-0.700E-0.500Z-cutoff.txt","r");
//   unrinput=fopen("weighted-unreactive-waters.txt","r");
   output=fopen("indiv-eigen-0.700E-0.500Z-weights.txt","a");
   zoutput=fopen("indiv-zundel-0.700E-0.500Z-weights.txt","a");

   while(fscanf(unrinput,"%d %d %d %d %d %d\n",&Hindex,&Oindex,&flag,&initialsnap,&finalsnap,&diff)==6){
//      if(diff>=500 && flag==55){ //for waters (55) greater than 500 and 30 for Eigens (1)
     if(flag!=55){ //Just store eigens and zundels into matrix
         unrmatrix[ulines][0]=Hindex;
         unrmatrix[ulines][1]=Oindex;
         unrmatrix[ulines][2]=flag;
         unrmatrix[ulines][3]=diff;
      }

//         tenpercent=(diff*0.10); //The first and last 10% will be excluded from analysis
//         unrmatrix[ulines][4]=initialsnap+tenpercent;
//         unrmatrix[ulines][5]=(finalsnap+1)-tenpercent;
         ulines++;
      }
   }
   fclose(unrinput);
/*
   for(z=0; z<ulines; z++){
      printf("%d %d %d %d %d %d\n",
      unrmatrix[z][0],
      unrmatrix[z][1],
      unrmatrix[z][2],
      unrmatrix[z][3],
      unrmatrix[z][4],
      unrmatrix[z][5]);
   }//end z-loop
*/

//***********************************************************************
//Section to read in weighted GraphGeod files.
//***********************************************************************

   for(urow=0; urow<ulines; urow++){
      for(snap=1; snap<=maxsnap; snap++){
         sprintf(ifile,"0.700E-0.500Z-EZW.O%d.H%d.wGraphGeod",snap,snap);
         input=fopen(ifile,"r");

         while(fscanf(input,"%d %d %d %d %d %f %f %f\n",&flag,&Oxy,&H1,&H2,&H3,&weight1,&weight2,&weight3)==8){
            ezwmatrix[elines][0]=flag;
	    ezwmatrix[elines][1]=Oxy;
	    ezwmatrix[elines][2]=H1;
	    ezwmatrix[elines][3]=H2;
	    ezwmatrix[elines][4]=H3;

	    weightmatrix[elines][2]=weight1;
	    weightmatrix[elines][3]=weight2;
	    weightmatrix[elines][4]=weight3; 
            elines++;
         }
         fclose(input); 
 
/*
         for(y=0; y<elines; y++){
	    printf("%d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	    snap,
            ezwmatrix[y][0],
	    ezwmatrix[y][1],
	    ezwmatrix[y][2],
	    ezwmatrix[y][3],
	    ezwmatrix[y][4],
	    weightmatrix[y][2],
	    weightmatrix[y][3],
	    weightmatrix[y][4]);
         }//end y-loop 
*/

//***********************************************************************
//Section to extract weights from unreactive structures.
//***********************************************************************
      for(erow=0; erow<elines; erow++){
         for(ecoln=2; ecoln<5; ecoln++){
	    if(ezwmatrix[erow][0]==1){ //If flag is 55 (water) or 1 (eigen)
	       if(unrmatrix[urow][1]==ezwmatrix[erow][1] && unrmatrix[urow][0]==ezwmatrix[erow][ecoln]){
	          if(snap >= unrmatrix[urow][4] && snap <= unrmatrix[urow][5]){ //if snapshot is within time window
//If O and H index matches list and GraphGeod
		      fprintf(output,"%d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
		      snap,
		      ezwmatrix[erow][0], //flag
		      ezwmatrix[erow][1], //Oindex
		      ezwmatrix[erow][2], //H1
		      ezwmatrix[erow][3], //H2
		      ezwmatrix[erow][4], //H3
		      weightmatrix[erow][2], //weight 1
		      weightmatrix[erow][3], //weight 2
		      weightmatrix[erow][4]); //weight 3
		  }
	       }
	    }
	    if(ezwmatrix[erow][0]==2){ //If flag is 55 (water) or 1 (eigen)
	       if(unrmatrix[urow][1]==ezwmatrix[erow][1] && unrmatrix[urow][0]==ezwmatrix[erow][ecoln]){
	          if(snap >= unrmatrix[urow][4] && snap <= unrmatrix[urow][5]){ //if snapshot is within time window
//If O and H index matches list and GraphGeod
		      fprintf(zoutput,"%d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
		      snap,
		      ezwmatrix[erow][0], //flag
		      ezwmatrix[erow][1], //Oindex
		      ezwmatrix[erow][2], //H1
		      ezwmatrix[erow][3], //H2
		      ezwmatrix[erow][4], //H3
		      weightmatrix[erow][2], //weight 1
		      weightmatrix[erow][3], //weight 2
		      weightmatrix[erow][4]); //weight 3
		  }
	       }
	    }

	 }//end ecoln-loop
      }//end erow-loop
		
//Initialize ezwmatrix and weightmatrix
         for(y=0; y<elines; y++){
            for(z=0; z<5; z++){
	       ezwmatrix[y][z]=0;
	       weightmatrix[y][z]=0.0;
            }//end z-loop
         }//end y-loop
         elines=0; 
      }//end snap-loop
   }//end urow-loop
   fclose(output);
   fclose(zoutput);
}//end main-loop
