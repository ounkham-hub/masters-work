#include<stdio.h>
#include<stdlib.h>

//*************************************************************************
//Program created to calculate the coordination number for an O atom and 
//OH weight by employing a user-input cutoff.
//
//Input: R-avg-O#.H#.wGraphGeod
//Output: 
//      coln 1: O index
//      coln 2: Coordination #
//      coln 3: avg. weight
//      (optional) colns 4+: H indices and OH weights
//     
//
//Created by Lelee Ounkham
//Last modified September 13, 2018
//
//*************************************************************************

int main(){
   int snap,Oxy,Hyd,rep,sindex,Oindex,Hhold;
   int x,y,z,maxsnap,countOindex,clines,colns;
   int bcount,blines,crow,bcoln,CNcount,nlines;
   int nrow,ncoln,ecount,zcount,jcount,wcount;
   int twocount,threecount;
   float weight,whold,water,eigen,zundel;
   FILE *input,*output,*output2,*output3,*output4;
   char ifile[100],ofile4[100];

   int CBmatrix[1000][20];
   int newmatrix[1000][7];
   float weightmatrix[1000][20];
   float newweight[1000][7];

   water=0.900;
   eigen=0.700;
   zundel=0.500; //also the defined as the smallest OH weight

//   output=fopen("0.900W-0.700E-0.500Z-count-EZW.txt","a");
//   output2=fopen("0.900W-0.700E-0.500Z-CN-2.txt","a");
//   output3=fopen("0.900W-0.700E-0.500Z-CN-3.txt","a");

   for(snap=1; snap<=61555; snap++){
//Intialize matrices and counters
      for(y=0; y<1000; y++){
	 for(z=0; z<20; z++){
	    CBmatrix[y][z]=0;
	    weightmatrix[y][z]=0.0;
	 }
	 for(x=0; x<7; x++){
	    newmatrix[y][x]=0;
	    newweight[y][x]=0;
	 }
      }
      clines=0;
      colns=2; //starts at to because of position in matrix
      countOindex=1;

      sprintf(ifile,"NZ-avg-O%d.H%d.wGraphGeod",snap,snap); //open current snap
      sprintf(ofile4,"0.900W-0.700E-0.500Z-avg-O%d.H%d.wGraphGeod",snap,snap); //open current snap
      input=fopen(ifile,"r");
      output4=fopen(ofile4,"w");

//Read in input
      while(fscanf(input,"%d %d %d %f\n",&sindex,&Oxy,&Hyd,&weight)==4){
	 if(countOindex==1){ //New line, store info on O atom
	     Oindex=Oxy;
	     Hhold=Hyd;
	     whold=weight;
         }
	 else if(countOindex>1){
	    if(Oindex==Oxy){ //if O index matches the previous one, store info
	       if(Hhold!=Hyd){
	          if(countOindex==2){
	             CBmatrix[clines][0]=snap;
	             CBmatrix[clines][1]=Oindex;
	             CBmatrix[clines][colns]=Hhold;
	             weightmatrix[clines][colns]=whold;
		     colns++;
		   }
		   CBmatrix[clines][0]=snap;
		   CBmatrix[clines][1]=Oxy;
		   CBmatrix[clines][colns]=Hyd;
		   weightmatrix[clines][colns]=weight;
		   Oindex=Oxy;
		   Hhold=Hyd;
		   whold=weight;
		   colns++;
	       }	
	       else if(Hhold==Hyd){ 
//For instances when the duplicates are the first two bonds of O, only count 1
	          CBmatrix[clines][0]=snap;
	          CBmatrix[clines][1]=Oindex;
	          CBmatrix[clines][colns]=Hhold;
	          weightmatrix[clines][colns]=whold;
		  colns++;
	       }
	    }
	    else if(Oindex!=Oxy){ //If O's don't match, then new O appeared - restart loop
	       countOindex=1;
	       colns=2;
	       clines++;

  	       Oindex=Oxy;
	       Hhold=Hyd;
	       whold=weight;
	    }
	 }
	 countOindex++;
      }//end while-loop
      fclose(input);
 
/*    bcount=0;
//Cross check to see if H indices and associated weights matches
      for(y=0; y<clines+1; y++){
	 for(z=0; z<20; z++){
	    if(CBmatrix[y][z]!=0){
//	       printf("%d ",CBmatrix[y][z]);
	       bcount++;
	    }
	 }
//Print out ONLY O index and coordination number
fprintf(output,"%d %d %d\n",snap,CBmatrix[y][1],bcount-2);
	 blines=bcount;
//	 printf("\n");
	 for(z=2; z<blines; z++){
//	    printf("%0.3f ",weightmatrix[y][z]);
	}
	bcount=0;
//	printf("\n");
      }
*/

//*************************************************************************
//Section II: Employ user-input weight cutoff to reduce matrix size.
//*************************************************************************

   for(crow=0; crow<clines+1; crow++){
      for(bcoln=2; bcoln<20; bcoln++){
	if(weightmatrix[crow][bcoln]!=0.000 && CBmatrix[crow][bcoln]!=0){
           if(weightmatrix[crow][bcoln] >= zundel){
//If ALL OH weight is greater than 0.700, could be E,Z,W
/*	       for(z=0; z<20; z++){ //Recycle through each O index
	       if(CBmatrix[crow][z]!=0){
	       if(weightmatrix[crow][bcoln] >=0.990 && weightmatrix[crow][z] >= 0.990){
//Water OH weight cutoff, all OH weights need to satisfy this!
		   printf("water: %d %d %d %0.3f ",
		   CBmatrix[crow][0],
		   CBmatrix[crow][1],
		   CBmatrix[crow][z],
		   weightmatrix[crow][z]);
		   printf("\n");

  	           CBmatrix[crow][bcoln]=0;
		   weightmatrix[crow][bcoln]=0.000; 
  	           CBmatrix[crow][z]=0;
		   weightmatrix[crow][z]=0.000; 

		}
		else if(weightmatrix[crow][bcoln] < 0.990 && weightmatrix[crow][z] < 0.990){
		   if(weightmatrix[crow][bcoln] >= 0.900 && weightmatrix[crow][z] >= 0.900){ //Could only be what we define as an Eigen
		      printf("eigen: %d %d %d %0.3f ",
		      CBmatrix[crow][0],
		      CBmatrix[crow][1],
		      CBmatrix[crow][z],
		      weightmatrix[crow][z]);
		      printf("\n");

  	              CBmatrix[crow][bcoln]=0;
		      weightmatrix[crow][bcoln]=0.000; 
  	              CBmatrix[crow][z]=0;
		      weightmatrix[crow][z]=0.000; 
		   }
		   else if(weightmatrix[crow][bcoln] >= 0.700 && weightmatrix[crow][bcoln] < 0.900 && weightmatrix[crow][z] >= 0.700 && weightmatrix[crow][z] < 0.900){
//Could only be what we define as a Zundel
	              printf("zundel: %d %d %d %0.3f ",
		      CBmatrix[crow][0],
		      CBmatrix[crow][1],
		      CBmatrix[crow][z],
		      weightmatrix[crow][z]);
		      printf("\n");

  	              CBmatrix[crow][bcoln]=0;
		      weightmatrix[crow][bcoln]=0.000; 
  	              CBmatrix[crow][z]=0;
		      weightmatrix[crow][z]=0.000; 
		   }
		   }
		   }
		}//end z-loop
*/	    }
            else if(weightmatrix[crow][bcoln] <= zundel){
//Otherwise, remove bond for the matrix
	        CBmatrix[crow][bcoln]=0;
		weightmatrix[crow][bcoln]=0.000; 
	    }
	 }
      }//end coln-loop
   }//end crow-loop

//*************************************************************************
//Section III: After reducing matrix size, create a new matrix to identify 
//and associated coordination #.
//*************************************************************************
   CNcount=3;
   nlines=0;
   for(crow=0; crow<clines+1; crow++){
      for(bcoln=2; bcoln<20; bcoln++){
	 if(CBmatrix[crow][bcoln]!=0 && weightmatrix[crow][bcoln]!=0.000){ //for nonzero elements
	    newmatrix[nlines][0]=snap;
	    newmatrix[nlines][1]=CBmatrix[crow][1];
	    newmatrix[nlines][CNcount]=CBmatrix[crow][bcoln];
	    newweight[nlines][CNcount]=weightmatrix[crow][bcoln];
	    CNcount++;
	 }
      }//end bcoln-loop
      newmatrix[nlines][2]=CNcount-3;
      CNcount=3;
      nlines++;
   }//end crow-loop 

//newmatrix[ ][0] = snapshot
//newmatrix[ ][1] = O index
//newmatrix[ ][2] = Coord. #
//newmatrix[ ][3] to [5] = H index
//newweight[ ][3] to [5] = OH weight

/*
   for(x=0; x<nlines; x++){
      for(y=0; y<6; y++){
         if(newmatrix[x][y]!=0){
	    printf("%d ",newmatrix[x][y]); 
         }
      }
      for(y=0; y<6; y++){
         if(newweight[x][y]!=0.000){
	    printf("%0.3f ",newweight[x][y]);
         }
      }
      printf("\n");
   }
*/
//*************************************************************************
//Section IV: Identify Eigens and Water molecule based on 
//user-input definitions (weight and CN) across all OH bonds.
//*************************************************************************
   for(nrow=0; nrow<nlines; nrow++){
      for(ncoln=3; ncoln<6; ncoln++){
         if(newmatrix[nrow][0]!=0){
            if(newmatrix[nrow][2]==2){ //Should only be the case of a water molecule
               if(newweight[nrow][3] >= water && newweight[nrow][4] >= water){
/*	          printf("water: %d %d %d %d %d %0.3f %0.3f\n",
	          newmatrix[nrow][0],
	          newmatrix[nrow][2],
	          newmatrix[nrow][1],
	          newmatrix[nrow][3],
	          newmatrix[nrow][4],
	          newweight[nrow][3],	
	          newweight[nrow][4]); 
*/	
	          newmatrix[nrow][6]=55;

//	          newmatrix[nrow][0]=0;
	       } 
            }
            else if(newmatrix[nrow][2]==3){ //Distinguishing Eigen from Zundel
	       for(y=0; y<nlines; y++){
	          for(z=3; z<6; z++){
		     if(newmatrix[y][0]!=0){
		        if(newmatrix[nrow][ncoln]==newmatrix[y][z] && newmatrix[nrow][1]!=newmatrix[y][1]){ 
//If any H indices matches and O's do NOT match - Zundel based on CN
/*	          	   printf("Zundel 1: %d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	          	   newmatrix[nrow][0],
	          	   newmatrix[nrow][2],
	          	   newmatrix[nrow][1],
	          	   newmatrix[nrow][3],
	          	   newmatrix[nrow][4],
	          	   newmatrix[nrow][5],
	          	   newweight[nrow][3],	
	          	   newweight[nrow][4],
	          	   newweight[nrow][5]);
 
	          	   printf("Zundel 2: %d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	          	   newmatrix[y][0],
	          	   newmatrix[y][2],
	          	   newmatrix[y][1],
	          	   newmatrix[y][3],
	          	   newmatrix[y][4],
	          	   newmatrix[y][5],
	          	   newweight[y][3],	
	          	   newweight[y][4],
	          	   newweight[y][5]); 
*/
	          	   newmatrix[nrow][6]=2;
	          	   newmatrix[y][6]=2;

//	          	   newmatrix[nrow][0]=0;
		        }  
		     }
		  }//end z-loop
	       }//end y-loop
               if(newweight[nrow][3] >= eigen && newweight[nrow][4] >= eigen && newweight[nrow][5] >= eigen){
		  if(newmatrix[nrow][6]!=2){ //if Zundel not flagged!
/*	          printf("Eigen: %d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	          newmatrix[nrow][0],
	          newmatrix[nrow][2],
	          newmatrix[nrow][1],
	          newmatrix[nrow][3],
	          newmatrix[nrow][4],
	          newmatrix[nrow][5],
	          newweight[nrow][3],	
	          newweight[nrow][4],
	          newweight[nrow][5]); 
*/
	             newmatrix[nrow][6]=1;

//	             newmatrix[nrow][0]=0;
		  }
	       }
	    }
	 }
      }//end ncoln-loop
   }//end nrow-loop

//*************************************************************************
//Section V: Identify whether a structure is Zundel or not.
//*************************************************************************
      for(nrow=0; nrow<nlines; nrow++){
	 if(newmatrix[nrow][0]!=0 && newmatrix[nrow][6]==0){ //for non-classified species
            if(newmatrix[nrow][2]==3){
/*	       fprintf(output3,"%d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
	       newmatrix[nrow][0],
	       newmatrix[nrow][2],
	       newmatrix[nrow][1],
	       newmatrix[nrow][3],
	       newmatrix[nrow][4],
	       newmatrix[nrow][5],
	       newweight[nrow][3],	
	       newweight[nrow][4],
	       newweight[nrow][5]); 
*/	   }
	   else if(newmatrix[nrow][2]==2){
/*	       fprintf(output2,"%d %d %d %d %d %0.3f %0.3f\n",
	       newmatrix[nrow][0],
	       newmatrix[nrow][2],
	       newmatrix[nrow][1],
	       newmatrix[nrow][3],
	       newmatrix[nrow][4],
	       newweight[nrow][3],	
	       newweight[nrow][4]);
*/	    }
	       newmatrix[nrow][6]=7;

	       newmatrix[nrow][0]=0;
	 }
      }
//*************************************************************************
//Section VI: Quantify the number of Eigens, Zundels, waters and everything
//that isn't categorized as one of those structures. Will help
//elucidate weight definitions for future studies.
//*************************************************************************
      ecount=0;
      zcount=0;
      jcount=0;
      wcount=0;
      twocount=0;
      threecount=0;

      for(nrow=0; nrow<nlines; nrow++){
	 if(newmatrix[nrow][6]==1){ //Eigen flag
	    ecount++;
	 }
	 else if(newmatrix[nrow][6]==2){ //Zundel flag	
	    zcount++;
	 } 
	 else if(newmatrix[nrow][6]==7){ //Others flag
	    jcount++;
	    if(newmatrix[nrow][2]==2){ //CN 2
               twocount++;
            }
	    else if(newmatrix[nrow][2]==3){ //CN 3
	       threecount++;
	    }
	 }
	 else if(newmatrix[nrow][6]==55){ //Water flag
	    wcount++;
	 }
      }

//Section to print out total number of Zundels,Eigens,waters, and everthing that didn't satisfy criteria
/*     fprintf(output,"%d %d %d %d %d %d %d\n",
     snap,
     ecount,
     zcount,
     wcount,
     jcount,
     twocount,
     threecount); 
*/

//Section to print out GraphGeod files and flag Zundels and Eigens and necessary
      for(nrow=0; nrow<nlines; nrow++){
	 if(newmatrix[nrow][0]!=0 && newmatrix[nrow][6]!=0){ //for non-classified species
	       fprintf(output4,"%d %d %d %d %d %0.3f %0.3f %0.3f\n",
	       newmatrix[nrow][6], //flag
	       newmatrix[nrow][1], //O index
	       newmatrix[nrow][3], //H1
	       newmatrix[nrow][4], //H2
	       newmatrix[nrow][5], //H3
	       newweight[nrow][3], //wgt1	
	       newweight[nrow][4], //wgt2
	       newweight[nrow][5]); //wg3
        }
     }
     fclose(output4);

   }//end snap-loop
//   fclose(output);
//   fclose(output2);
//   fclose(output3);
}//end main-loop
