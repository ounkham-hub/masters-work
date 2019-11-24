#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//*******************************************************************************
//Program created to take average weights from connections among 32 beads 
//to construct a new set of weighted networks.
//
//Input: #.O#.H#.wGraph
//Output: avg.O#.H#.WGraphGeod
//
//Created by Lelee Ounkham
//Last modified September 5, 2019
//
//*******************************************************************************

int main(){
   int snap,bead,Oxy,Hyd,y,z,num;
   int blines,mlines,brow,mrow,mcoln;
   float weight,sum;
   FILE *input,*output;
   char ifile[100],ofile[100];

   int mastermatrix[2000][3];
   int beadmatrix[1000][2];
   float masterweight[2000][1][33];
   float beadweight[1000][2];

//*******************************************************************************
//Section I: Initialize mastermatrix and beadmatrix prior to analysis.
//*******************************************************************************
   for(y=0; y<2000; y++){
      for(z=0; z<3; z++){
	 mastermatrix[y][z]=0;
      }
   }
   for(y=0; y<2000; y++){
      for(z=0; z<33; z++){
	  masterweight[y][0][z]=0.0;
      }
   }
   for(y=0; y<1000; y++){
      for(z=0; z<3; z++){
	 beadmatrix[y][z]=0;
	 beadweight[y][0]=0.0;
      }
   }
   blines=0;
   mlines=0;

//*******************************************************************************
//Section II.a: Store connections into beadmatrix and compare to mastermatrix.
//*******************************************************************************

   for(snap=1; snap<=2; snap++){
      for(bead=1; bead<=3; bead++){
         sprintf(ifile,"%d.O%d.H%d.wGraph",bead,snap,snap); //open current snap
         input=fopen(ifile,"r");

         while(fscanf(input,"%d %d %f\n",&Oxy,&Hyd,&weight)==3){
	    if(bead==1){
	        mastermatrix[mlines][0]=Oxy;
	        mastermatrix[mlines][1]=Hyd;
	        mastermatrix[mlines][2]=1; //Count starts at 1
	        masterweight[mlines][0][1]=weight; //store weight for bead
	        mlines++; 

//              printf("%d %d %0.3f\n",Oxy,Hyd,weight);     
            } 
	    else if(bead!=1){
	       beadmatrix[blines][0]=Oxy;
	       beadmatrix[blines][1]=Hyd;
	       beadweight[blines][0]=weight;
	       blines++;
	    }
	 }//end while-loop
         fclose(input);

         for(brow=0; brow<blines; brow++){
	    for(mrow=0; mrow<mlines; mrow++){
	       if(beadmatrix[brow][0]!=0){
	          if(beadmatrix[brow][0]==mastermatrix[mrow][0] && beadmatrix[brow][1]==mastermatrix[mrow][1]){
//If O and H atom matches, connection is agrees
		      mastermatrix[mrow][2]=mastermatrix[mrow][2]+1; //increment count
		      num=bead; //store bead number 
		      masterweight[mrow][0][num]=beadweight[brow][0];
		      beadmatrix[brow][0]=0;  
		  }  
	       }
	    } 
	 }//end brow-loop
     
//*******************************************************************************
//Section II.b: Store new connections into mastermatrix
//*******************************************************************************
	 for(brow=0; brow<blines; brow++){
	    if(beadmatrix[brow][0]!=0){
	       mastermatrix[mlines][0]=beadmatrix[brow][0]; //Store O index
	       mastermatrix[mlines][1]=beadmatrix[brow][1]; //Store H index
	       mastermatrix[mlines][2]=1;
	       num=bead;
	       masterweight[mlines][0][num]=beadweight[brow][0]; //Store weight
	       mlines++;
	    }
	 }
//*******************************************************************************
//Section III: Print out averaged weighted networks
//*******************************************************************************
	   for(mrow=0; mrow<mlines; mrow++){
	      for(mcoln=1; mcoln<33; mcoln++){
	   	 sum=sum+masterweight[mrow][0][mcoln]; //sum all weights 
	      }
	      masterweight[mrow][0][0]=sum;
	      sum=0;
	   }
	   
           for(y=0; y<blines; y++){
              for(z=0; z<2; z++){
	         beadmatrix[y][z]=0;
	         beadweight[y][0]=0.0;   
           }
        } 
	blines=0;
     }//end bead-loop 
	 
     for(mrow=0; mrow<mlines; mrow++){
	 printf("%d %d %d %d %0.3f\n",
	 snap,
	 mastermatrix[mrow][0],
         mastermatrix[mrow][1],
   	 mastermatrix[mrow][2],
	 masterweight[mrow][0][0]/32);
     }
   
     for(y=0; y<mlines; y++){
        for(z=0; z<3; z++){
  	       mastermatrix[y][z]=0;
        }
     }
	 
     for(y=0; y<mlines; y++){
        for(z=0; z<33; z++){
  	       masterweight[y][0][z]=0.0;
        }
     } 
     mlines=0;
   }//end snap-loop

}//end main-loop
