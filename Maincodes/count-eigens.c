#include<stdio.h>
#include<stdlib.h>

//***********************************************************************
//Program created to count the number of Eigens and Zundels per
//snapshot. 
//
//Input: E- and Z-HCl.input.O#.xyz.H#.xyz.GraphGeods
//Output: 
//   Contains two columns:
//      First coln: Snapshot #
//      Second coln: # of Eigens/Zundels
//
//Made by Lelee Ounkham
//Last Updated: Novemeber 20, 2017
//
//***********************************************************************


int main(){
//Setting up parameters

   int Oindex,Hindex,snap;
   float ang,dist;
   int PBC1,PBC2,PBC3,PBC4,PBC5;
   int lines,numeigens;
   FILE *input,*output;
   char ifile[100]; 

      output=fopen("total-eigen-vs-time.txt","a");

   for(snap=1; snap<=50001; snap++){

      sprintf(ifile,"E-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",snap,snap);
      input=fopen(ifile,"r");
      lines=0;
    
      while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oindex,&Hindex,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&ang)==9){
         lines++; //Count number of lines in GG, corresponds to one O-H bond

      }//end while
      fclose(input);

      if(lines==0){//if no Eigens
          fprintf(output,"%d 0\n",snap);
      }
      else if(lines!=0){
          numeigens=lines/3; //3 O-H bonds per Eigen
          fprintf(output,"%d %d\n",snap,numeigens);
      }	
    }//end snap-loop
    fclose(output);

}//end main
