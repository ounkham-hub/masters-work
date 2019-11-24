#include<stdio.h>
#include<stdlib.h>

//*************************************************************
//Modified: JUST SEPERATE EIGENS FROM ZUNDELS
//Program written to identify Eigens and Zundels in a single
//snapshot in i(organized) covalent bond GraphGeods.
//	- needed to calculate associated lifetimes 
//
//Note the number prefixes:
//	- 70(real Oindex)
//	- 80(real Oindex)
//	- 9(shared H)
//	   example: O(40) and O(97) in a Zundel -> 7040 and 7097
//	   example: O(39) in an Eigen -> 8039
//	   example: H(132) in E/Z -> 9132
//
//Zundels identified when two O's share an H, otherwise it's
//an Eigen.
//
//Made by Lelee Ounkham
//Last Updated: November 13, 2017
//**************************************************************

int main (){
//Setting up parameters for indexarray
   int Oindex,Hindex,numwaters;
   float ang,dist;
   int PBC1,PBC2,PBC3,PBC4,PBC5;
   int snap,index,lines,countOindex,countwater;
   int row,row2,coln1,coln2,i,y,z;
   int origGraphGeod[400][7];
   int OCOmatrix[400][8];
   float origdistangle[400][2];
   float OCOdistangle[400][2];
   FILE *input,*eoutput,*zoutput;
   char ifile[100],efile[100],zfile[100];

   countwater=0;
   for(snap=1; snap<=50001; snap++){

      sprintf(ifile,"R-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",snap,snap);
      sprintf(efile,"E-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",snap,snap);
      sprintf(zfile,"Z-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",snap,snap);
      input=fopen(ifile,"r");
      eoutput=fopen(efile,"w");
      zoutput=fopen(zfile,"w");
      countOindex=1; //Count the Oindex
      lines=0;
     

      for(y=0; y<400; y++){
	for(z=0; z<8; z++){
	   origGraphGeod[y][z]=0;
	   origdistangle[y][0]=0.0;
	   origdistangle[y][1]=0.0;

           OCOmatrix[y][z]=0;
	   OCOdistangle[y][0]=0.0;
	   OCOdistangle[y][1]=0.0;
        }//end z-loop
      }//end y-loop

      while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oindex,&Hindex,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&ang)==9){
	 origGraphGeod[lines][0]=Oindex;
	 origGraphGeod[lines][1]=Hindex;
	 origGraphGeod[lines][2]=PBC1;
	 origGraphGeod[lines][3]=PBC2;
	 origGraphGeod[lines][4]=PBC3;
	 origGraphGeod[lines][5]=PBC4;
	 origGraphGeod[lines][6]=PBC5;
	 origdistangle[lines][0]=dist;
	 origdistangle[lines][1]=ang;
	 lines++;
      }//end while-loop
      fclose(input);

//Section 1: Identifying all OCO's and storing into OCOmatrix
      for(row=0; row<lines+1; row++){
	 if(countOindex==1){
	    OCOmatrix[row][0]=origGraphGeod[row][0];
	    OCOmatrix[row][1]=origGraphGeod[row][1];
	    OCOmatrix[row][2]=origGraphGeod[row][2];
	    OCOmatrix[row][3]=origGraphGeod[row][3];
	    OCOmatrix[row][4]=origGraphGeod[row][4];
	    OCOmatrix[row][5]=origGraphGeod[row][5];
	    OCOmatrix[row][6]=origGraphGeod[row][6];
	    OCOdistangle[row][0]=origdistangle[row][0];
	    OCOdistangle[row][1]=origdistangle[row][1];
	 }
	 else if(countOindex==2){
	    if(origGraphGeod[row][0]==OCOmatrix[row-1][0]){ //if the next O matches the proceeding one
 
	       OCOmatrix[row][0]=origGraphGeod[row][0];
	       OCOmatrix[row][1]=origGraphGeod[row][1];
	       OCOmatrix[row][2]=origGraphGeod[row][2];
	       OCOmatrix[row][3]=origGraphGeod[row][3];
	       OCOmatrix[row][4]=origGraphGeod[row][4];
	       OCOmatrix[row][5]=origGraphGeod[row][5];
	       OCOmatrix[row][6]=origGraphGeod[row][6];
	       OCOdistangle[row][0]=origdistangle[row][0];
	       OCOdistangle[row][1]=origdistangle[row][1];
	    }
	    else{ //Only occurs if water is identified
	       countOindex=1;
               countwater++;

	       OCOmatrix[row-1][0]=0;
	       OCOmatrix[row-1][1]=0;
	       OCOmatrix[row-1][2]=0;
	       OCOmatrix[row-1][3]=0;
	       OCOmatrix[row-1][4]=0;
	       OCOmatrix[row-1][5]=0;
	       OCOmatrix[row-1][6]=0;
	       OCOdistangle[row-1][0]=0.0;
	       OCOdistangle[row-1][1]=0.0;

	       OCOmatrix[row-2][0]=0;
	       OCOmatrix[row-2][1]=0;
	       OCOmatrix[row-2][2]=0;
	       OCOmatrix[row-2][3]=0;
	       OCOmatrix[row-2][4]=0;
	       OCOmatrix[row-2][5]=0;
	       OCOmatrix[row-2][6]=0;
	       OCOdistangle[row-2][0]=0.0;
	       OCOdistangle[row-2][1]=0.0;

    	       OCOmatrix[row][0]=origGraphGeod[row][0];
	       OCOmatrix[row][1]=origGraphGeod[row][1];
	       OCOmatrix[row][2]=origGraphGeod[row][2];
	       OCOmatrix[row][3]=origGraphGeod[row][3];
	       OCOmatrix[row][4]=origGraphGeod[row][4];
	       OCOmatrix[row][5]=origGraphGeod[row][5];
	       OCOmatrix[row][6]=origGraphGeod[row][6];
	       OCOdistangle[row][0]=origdistangle[row][0];
	       OCOdistangle[row][1]=origdistangle[row][1];
	    }
	 }
	 else if (countOindex==3){
	    if(origGraphGeod[row][0]==OCOmatrix[row-1][0]){ //OCO identified, info saved in OCOmatrix
	       OCOmatrix[row][0]=origGraphGeod[row][0];
	       OCOmatrix[row][1]=origGraphGeod[row][1];
	       OCOmatrix[row][2]=origGraphGeod[row][2];
	       OCOmatrix[row][3]=origGraphGeod[row][3];
	       OCOmatrix[row][4]=origGraphGeod[row][4];
	       OCOmatrix[row][5]=origGraphGeod[row][5];
	       OCOmatrix[row][6]=origGraphGeod[row][6];
	       OCOdistangle[row][0]=origdistangle[row][0];
	       OCOdistangle[row][1]=origdistangle[row][1];
	       countOindex=0; //restart counter
	    }
	    else{ //Water identified - remove from OCOmatrix
	       countOindex=1;
               countwater++;

	       OCOmatrix[row-1][0]=0;
	       OCOmatrix[row-1][1]=0;
	       OCOmatrix[row-1][2]=0;
	       OCOmatrix[row-1][3]=0;
	       OCOmatrix[row-1][4]=0;
	       OCOmatrix[row-1][5]=0;
	       OCOmatrix[row-1][6]=0;
	       OCOdistangle[row-1][0]=0.0;
	       OCOdistangle[row-1][1]=0.0;

	       OCOmatrix[row-2][0]=0;
	       OCOmatrix[row-2][1]=0;
	       OCOmatrix[row-2][2]=0;
	       OCOmatrix[row-2][3]=0;
	       OCOmatrix[row-2][4]=0;
	       OCOmatrix[row-2][5]=0;
	       OCOmatrix[row-2][6]=0;
	       OCOdistangle[row-2][0]=0.0;
	       OCOdistangle[row-2][1]=0.0;

	       OCOmatrix[row][0]=origGraphGeod[row][0];
	       OCOmatrix[row][1]=origGraphGeod[row][1];
	       OCOmatrix[row][2]=origGraphGeod[row][2];
	       OCOmatrix[row][3]=origGraphGeod[row][3];
	       OCOmatrix[row][4]=origGraphGeod[row][4];
	       OCOmatrix[row][5]=origGraphGeod[row][5];
	       OCOmatrix[row][6]=origGraphGeod[row][6];
	       OCOdistangle[row][0]=origdistangle[row][0];
	       OCOdistangle[row][1]=origdistangle[row][1];
	     }
	 }
	 countOindex++;
      }//end row-loop

/*
//Print out O-H's participating in OCO
      for(z=0; z<lines; z++){
	 if(OCOmatrix[z][0]!=0){ //nonzero element
	    printf("%d %d %d %d %d %d %d %.3f %.2f\n",
	    OCOmatrix[z][0],
	    OCOmatrix[z][1],
	    OCOmatrix[z][2],
	    OCOmatrix[z][3],
	    OCOmatrix[z][4],
	    OCOmatrix[z][5],
	    OCOmatrix[z][6],
	    OCOdistangle[z][0],
	    OCOdistangle[z][1]);
	 }
      }//end z-loop
*/

//Section 2: Distinguishing Eigens from Zundels
      for(row=0; row<lines-1; row++){
         for(row2=1; row2<lines; row2++){  
            if(OCOmatrix[row][0]!=0 && OCOmatrix[row2][0]!=0){ //if matrix element is not 0
	       if(OCOmatrix[row][1]==OCOmatrix[row2][1]){//If H's match
		  if(OCOmatrix[row][0]!=OCOmatrix[row2][0]){//If O's DO NOT match - Zundel identified!
//Flag O's that participate in Zundel with a 1
  		     for(i=0; i<lines; i++){
			if(OCOmatrix[row][0]==OCOmatrix[i][0]){
		           OCOmatrix[i][7]=1;
			}
		 	else if(OCOmatrix[row2][0]==OCOmatrix[i][0]){
			   OCOmatrix[i][7]=1;
			}
		     }//end i-loop
		  }
	       }
	    }
	 }//end row2 
      }//end row
/*
//Relabel Zundels and Eigens with appropriate prefixes (i.e. +7000 for Z, +8000 for E, +9000 for H)

      for(row=0; row<lines; row++){
	 if(OCOmatrix[row][0]!=0){//if matrix element is nonzero
	    if(OCOmatrix[row][7]==1){ //Zundel flag
	        OCOmatrix[row][0]=OCOmatrix[row][0]+7000; //relabel Zundel O index
		OCOmatrix[row][1]=OCOmatrix[row][1]+9000; //relabel H index
	    }
	    else{ //Eigen flag
	        OCOmatrix[row][0]=OCOmatrix[row][0]+8000; //relabel Eigen O index
		OCOmatrix[row][1]=OCOmatrix[row][1]+9000; //relabel H
	    }
	 }
      }//end row-loop
*/
//Print out O-H's participating in OCO
      for(z=0; z<lines; z++){
	 if(OCOmatrix[z][0]!=0){ //nonzero element
 	    if(OCOmatrix[z][7]==1){ //O-H is apart of Zundel
	       fprintf(zoutput,"%d %d %d %d %d %d %d %.3f %.2f\n",
	       OCOmatrix[z][0],
	       OCOmatrix[z][1],
	       OCOmatrix[z][2],
	       OCOmatrix[z][3],
	       OCOmatrix[z][4],
	       OCOmatrix[z][5],
	       OCOmatrix[z][6],
	       OCOdistangle[z][0],
	       OCOdistangle[z][1]);
	    }
	    else if(OCOmatrix[z][7]==0){ //O-H is a part of Eigen
	       fprintf(eoutput,"%d %d %d %d %d %d %d %.3f %.2f\n",
	       OCOmatrix[z][0],
	       OCOmatrix[z][1],
	       OCOmatrix[z][2],
	       OCOmatrix[z][3],
	       OCOmatrix[z][4],
	       OCOmatrix[z][5],
	       OCOmatrix[z][6],
	       OCOdistangle[z][0],
	       OCOdistangle[z][1]);
	    }
	 }
      }//end z-loop
      fclose(eoutput);  
      fclose(zoutput); 

   }//end snap-loop
//   numwaters=countwater/2; //there are two covalent bonds in a water
   printf("%d\n",countwater);
}//end main
