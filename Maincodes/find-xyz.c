#include<stdio.h>
#include<stdlib.h>

//*****************************************************************
//Program created to read O-O GraphGeods and find corresponding
//xyz files per index.
//
//Created by Lelee Ounkham
//Last Modified on Decemeber 4, 2017
//
//*****************************************************************  

int main(){
//Setting up parameters
   int O1index,O2index,PBC1,PBC2,PBC3,PBC4,PBC5;
   float dist,ang,xcoord,ycoord,zcoord;
   int label,atomID,snap,lines,xlines,xlines2;
   int orow,xrow,y,z;
   int origGraphGeod[150][2];
   int xyzinteger[150][1];
//   int xyzinteger2[150][1];
   float xyzfloat[150][3];
   float xyzfloat2[150][3];
   char xyzchar[150][1];
   FILE *input,*input2,*xyzfile,*output,*output2;
   char ifile[100],ifile2[100],xyz[100],ofile[100],ofile2[100];


   for(snap=1; snap<=50001; snap++){
      sprintf(ifile,"OCO-HCl.input.O%d.xyz.O%d.xyz.GraphGeod",snap,snap);
      sprintf(xyz,"num-O%d.xyz",snap);
      sprintf(ofile,"Zpair-OO-%d.xyz",snap);
      input=fopen(ifile,"r");
      xyzfile=fopen(xyz,"r");
      output=fopen(ofile,"a");

//Initialize parameters
      for(y=0; y<=150; y++){
         for(z=0; z<4; z++){
            origGraphGeod[y][z]=0;
	    xyzinteger[y][z]=0;
	    xyzfloat[y][z]=0.0;
         }
      }
      lines=0;
      xlines=0;

//Store O indices from Z-GraphGeod files   
      while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&O1index,&O2index,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&ang)==9){
         origGraphGeod[lines][0]=O1index;
         origGraphGeod[lines][1]=O2index;
         lines++;
      }//end while-loop
      fclose(input);

//for(z=1; z<=lines; z++){
//printf("%d %d %d\n",origGraphGeod[z][0],origGraphGeod[z][1],lines);
//}

      while(fscanf(xyzfile,"%d %s %f %f %f\n",&label,&atomID,&xcoord,&ycoord,&zcoord)==5){ 
         xyzinteger[xlines][0]=label;
	 xyzfloat[xlines][0]=xcoord;
	 xyzfloat[xlines][1]=ycoord;
	 xyzfloat[xlines][2]=zcoord;
         xlines++;
      }//end while-loop
      fclose(xyzfile);

      for(orow=0; orow<lines; orow++){
	 for(xrow=0; xrow<xlines; xrow++){
	    if(origGraphGeod[orow][0]==xyzinteger[xrow][0]){//if Zundel O1 matches the label in xyz file
//printf("coln1: %d\n",origGraphGeod[orow][0]);
	       if(xyzfloat[xrow][0]<0.0){// if x-coord is negative
                  if(xyzinteger[xrow][0]<10){
	             fprintf(output,"%d    Z    %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
		  }
		  else if(xyzinteger[xrow][0]>=10 && xyzinteger[xrow][0]<100){
		     fprintf(output,"%d   Z    %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);

		  }
		  else if(xyzinteger[xrow][0]>=100){
		     fprintf(output,"%d  Z    %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
		  }
	       }
	       else if(xyzfloat[xrow][0]>0.0){
		  if(xyzinteger[xrow][0]<10){
	             fprintf(output,"%d    Z     %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
		  }
		  else if(xyzinteger[xrow][0]>=10 && xyzinteger[xrow][0]<100){
	             fprintf(output,"%d   Z     %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
		  }
		  else if(xyzinteger[xrow][0]>=100){
	             fprintf(output,"%d  Z     %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
                  }
	       }
	       xyzinteger[xrow][0]=0;
            }
            else if(origGraphGeod[orow][1]==xyzinteger[xrow][0]){ //if Zundel O2 matches the label in xyz file
//printf("coln2: %d %d\n",origGraphGeod[orow][1],orow);
	       if(xyzfloat[xrow][0]<0.0){
		  if(xyzinteger[xrow][0]<10){
	             fprintf(output,"%d    Z    %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
		  }
		  else if(xyzinteger[xrow][0]>=10 && xyzinteger[xrow][0]<100){
	             fprintf(output,"%d   Z    %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
		  }
		  else if(xyzinteger[xrow][0]>=100){
	             fprintf(output,"%d  Z    %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
                  }
	       }
	       else if(xyzfloat[xrow][0]>0.0){
		  if(xyzinteger[xrow][0]<10){
	             fprintf(output,"%d    Z     %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
		  }
		  else if(xyzinteger[xrow][0]>=10 && xyzinteger[xrow][0]<100){
	             fprintf(output,"%d   Z     %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
		  }
		  else if(xyzinteger[xrow][0]>=100){
	             fprintf(output,"%d  Z     %f %f %f %d\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2],
		     orow);
                  }
	       }	
    	       xyzinteger[xrow][0]=0;
	    }
         }//end xrow-loop
      }//end orow-loop

//Print the non-participating O's
      for(xrow=0; xrow<xlines; xrow++){
         if(xyzinteger[xrow][0]!=0){
//if O doesn't match any Zundel O indice, leave unchanged
   	    if(xyzfloat[xrow][0]<0.0){
		 if(xyzinteger[xrow][0]<10){
	             fprintf(output,"%d    O    %f %f %f 77\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2]);
		  }
		  else if(xyzinteger[xrow][0]>=10 && xyzinteger[xrow][0]<100){
	             fprintf(output,"%d   O    %f %f %f 77\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2]);
		  }
		  else if(xyzinteger[xrow][0]>=100){
	             fprintf(output,"%d  O    %f %f %f 77\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2]);
                  }
	    }
	    else if(xyzfloat[xrow][0]>0.0){
		  if(xyzinteger[xrow][0]<10){
	             fprintf(output,"%d    O     %f %f %f 77\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2]);
		  }
		  else if(xyzinteger[xrow][0]>=10 && xyzinteger[xrow][0]<100){
	             fprintf(output,"%d   O     %f %f %f 77\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2]);
		  }
		  else if(xyzinteger[xrow][0]>=100){
	             fprintf(output,"%d  O     %f %f %f 77\n",
	             xyzinteger[xrow][0],
	             xyzfloat[xrow][0],
                     xyzfloat[xrow][1],
	             xyzfloat[xrow][2]);
                  } 
	      }
	   }
      }//end xrow-loop
      fclose(output);
   }//end snap-loop

/*
    for(snap=1; snap<=1; snap++){
      sprintf(ifile2,"temp-%d.xyz",snap);
      sprintf(ofile2,"zundel-%d.xyz",snap);
      input2=fopen(ifile2,"r");
      output2=fopen(ofile2,"a");

//Initialize parameters
      for(y=0; y<=150; y++){
         for(z=0; z<3; z++){
	    xyzinteger2[y][z]=0;
	    xyzfloat2[y][z]=0.0;
         }
      }

      xlines2=0;

      while(fscanf(input2,"%d %s %f %f %f\n",&label,&atomID,&xcoord,&ycoord,&zcoord)==5){
         if(xcoord<0.0){ //if xcoord is negative
	    fprintf(output2,"%s    %f %f %f\n", 
	    atomID,
	    xcoord,
	    ycoord,
	    zcoord);
	 }
         else if(xcoord>0.0){
	    fprintf(output2,"%s     %f %f %f\n", 
	    atomID,
	    xcoord,
	    ycoord,
	    zcoord);
         } 
         xlines2++;
      }//end while-loop
      fclose(input2);
      fclose(output2);
    }//end snap-loop
*/
}//end main

