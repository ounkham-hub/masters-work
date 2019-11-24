#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//***************************************************
//Program created to calculate the O-O distance
//for previously identified O's participating in
//Zundel. 
//
//Input: labeled-Z-#.xyz
//Output: total-OO-distance.txt
//
//Made by Lelee Ounkham
//Last modified December 12, 2017
//***************************************************

int main(){
   int snap,atomID,pairno,lines,zlines;
   int y,z,row,row2,zrow,coln;
   float xcoord,ycoord,zcoord,distO1O2;
   float O1x,O1y,O1z,O2x,O2y,O2z;
   char type;
   int origGraphGeod[350][3];
   int zundelindex[350][3];
   char origtype[350];
   float origxyz[350][3];
   float zundelxyz[350][2][3];
   FILE *input,*output;
   char ifile[100],ofile[100];

   for(snap=1; snap<=50001; snap++){ 
      sprintf(ifile,"Zpairs-OO-%d.xyz",snap);
      input=fopen(ifile,"r");
      output=fopen("Zpair-distance.txt","a");

      lines=0;
      zlines=0;

      for(y=0; y<350; y++){
	 for(z=0; z<4; z++){
	    origGraphGeod[y][z]=0;
	    origxyz[y][z]=0.0;
	    zundelindex[y][z]=0;
	    zundelxyz[y][0][z]=0.0;
	    zundelxyz[y][1][z]=0.0;
	 }//end z-loop
      }//end y-loop

      while(fscanf(input,"%d %c %f %f %f %d\n",&atomID,&type,&xcoord,&ycoord,&zcoord,&pairno)==6){
	 origGraphGeod[lines][0]=snap;
         origGraphGeod[lines][1]=atomID;
	 origGraphGeod[lines][2]=pairno;
	 origtype[lines]=type;
	 origxyz[lines][0]=xcoord;
	 origxyz[lines][1]=ycoord;
	 origxyz[lines][2]=zcoord;
	 lines++;
       }
       //end while-loop
       fclose(input);
     
/*     for(z=0; z<lines; z++){ 
	printf("%d %d %d %c %f %f %f\n",
	origGraphGeod[z][0],
        origGraphGeod[z][1],
	origGraphGeod[z][2],
	origtype[z],
	origxyz[z][0],
	origxyz[z][1],
	origxyz[z][2]);
     }
*/

//Create Zundel array
       for(row=0; row<lines-1; row++){
	  for(row2=1; row2<lines; row2++){
	     if(origGraphGeod[row][2]!=77 && origGraphGeod[row2][2]!=77){
		if(origxyz[row][0]!=origxyz[row2][0]){
	           if(origGraphGeod[row][2]==origGraphGeod[row2][2]){//If Zundel OO pairs identified
		      zundelindex[zlines][0]=origGraphGeod[row][0]; //snapno
	              zundelindex[zlines][1]=origGraphGeod[row][1]; //O1 index
		      zundelindex[zlines][2]=origGraphGeod[row2][1]; //O2 index
	  	      zundelxyz[zlines][0][0]=origxyz[row][0]; //O1 x-coord
		      zundelxyz[zlines][0][1]=origxyz[row][1]; //O1 y-coord
		      zundelxyz[zlines][0][2]=origxyz[row][2]; //O1 z-coord
		      zundelxyz[zlines][1][0]=origxyz[row2][0]; //O2 x-coord
		      zundelxyz[zlines][1][1]=origxyz[row2][1]; //O2 y-coord
		      zundelxyz[zlines][1][2]=origxyz[row2][2]; //O2 z-coord
		      zlines++;

	              origGraphGeod[row][2]=77; //O1 index
		      origGraphGeod[row2][2]=77; //O2 index
		   }
		}
	     }//
	  }//end row2-loop
       }//end row-loop

/*       for(z=0; z<zlines; z++){
          printf("%d %d %f %f %f %f %f %f\n",
	  zundelindex[z][0],
	  zundelindex[z][1],
	  zundelxyz[z][0][0],
	  zundelxyz[z][0][1],
	  zundelxyz[z][0][2],
	  zundelxyz[z][1][0],
	  zundelxyz[z][1][1],
	  zundelxyz[z][1][2]);
       }
*/
	   O1x=0.0;
	   O1y=0.0;
	   O1z=0.0;
	   O2z=0.0;
	   O2y=0.0;
	   O2z=0.0;

 
        for(zrow=0; zrow<zlines; zrow++){
           O1x=zundelxyz[zrow][0][0];
 	   O1y=zundelxyz[zrow][0][1];
           O1z=zundelxyz[zrow][0][2];
	   O2x=zundelxyz[zrow][1][0];
	   O2y=zundelxyz[zrow][1][1];
	   O2z=zundelxyz[zrow][1][2];

           distO1O2 = sqrt( pow(O2x-O1x,2) + pow(O2y-O1y,2) + pow(O2z-O1z,2) );

           fprintf(output,"%d %d %d %f\n",
	   zundelindex[zrow][0],
	   zundelindex[zrow][1],
	   zundelindex[zrow][2],
	   distO1O2);

        }//end zrow-loop
       fclose(output);

   }//end snap-loop

}//end main-loop
