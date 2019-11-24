#include<stdio.h>
#include<stdlib.h>

//*****************************************************************************
//Program created to identify weights of unreactive eigens and zundels at
//every snapshot.
//
//Input: HCl.covalent.O#.xyz.H#.xyz.wGraphGeod
//       UNR-EZ.O#.xyz.H#.xyz.index
//
//Created by Lelee Ounkham
//Last Modified March 3, 2019
//
//******************************************************************************

int main(){

  int Oxy,Hyd,index,crow,row,snap;
  int y,z,clines,lines,O,H,flag;
  float weight,dist;
  FILE *cinput,*winput,*eoutput,*zoutput,*spoutput;
  char cfile[100],wfile[100];

  int centmatrix[480][4];
  int indexmatrix[480][4];
  float weightmatrix[480][1];


//Section to intialize matrices and counters

    for(y=0; y<480; y++){
       for(z=0; z<4; z++){
	   centmatrix[y][z]=0;
	   indexmatrix[y][z]=0;
	   weightmatrix[y][0]=0.000;
       }
    }
    clines=0;
    lines=0;


    eoutput=fopen("complete-UNR-CBs-in-Eigens.txt","a");
    spoutput=fopen("complete-UNR-CBs-in-SP.txt","a");
    zoutput=fopen("complete-UNR-CBs-in-Zundel.txt","a");
    for(snap=1; snap<=61555; snap++){
       sprintf(wfile,"HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap,snap); //open current snap
       winput=fopen(wfile,"r");

       while(fscanf(winput,"%d %d %f %f\n",&O,&H,&dist,&weight)==4){

           indexmatrix[lines][0]=snap;
	   indexmatrix[lines][1]=O;
	   indexmatrix[lines][2]=H;
	   weightmatrix[lines][0]=weight;
	   lines++;

       }//end while-loop
       fclose(winput);
/*
       for(y=0; y<lines; y++){
          printf("%d %d %d %0.3f\n",
	  indexmatrix[y][0],
	  indexmatrix[y][1],
	  indexmatrix[y][2],
	  weightmatrix[y][0]);
       }
*/
       sprintf(cfile,"UNR-EZ.O%d.xyz.H%d.xyz.index",snap,snap); //open current snap
       cinput=fopen(cfile,"r");

//Store indices of unreactive O and H indices from centroid data into centmatrix
       while(fscanf(cinput,"%d %d %d %d\n",&index,&flag,&Oxy,&Hyd)==4){

           centmatrix[clines][0]=index;
	   centmatrix[clines][1]=Oxy;
	   centmatrix[clines][2]=Hyd;
	   centmatrix[clines][3]=flag;
	   clines++;

       }//end while-loop
       fclose(cinput);
/*
       for(y=0; y<clines; y++){
	  printf("%d %d %d %d\n",
	  centmatrix[y][0],
	  centmatrix[y][1],
	  centmatrix[y][2],
	  centmatrix[y][3]);
       }
*/
//******************************************************************************
//Section to identify matching O and H indices then output weights into a 
//corresponding Zundel or Eigen CBs file.
//******************************************************************************

       for(crow=0; crow<clines; crow++){
          for(row=0; row<lines; row++){
	     if(centmatrix[crow][0]==indexmatrix[row][0]){ //if snapshots matches
		if(centmatrix[crow][1]==indexmatrix[row][1]){ //if O index matches
		   if(centmatrix[crow][2]==indexmatrix[row][2]){ //if H index matches
//Same bond has been identified, output weights
			if(centmatrix[crow][3]==2){ //CBs in special pair
			   fprintf(spoutput,"%d %d %d %d %0.3f\n",
			   centmatrix[crow][0], //snapshot
			   centmatrix[crow][3], //flag
			   centmatrix[crow][1], //Oindex
			   centmatrix[crow][2], //H index
			   weightmatrix[row][0]); //OH weight
			}
			else if(centmatrix[crow][3]==1){ //CBs in eigens
			   fprintf(eoutput,"%d %d %d %d %0.3f\n",
			   centmatrix[crow][0], //snapshot
			   centmatrix[crow][3], //flag
			   centmatrix[crow][1], //Oindex
			   centmatrix[crow][2], //H index
			   weightmatrix[row][0]); //OH weight
			}
			else if(centmatrix[crow][3]==4){ //outer bonds in Zundel
			   fprintf(zoutput,"%d %d %d %d %0.3f\n",
			   centmatrix[crow][0], //snapshot
			   centmatrix[crow][3], //flag
			   centmatrix[crow][1], //Oindex
			   centmatrix[crow][2], //H index
			   weightmatrix[row][0]); //OH weight
		       }
		   }
		}
	     }
	  }//end row-loop
       }//end crow-loop


//Intialize matrices and counters for next snapshot
       for(y=0; y<480; y++){
          for(z=0; z<4; z++){
	     centmatrix[y][z]=0;
	     indexmatrix[y][z]=0;
	     weightmatrix[y][0]=0.000;
         }
      }
      clines=0;
      lines=0;
   }//end snap-loop
   fclose(eoutput);
   fclose(zoutput);
   fclose(spoutput);
}//end main-loop
