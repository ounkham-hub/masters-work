#include<stdio.h>
#include<stdlib.h>

//***********************************************************************
//Program created to read file containing all the list of H atoms
//that remain unreactive for more than 10k snapshots and output
//the weights associated with that particular bond.
//
//Input: centroid-all-H-PT-events.txt
//       EZW-O#.H#.wGraphGeod
//Output: PTs-weights-using-#window.txt
//
//Created by Lelee Ounkham
//Last modified May 23, 2018
//
//***********************************************************************

int main(){

  int Hindex,Oindex,flag,initialsnap,finalsnap;
  int H1,H2,H3,maxsnap,elines,plines,y,z,snap;
  int erow,ecoln,prow,Oxy,Opartner,PTsnap;
  int nlines,olines;
  float weight1,weight2,weight3;
  FILE *input,*PTinput,*output,*noutput;
  char ifile[100]; 

  int PTmatrix[10000][6]; 
  int ezwmatrix[103][6];
  int newOatom[25][3];
  int origOatom[25][3];
  float weightmatrix[103][6];

//Intialize counters 
   plines=0;
   maxsnap=61555;
   elines=0;
   nlines=0;
   olines=0;

//Read in list of unreactive H atoms and store into matrix

   PTinput=fopen("centroid-all-H-PT-snaps.txt","r");
   output=fopen("PTs-orig-Oatom-weights-5-redo.txt","a");
   noutput=fopen("PTs-new-Oatom-weights-5-redo.txt","a");
   while(fscanf(PTinput,"%d %d %d %d\n",&Hindex,&Oindex,&Opartner,&PTsnap)==4){
      PTmatrix[plines][0]=Hindex;
      PTmatrix[plines][1]=Oindex;
      PTmatrix[plines][2]=Opartner;
      PTmatrix[plines][3]=PTsnap;
      PTmatrix[plines][4]=PTsnap-5;
      PTmatrix[plines][5]=PTsnap+5;
      plines++;
   }
   fclose(PTinput);

/*
   for(z=0; z<plines; z++){
      printf("%d %d %d %d %d %d\n",
      PTmatrix[z][0],
      PTmatrix[z][1],
      PTmatrix[z][2],
      PTmatrix[z][3],
      PTmatrix[z][4],
      PTmatrix[z][5]);
   }//end z-loop
*/

for(prow=0; prow<plines; prow++){
   for(snap=1; snap<=maxsnap; snap++){
      sprintf(ifile,"EZW.O%d.H%d.wGraphGeod",snap,snap);
      input=fopen(ifile,"r");

      while(fscanf(input,"%d %d %d %d %d %f %f %f\n",&flag,&Oxy,&H1,&H2,&H3,&weight1,&weight2,&weight3)==8){
         ezwmatrix[elines][0]=flag;
	 ezwmatrix[elines][1]=Oxy;
	 ezwmatrix[elines][2]=H1;
	 ezwmatrix[elines][3]=H2;
	 ezwmatrix[elines][4]=H3;
	 ezwmatrix[elines][5]=0;

	 weightmatrix[elines][2]=weight1;
	 weightmatrix[elines][3]=weight2;
	 weightmatrix[elines][4]=weight3; 
         elines++;
      }
      fclose(input); 
 
/*
      for(y=0; y<elines; y++){
	 printf("%d %d %d %d %d %0.3f %0.3f %0.3f\n",
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
      for(erow=0; erow<elines; erow++){
         for(ecoln=2; ecoln<5; ecoln++){
	    if(snap >= PTmatrix[prow][4] && snap <= PTmatrix[prow][5]){ //if snapshot is within time window
//If O and H index matches list and GraphGeod
/*		  printf("%d %d %d %d %d %d %0.3f %0.3f %0.3f\n",
		  snap,
		  ezwmatrix[erow][0], //flag
		  ezwmatrix[erow][1], //Oindex
		  ezwmatrix[erow][2], //H1
		  ezwmatrix[erow][3], //H2
		  ezwmatrix[erow][4], //H3
		  weightmatrix[erow][2], //weight 1
		  weightmatrix[erow][3], //weight 2
		  weightmatrix[erow][4]); //weight 3
*/
//***********************************************************************
//Section to calculate store weights for each O atom per snap index
//participating in a PT event.
//***********************************************************************
		  if(PTmatrix[prow][1]==ezwmatrix[erow][1]){ //if O atom matches then it's the new O atom
	             if(PTmatrix[prow][0]==ezwmatrix[erow][ecoln]){ //If H index matches
/*
		      printf("%d %d %d %d %d %0.3f\n",
		      nlines,
		      snap,
		      ezwmatrix[erow][0],
		      ezwmatrix[erow][1],
		      ezwmatrix[erow][ecoln],
		      weightmatrix[erow][ecoln]);	
*/
		      fprintf(noutput,"%0.3f ",weightmatrix[erow][ecoln]);
		      ezwmatrix[erow][5]=1; //H index has been printed already

		    }
		    else if(ezwmatrix[erow][ecoln]==0 && ezwmatrix[erow][5]==0){ //If H index isn't bonded to O anymore, unmarked
		       fprintf(noutput,"0.000 "); //print 0.000 weight
		    }
		  }
		  else if(PTmatrix[prow][2]==ezwmatrix[erow][1]){ //If O atom doesnt match then it's the original O atom
	             if(PTmatrix[prow][0]==ezwmatrix[erow][ecoln]){ //If H index matches

/*
		      printf("%d %d %d %d %d %0.3f\n",
		      olines,
		      snap,
		      ezwmatrix[erow][0],
		      ezwmatrix[erow][1],
		      ezwmatrix[erow][ecoln],
		      weightmatrix[erow][ecoln]);
*/
		        fprintf(output,"%0.3f ",weightmatrix[erow][ecoln]);
		        ezwmatrix[erow][5]=1; //H index has been printed already
			}
		        else if(ezwmatrix[erow][ecoln]==0 && ezwmatrix[erow][5]==0){ //If H index is 0, structure is no longer Eigen or Zundel
		           fprintf(output,"0.000 "); //print 0.000 wewight
		       }
	           }
	       }
	       else{
	          ezwmatrix[erow][0]=0; //flag
		  ezwmatrix[erow][1]=0; //Oindex
		  ezwmatrix[erow][2]=0; //H1
		  ezwmatrix[erow][3]=0; //H2
		  ezwmatrix[erow][4]=0; //H3
		  weightmatrix[erow][2]=0.0; //weight 1
		  weightmatrix[erow][3]=0.0; //weight 2
		  weightmatrix[erow][4]=0.0; //weight 3
	          nlines=0;
   	          olines=0;	
	       }
	       ezwmatrix[erow][5]=0; //reset trigger
         }//end ecoln-loop
      }//end erow-loop

//Initialize ezwmatrix and weightmatrix
      for(y=0; y<elines; y++){
         for(z=0; z<6; z++){
	    ezwmatrix[y][z]=0;
	    weightmatrix[y][z]=0.0;
         }//end z-loop
      }//end y-loop
      elines=0; 
    
      for(y=0; y<olines; y++){
	for(z=0; z<3; z++){
	   origOatom[y][z]=0;
	}//end z-loop
      }//end y-loop
    
     for(y=0; y<nlines; y++){
	for(z=0; z<3; z++){
	   newOatom[y][z]=0;
	}//end z-loop
      }//end y-loop

    }//end snap-loop
    fprintf(output,"\n");
    fprintf(noutput,"\n");
  }//end prow-loop
  
   fclose(output);
   fclose(noutput);
}//end main-loop
