#include<stdio.h>
#include<stdlib.h>

//*************************************************************************************
//Program created to track the persistence of individual covalent bonds in centroid
//datases. This will be from the H's perspective and dependent on the change in 
//O parter irrespective of its molecular type (i.e. Eigen, Zundel, Water).
//
//Input: H-HCl.input.solB3.xyz.solC#.xyz.GraphGeod
//
//Created by Lelee Ounkham
//Last modified on January 22, 2019
//
//*************************************************************************************

int main(){

   int snap,Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5;
   int mrow,y,z,indexno,clines,maxsnap,crow;
   int Hhold1,Oindex,countHindex;
   float dist,angle;
   FILE *input,*output,*input1,*output1,*output2;
   char ifile[100],ifile1[100],ofile[100];

   int mastermatrix[300][9]; 
   int currentsnap[300][9];

//Initlaize matrices and counters
   for(y=0; y<300; y++){
      for(z=0; z<9; z++){
	 mastermatrix[y][z]=0;
	 currentsnap[y][z]=0;
      }
   }
   indexno=103;
   clines=0;
   countHindex=1;
   maxsnap=61555;

   output=fopen("total-O-partner-list-noMS.txt","a");
//   output1=fopen("tunneling-events.txt","a");
   for(snap=1; snap<=maxsnap; snap++){
      if(snap==1){
         sprintf(ifile,"H-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
         input=fopen(ifile,"r");

         while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
	    if(countHindex==1){
                Oindex=Oxy;
                Hhold1=Hyd;

                currentsnap[clines][0]=Hyd; //H index
                currentsnap[clines][1]=Oxy;
                currentsnap[clines][2]=0; //Only has 1 O partner, it's a water
                currentsnap[clines][5]=1; //CN
                clines++;

             }
             else if(countHindex==2){
                if(Hyd==Hhold1){ //when H index matches
                   currentsnap[clines-1][0]=0;
                   currentsnap[clines-1][1]=0;
                   currentsnap[clines-1][2]=0;
                   currentsnap[clines-1][5]=0;

                   currentsnap[clines][0]=Hyd; //indicative of Zundel
                   currentsnap[clines][1]=Oindex; //O1
                   currentsnap[clines][2]=Oxy; //O2
                   currentsnap[clines][5]=2; //CN
                   clines++;
                   countHindex=1; //restart loop
                }
                else{
                   Oindex=Oxy;
                   Hhold1=Hyd;

                   currentsnap[clines][0]=Hyd; //H index
                   currentsnap[clines][1]=Oxy;
                   currentsnap[clines][2]=0; //Only has 1 O partner, it's a water
                   currentsnap[clines][5]=1;
                   clines++;
                   countHindex=1; //restart loop
                }
             }
             countHindex++;
         }
	 fclose(input);
/*
	 for(y=0; y<clines; y++){
	    if(currentsnap[y][0]!=0){
	       printf("%d %d %d %d %d\n",
	       snap,
	       currentsnap[y][0],
	       currentsnap[y][1],
	       currentsnap[y][2],
	       currentsnap[y][5]);
	    }
	 }//end y-loop
*/
//Allocate H indices in mastermatrix 
         for(mrow=1; mrow<213; mrow++){
            mastermatrix[mrow][0]=indexno;
            indexno++;
         }//end mrow-loop

/*	 for(y=1; y<213; y++){
	    printf("%d %d\n",
	    mastermatrix[y][0],
	    mastermatrix[y][1]);
	 }//end crow-loop
*/

//*************************************************************************************
//Section I: Populate mastermatrix with initial H indices and corresponding O partners,
//
//*************************************************************************************

         for(crow=0; crow<clines; crow++){
            for(mrow=1; mrow<213; mrow++){
	       if(mastermatrix[mrow][0]==currentsnap[crow][0]){
//If H indices matches
		  mastermatrix[mrow][1]=currentsnap[crow][1]; //O1
		  mastermatrix[mrow][2]=currentsnap[crow][2]; //O2
		  mastermatrix[mrow][3]=snap; //initial snap	
		  mastermatrix[mrow][4]=snap; //final snap
		  mastermatrix[mrow][5]=currentsnap[crow][5]; //CN
		  mastermatrix[mrow][6]=snap; //for second bond in zundel
		  mastermatrix[mrow][7]=snap; //for second bond in zundel
	       }
	    }//end mrow-loop
	 }//end crow-loop
/*
	 for(y=1; y<213; y++){
	    printf("%d %d %d %d %d %d\n",
	    mastermatrix[y][0],
	    mastermatrix[y][1],
            mastermatrix[y][2],
	    mastermatrix[y][3],
	    mastermatrix[y][4],
	    mastermatrix[y][5]);
	 }//end y-loop
*/

//for(y=1; y<213; y++){
//fprintf(output,"%d %d %d\n",mastermatrix[y][5],mastermatrix[y][1],mastermatrix[y][2]);


//Initialize currentsnap for next snapshot
         for(y=0; y<clines; y++){
            for(z=0; z<9; z++){
	       currentsnap[y][z]=0;
            }
         }
	 clines=0;
	 countHindex=1;

      }//end snap==1

//Reading next round of snapshots
      else if(snap>1){
         sprintf(ifile1,"H-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
         input1=fopen(ifile1,"r");

         while(fscanf(input1,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
	    if(countHindex==1){
                Oindex=Oxy;
                Hhold1=Hyd;

                currentsnap[clines][0]=Hyd; //H index
                currentsnap[clines][1]=Oxy;
                currentsnap[clines][2]=0; //Only has 1 O partner, it's a water
                currentsnap[clines][5]=1; //CN
                clines++;

             }
             else if(countHindex==2){
                if(Hyd==Hhold1){ //when H index matches
                   currentsnap[clines-1][0]=0;
                   currentsnap[clines-1][1]=0;
                   currentsnap[clines-1][2]=0;
                   currentsnap[clines-1][5]=0;

                   currentsnap[clines][0]=Hyd; //indicative of Zundel
                   currentsnap[clines][1]=Oindex; //O1
                   currentsnap[clines][2]=Oxy; //O2
                   currentsnap[clines][5]=2; //CN
                   clines++;
                   countHindex=1; //restart loop
                }
                else{
                   Oindex=Oxy;
                   Hhold1=Hyd;

                   currentsnap[clines][0]=Hyd; //H index
                   currentsnap[clines][1]=Oxy;
                   currentsnap[clines][2]=0; //Only has 1 O partner, it's a water
                   currentsnap[clines][5]=1;
                   clines++;
                   countHindex=1; //restart loop
                }
             }
             countHindex++;
         }//end while-loop
	 fclose(input1);

/*
	 for(y=0; y<clines; y++){
	    if(currentsnap[y][0]!=0){
	       printf("%d %d %d %d %d\n",
  	       snap,
	       currentsnap[y][0],
	       currentsnap[y][1],
               currentsnap[y][2],
	       currentsnap[y][5]);
	    }
	 }//end y-loop
*/
//*************************************************************************************
//Section II: Compare O partners per H index every consecutive snapshot and update
//as necessary. The output will yield information about persistences.
//
//*************************************************************************************
         for(crow=0; crow<clines; crow++){
            for(mrow=1; mrow<213; mrow++){
	       if(currentsnap[crow][0]==mastermatrix[mrow][0]){
//If H atoms matches
	          if(mastermatrix[mrow][5]==1 && currentsnap[crow][5]==1){ //if flag matches
	             if(currentsnap[crow][1]==mastermatrix[mrow][1]){
//Case 1a: If O partner is unchanged for Eigen or water, structure is the same
//			   printf("1a: No change in Eigen/Water! %d\n",snap);

		        mastermatrix[mrow][4]=snap; //update final snapshot
		        mastermatrix[mrow][8]=1; 		    
		     }//end Case 1a
		     else if(currentsnap[crow][1]!=mastermatrix[mrow][1] && mastermatrix[mrow][1]!=0){
//Case 1b, tunneling: If O partner changes and not reappearing bond, update mastermatrix 
//		           printf("1b: Change in O partner, Eigen/Water! %d\n",snap);

/*		           fprintf(output1,"%d %d %d %d %d %d\n",
		           mastermatrix[mrow][0],//H index
		           mastermatrix[mrow][1],//O1
		           mastermatrix[mrow][2],//O2 - should be zero 
		           mastermatrix[mrow][3], //initial snap
			   mastermatrix[mrow][4], //final snap
			   mastermatrix[mrow][5]); //flag
*/
		        mastermatrix[mrow][1]=currentsnap[crow][1]; //update O index
		        mastermatrix[mrow][2]=currentsnap[crow][2]; //update O index, should be 0
		        mastermatrix[mrow][3]=snap; //restart final snap
			mastermatrix[mrow][4]=snap;
		        mastermatrix[mrow][8]=1;
	             }//end Case 1b
		  }//end Eigen/Water case
		  else if(mastermatrix[mrow][5]==2 && currentsnap[crow][5]==2){//For Zundels
		     if(currentsnap[crow][1]==mastermatrix[mrow][1] && currentsnap[crow][2]==mastermatrix[mrow][2]){
//Case 1c: When proceeding snapshot remains Zundel with the O1 and O2, update mastermatrix
//		           printf("1c: No change in Zundel! %d\n",snap);
/*
		           printf("%d %d %d %d %d %d\n",
		           mastermatrix[mrow][0],//H index
		           mastermatrix[mrow][1],//O1
		           mastermatrix[mrow][2],//O2 - should be nonzero
		           mastermatrix[mrow][3], //initial snap
			   mastermatrix[mrow][4], //final snap
			   mastermatrix[mrow][5]); //flag
*/
		        mastermatrix[mrow][4]=snap; //update final snap
		        mastermatrix[mrow][7]=snap; //update final snap for 2nd bond
		        mastermatrix[mrow][8]=1;
		     }//end Case 1c
		     else if(currentsnap[crow][1]==mastermatrix[mrow][2] && currentsnap[crow][2]==mastermatrix[mrow][1]){
//Case 1c: When proceeding snapshot remains Zundel with the O1 and O2, update mastermatrix
//		           printf("1d: No change in Zundel! %d\n",snap);
/*
		           printf("%d %d %d %d %d %d\n",
		           mastermatrix[mrow][0],//H index
		           mastermatrix[mrow][1],//O1
		           mastermatrix[mrow][2],//O2 - should be nonzero
		           mastermatrix[mrow][3], //initial snap
			   mastermatrix[mrow][4], //final snap
			   mastermatrix[mrow][5]); //flag
*/
		        mastermatrix[mrow][4]=snap; //update final snap
		        mastermatrix[mrow][7]=snap; //update final snap for 2nd bond
		        mastermatrix[mrow][8]=1;
		     }//end Case 1d

		  }
		  else if(mastermatrix[mrow][5]==2 && currentsnap[crow][5]==1){
//Case 2: Zundel dissociates into Eigen or Water
		     if(mastermatrix[mrow][1]==currentsnap[crow][1]){
//When O matches O1 in mastermatrix, PT
	                printf("2a: Change from Zundel to Eigen/water! %d\n",snap);

		        fprintf(output,"%d %d %d %d\n",
		        mastermatrix[mrow][0],//H index
		        mastermatrix[mrow][2],//O2 - bond dissociated
		        mastermatrix[mrow][6], //initial snap
		        mastermatrix[mrow][7]); //final snap

		        mastermatrix[mrow][4]=snap; //update final snap	
	   	        mastermatrix[mrow][5]=1; //CN
		        mastermatrix[mrow][8]=1; //switch

			mastermatrix[mrow][2]=0; //remove partner O
			mastermatrix[mrow][6]=0;
		        mastermatrix[mrow][7]=0;
		     }
		     else if(mastermatrix[mrow][2]==currentsnap[crow][1]){
//When O matches O2 in mastermatrix, PT
	                printf("2b: Change from Zundel to Eigen/water! %d\n",snap);

		        fprintf(output,"%d %d %d %d\n",
		        mastermatrix[mrow][0],//H index
		        mastermatrix[mrow][1],//O1 - bond dissociated
		        mastermatrix[mrow][3], //initial snap
		        mastermatrix[mrow][4]); //final snap

		        mastermatrix[mrow][1]=mastermatrix[mrow][2]; //copy over to coln 1
		        mastermatrix[mrow][3]=mastermatrix[mrow][6]; //copy initial snap
			mastermatrix[mrow][4]=snap; //update snap
	   	        mastermatrix[mrow][5]=1;
		        mastermatrix[mrow][8]=1;

			mastermatrix[mrow][2]=0; //remove partner O
			mastermatrix[mrow][6]=0;
		        mastermatrix[mrow][7]=0;
		     }	
		  }//end of Case 2
		  else if(mastermatrix[mrow][5]==1 && currentsnap[crow][5]==2){
//Case 3: Eigen/water forms Zundel
		     if(mastermatrix[mrow][1]==currentsnap[crow][1]){
//When O matches O1 in mastermatrix, zundel formed

//	                printf("3: Eigen/water to Zundel! %d\n",snap);
/*
		        printf("%d %d %d %d %d\n",
		        mastermatrix[mrow][0],//H index
		        mastermatrix[mrow][1],//O1
		        mastermatrix[mrow][3], //initial snap
		        mastermatrix[mrow][4], //final snap
		        mastermatrix[mrow][5]); //flag
*/
		        mastermatrix[mrow][1]=currentsnap[crow][1];
		        mastermatrix[mrow][2]=currentsnap[crow][2];
		        mastermatrix[mrow][4]=snap; //update final snap, bond didn't break
			mastermatrix[mrow][6]=snap; //start new counter for new bond
			mastermatrix[mrow][7]=snap;

	   	        mastermatrix[mrow][5]=2;
		        mastermatrix[mrow][8]=1;	
		     }
		     else if(mastermatrix[mrow][1]==currentsnap[crow][2]){
//When O matches O2 in mastermatrix, zundel formed

//	                printf("3: Eigen/water to Zundel! %d\n",snap);
/*
		        printf("%d %d %d %d %d\n",
		        mastermatrix[mrow][0],//H index
		        mastermatrix[mrow][1],//O1
		        mastermatrix[mrow][3], //initial snap
		        mastermatrix[mrow][4], //final snap
		        mastermatrix[mrow][5]); //flag
*/
		        mastermatrix[mrow][2]=currentsnap[crow][1];
		        mastermatrix[mrow][4]=snap; //update final snap, bond didn't break
			mastermatrix[mrow][6]=snap; //start new counter for new bond
			mastermatrix[mrow][7]=snap;

	   	        mastermatrix[mrow][5]=2;
		        mastermatrix[mrow][8]=1;	
		     }	
		  } //end Case 3
		  else if(mastermatrix[mrow][5]==0){
//For reappearing bonds, update mastermatrix
		     if(currentsnap[crow][5]==1){
		        printf("Reappearing bond! %d %d\n",currentsnap[crow][0], snap);
	                mastermatrix[mrow][1]=currentsnap[crow][1]; //O1
		        mastermatrix[mrow][2]=0; //O2, not Zundel
		        mastermatrix[mrow][3]=snap; //initial snap
		        mastermatrix[mrow][4]=snap;
		        mastermatrix[mrow][5]=currentsnap[crow][5]; //flag
		        mastermatrix[mrow][6]=0;
		        mastermatrix[mrow][7]=0;
		        mastermatrix[mrow][8]=1;
		     }
		     else if(currentsnap[crow][5]==2){
		        printf("Reappearing bond! %d %d\n",currentsnap[crow][0], snap);
	                mastermatrix[mrow][1]=currentsnap[crow][1]; //O1
		        mastermatrix[mrow][2]=currentsnap[crow][2]; //O2
		        mastermatrix[mrow][3]=snap; //initial snap
		        mastermatrix[mrow][4]=snap;
		        mastermatrix[mrow][5]=currentsnap[crow][5]; //flag
		        mastermatrix[mrow][6]=snap;
		        mastermatrix[mrow][7]=snap;
		        mastermatrix[mrow][8]=1;
		     }
		  } 
	       }
	    }//end mrow-loop
         }//end crow-loop

//*************************************************************************************
//Section III: To update mastermatrix when a bond "disappears" from a graph.
//
//*************************************************************************************

	for(mrow=1; mrow<213; mrow++){
	   if(mastermatrix[mrow][8]==0){ //for non-updated connectivities
	      if(mastermatrix[mrow][5]==1){
	         printf("Connectivity lost %d: %d %d %d %d\n",
		 snap,
	         mastermatrix[mrow][0],
	         mastermatrix[mrow][1],
	         mastermatrix[mrow][3],
		 mastermatrix[mrow][4]);

	         mastermatrix[mrow][1]=0; //bond disappeared
	         mastermatrix[mrow][3]=0; 
	         mastermatrix[mrow][4]=0;
		 mastermatrix[mrow][5]=0;
	         mastermatrix[mrow][8]=1;
	      }
	      else if(mastermatrix[mrow][5]==2){
		 printf("Connectivity lost %d: %d %d %d %d %d %d %d\n",
		 snap,
	         mastermatrix[mrow][0],
	         mastermatrix[mrow][1],
		 mastermatrix[mrow][2],
	         mastermatrix[mrow][3],
		 mastermatrix[mrow][4],
		 mastermatrix[mrow][6],
		 mastermatrix[mrow][7]);

	         mastermatrix[mrow][1]=0; //bond disappeared
	         mastermatrix[mrow][2]=0; 
	         mastermatrix[mrow][3]=0; 
	         mastermatrix[mrow][4]=0;
		 mastermatrix[mrow][5]=0;
		 mastermatrix[mrow][6]=0;
		 mastermatrix[mrow][7]=0;
	         mastermatrix[mrow][8]=1;
	 
	      }
	   }
	}//end mrow-loop

//Print out connectivity # of O all H's per snapshot

//for(y=1; y<213; y++){
//fprintf(output2,"%d %d %d\n",mastermatrix[y][5],mastermatrix[y][1],mastermatrix[y][2]);

//Print out mastermatrix at during the at last snapshot
         if(snap==maxsnap){
            for(mrow=1; mrow<213; mrow++){
	       if(mastermatrix[mrow][5]==1){
	          fprintf(output,"%d %d %d %d %d\n",
	          mastermatrix[mrow][0], //H index
  	          mastermatrix[mrow][1], //O1 
 	          mastermatrix[mrow][3], //initial snap
 	          mastermatrix[mrow][4],
	          mastermatrix[mrow][5]); //final snap
	       }
	       else if(mastermatrix[mrow][5]==2){
	          fprintf(output,"%d %d %d %d %d %d %d %d\n",
	          mastermatrix[mrow][0], //H index
  	          mastermatrix[mrow][1], //O1 
		  mastermatrix[mrow][2], //O2
 	          mastermatrix[mrow][3], //initial snap
 	          mastermatrix[mrow][4],
		  mastermatrix[mrow][6],
		  mastermatrix[mrow][7],
	          mastermatrix[mrow][5]); //final snap
	       }
	    }
	 }

//Initialize currensnap for proceeding snapshot
	 for(y=1; y<213; y++){
	    mastermatrix[y][8]=0;
	 }

          for(y=0; y<clines; y++){
            for(z=0; z<9; z++){
	       currentsnap[y][z]=0;
            }
         }
	 clines=0;
	 countHindex=1;
      }//end if snap>1
   }//end snap-loop
   fclose(output);
//   fclose(output2);
}//end main-loop
