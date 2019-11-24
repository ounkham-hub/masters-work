#include<stdio.h>
#include<stdlib.h>

//***************************************************************************
//Program created to track the persistence of weights across the the
//entire trajectory.
//
//Input: H-HCl.covalent.O#.xyz.H#.xyz.wGraphGeod
//
//Created by Lelee Ounkham
//Last modified January 31, 2019
//
//***************************************************************************

int main(){
   int Oxy,Hyd,snap,crow,Oindex,Hhold1,clines,y,z,countHindex,maxsnap;
   int indexno,mrow;
   float dist,weight,whold1;
   FILE *input,*input2,*output;
   char ifile[100],ifile2[100];

   int currentmatrix[300][4];
   int mastermatrix[300][9];
   float currentweight[300][4];
   float masterweight[300][4];

//Initialize matrix and counters
   for(y=0; y<300; y++){
      for(z=0; z<4; z++){
	 currentmatrix[y][z]=0;
	 currentweight[y][z]=0.000;
      }
      for(z=0; z<9; z++){
	 mastermatrix[y][z]=0;
	 masterweight[y][z]=0.000;
      }//end z-loop
   }//end y-loop
   clines=0;
   countHindex=1;
   indexno=103;
   maxsnap=61555;
  
   output=fopen("total-connectivity-of-WN-NOweights.txt","a"); 
   for(snap=1; snap<=maxsnap; snap++){
      if(snap==1){
         sprintf(ifile,"H-HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap,snap); //open current snap
         input=fopen(ifile,"r"); 

         while(fscanf(input,"%d %d %f %f\n",&Oxy,&Hyd,&dist,&weight)==4){
            if(countHindex==1){
                Oindex=Oxy;
                Hhold1=Hyd;
		whold1=weight;

                currentmatrix[clines][0]=Hyd; //H index
                currentmatrix[clines][1]=Oxy;
                currentmatrix[clines][2]=0; //Only has 1 O partner, it's a water
                currentmatrix[clines][3]=1; //CN

		currentweight[clines][1]=weight;	
                clines++;

             }
             else if(countHindex==2){
                if(Hyd==Hhold1){ //when H index matches
                   currentmatrix[clines-1][0]=0;
                   currentmatrix[clines-1][1]=0;
                   currentmatrix[clines-1][2]=0;
                   currentmatrix[clines-1][3]=0;
		   currentweight[clines-1][1]=0.000;

                   currentmatrix[clines][0]=Hyd; //indicative of Zundel
                   currentmatrix[clines][1]=Oindex; //O1
                   currentmatrix[clines][2]=Oxy; //O2
                   currentmatrix[clines][3]=2; //CN

		   currentweight[clines][1]=whold1;
		   currentweight[clines][2]=weight;
                   clines++;

                   countHindex=1; //restart loop
                }
                else{
                   Oindex=Oxy;
                   Hhold1=Hyd;
		   whold1=weight;

                   currentmatrix[clines][0]=Hyd; //H index
                   currentmatrix[clines][1]=Oxy;
                   currentmatrix[clines][2]=0; //Only has 1 O partner, it's a water
                   currentmatrix[clines][3]=1;

		   currentweight[clines][1]=weight;
                   clines++;
		   countHindex=1;
                }
             }
             countHindex++;
          }//end while-loop
          fclose(input);

/*	 for(crow=0; crow<clines; crow++){
	    if(currentmatrix[crow][0]!=0){
               printf("%d %d %d %d %0.3f %0.3f %d\n",
	       snap,
	       currentmatrix[crow][0],
	       currentmatrix[crow][1],
	       currentmatrix[crow][2],
	       currentweight[crow][1],
	       currentweight[crow][2],
	       currentmatrix[crow][3]);
	    }
	 }//end crow-loop
*/
//Allocate H indices in mastermatrix 
         for(mrow=1; mrow<213; mrow++){
            mastermatrix[mrow][0]=indexno;
            indexno++;
         }//end mrow-loop

//        for(y=1; y<213; y++){
//            printf("%d\n",mastermatrix[y][0]);
//         }//end y-loop

//***************************************************************************
//Section to update mastematrix which contains persistence of each 
//covalent bond and associated weight. 
//***************************************************************************
        for(crow=0; crow<clines; crow++){
	   for(mrow=1; mrow<213; mrow++){
	      if(currentmatrix[crow][0]!=0){ //for nonzero elements
		 if(mastermatrix[mrow][0]==currentmatrix[crow][0]){
//If H indices matches, update info
		    if(currentmatrix[crow][3]==1){
//If Eigen/water. H bonded to single O
		       mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
		       mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
		       mastermatrix[mrow][3]=snap;
		       mastermatrix[mrow][4]=snap;
		       mastermatrix[mrow][7]=currentmatrix[crow][3]; //CN

		       masterweight[mrow][1]=currentweight[crow][1]; //O1H
		       masterweight[mrow][2]=currentweight[crow][2]; //O2H
		    }
		    else if(currentmatrix[crow][3]==2){
//If H has two O partners
		       mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
		       mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
		       mastermatrix[mrow][3]=snap; 
		       mastermatrix[mrow][4]=snap;
		       mastermatrix[mrow][5]=snap;
		       mastermatrix[mrow][6]=snap;
		       mastermatrix[mrow][7]=currentmatrix[crow][3]; //CN

		       masterweight[mrow][1]=currentweight[crow][1]; //O1H
		       masterweight[mrow][2]=currentweight[crow][2]; //O2H
		    }
		 }
	      }
	   }
	} 
/*
	for(mrow=1; mrow<213; mrow++){
	   printf("%d %d %d %d %d %d %d %d %0.3f %0.3f\n",
	   mastermatrix[mrow][0],
	   mastermatrix[mrow][1],
	   mastermatrix[mrow][2],
	   mastermatrix[mrow][3],
	   mastermatrix[mrow][4],
	   mastermatrix[mrow][5],
	   mastermatrix[mrow][6],
	   mastermatrix[mrow][7],
	   masterweight[mrow][1],
	   masterweight[mrow][2]);
	}
*/

         for(y=0; y<clines; y++){
	    for(z=0; z<4; z++){
	       currentmatrix[y][z]=0;
	       currentweight[y][z]=0.000;
	    }
         }//end y-loop
         clines=0;
         countHindex=1;

      }//end snap==1 loop
  
//***************************************************************************
//Section to transfer connectivity information from proceeding snapshot to
//currentmatrix.
//***************************************************************************
      else if(snap>1){
         sprintf(ifile2,"H-HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap,snap); //open current snap
         input2=fopen(ifile2,"r"); 

         while(fscanf(input2,"%d %d %f %f\n",&Oxy,&Hyd,&dist,&weight)==4){
            if(countHindex==1){
                Oindex=Oxy;
                Hhold1=Hyd;
		whold1=weight;

                currentmatrix[clines][0]=Hyd; //H index
                currentmatrix[clines][1]=Oxy;
                currentmatrix[clines][2]=0; //Only has 1 O partner, it's a water
                currentmatrix[clines][3]=1; //CN

		currentweight[clines][1]=weight;	
                clines++;

             }
             else if(countHindex==2){
                if(Hyd==Hhold1){ //when H index matches
                   currentmatrix[clines-1][0]=0;
                   currentmatrix[clines-1][1]=0;
                   currentmatrix[clines-1][2]=0;
                   currentmatrix[clines-1][3]=0;
		   currentweight[clines-1][1]=0.000;

                   currentmatrix[clines][0]=Hyd; //indicative of Zundel
                   currentmatrix[clines][1]=Oindex; //O1
                   currentmatrix[clines][2]=Oxy; //O2
                   currentmatrix[clines][3]=2; //CN

		   currentweight[clines][1]=whold1;
		   currentweight[clines][2]=weight;
                   clines++;

                   countHindex=1; //restart loop
                }
                else{
                   Oindex=Oxy;
                   Hhold1=Hyd;
		   whold1=weight;

                   currentmatrix[clines][0]=Hyd; //H index
                   currentmatrix[clines][1]=Oxy;
                   currentmatrix[clines][2]=0; //Only has 1 O partner, it's a water
                   currentmatrix[clines][3]=1;

		   currentweight[clines][1]=weight;
                   clines++;
		   countHindex=1;
                }
             }
             countHindex++;
          }//end while-loop
          fclose(input2);
/*
	 for(crow=0; crow<clines; crow++){
	    if(currentmatrix[crow][0]!=0){
               printf("%d %d %d %d %0.3f %0.3f %d\n",
	       snap,
	       currentmatrix[crow][0],
	       currentmatrix[crow][1],
	       currentmatrix[crow][2],
	       currentweight[crow][1],
	       currentweight[crow][2],
	       currentmatrix[crow][3]);
	    }
	 }//end crow-loop
*/

//***************************************************************************
//Section to compare connectivity patterns in consecutive snapshots.
//***************************************************************************
	 for(crow=0; crow<clines; crow++){
	    for(mrow=1; mrow<213; mrow++){
	       if(mastermatrix[mrow][0]==currentmatrix[crow][0]){
//If H indices matches
		  if(mastermatrix[mrow][7]==1 && currentmatrix[crow][3]==1){
//Case 1a: When H of Eigen/Water remains unchanged 
		     if(mastermatrix[mrow][1]==currentmatrix[crow][1]){
//O atoms matches
		        mastermatrix[mrow][4]=snap; //update initial snap
			mastermatrix[mrow][8]=1; //switch
		     }//end Case 1a
		     else if(mastermatrix[mrow][1]==currentmatrix[crow][1]){
//Case 1b: If H matches, but O partners do not - tunneling occurred

			printf("Tunneling occurred %d\n",snap);

		        mastermatrix[mrow][1]=currentmatrix[crow][1]; //update O index
			mastermatrix[mrow][3]=snap; //restart initial and final snap
			mastermatrix[mrow][4]=snap;
			mastermatrix[mrow][8]=1;
		     }//end case 1b
		  }//end of case 1 - H associated with Eigen/Waters
		  else if(mastermatrix[mrow][7]==2 && currentmatrix[crow][3]==2){
		     if(mastermatrix[mrow][1]==currentmatrix[crow][1] && mastermatrix[mrow][2]==currentmatrix[crow][2]){
//Case 2a: When Zundel remains unchanged
			mastermatrix[mrow][4]=snap; //update final snap for O1H
			mastermatrix[mrow][6]=snap; //update final snap for O2H
			mastermatrix[mrow][8]=1;
		     }
		     else if(mastermatrix[mrow][1]==currentmatrix[crow][2] && mastermatrix[mrow][2]==currentmatrix[crow][1]){
//Case 2b: When Zundel remains unchanged, but O partners swap - result of how files are read in
			mastermatrix[mrow][4]=snap; //update final snap for O1H
			mastermatrix[mrow][6]=snap; //update final snap for O2H
			mastermatrix[mrow][8]=1;
		     }
		  }//end of Case 2 - unchanged Zundels
		  else if(mastermatrix[mrow][7]==2 && currentmatrix[crow][3]==1){
//Case 3a: When Zundel dissociates into Eigen, O1 matches
		     if(mastermatrix[mrow][1]==currentmatrix[crow][1]){
                        printf("3a: Change from Zundel to Eigen/water! %d\n",snap);

                        fprintf(output,"%d %d %d %d\n",
                        mastermatrix[mrow][0],//H index
                        mastermatrix[mrow][2],//O2 - bond dissociated
                        mastermatrix[mrow][5], //initial snap
                        mastermatrix[mrow][6]); //final snap

                        mastermatrix[mrow][4]=snap; //update final snap 
                        mastermatrix[mrow][7]=1; //CN
                        mastermatrix[mrow][8]=1; //switch

                        mastermatrix[mrow][2]=0; //remove partner O
                        mastermatrix[mrow][6]=0;
                        mastermatrix[mrow][7]=0;
		     }
		     else if(mastermatrix[mrow][2]==currentmatrix[crow][1]){
//Case 3b: When Zundel dissociates into EIgen, O2 matches
                        printf("3b: Change from Zundel to Eigen/water! %d\n",snap);

                        fprintf(output,"%d %d %d %d\n",
                        mastermatrix[mrow][0],//H index
                        mastermatrix[mrow][1],//O1 - bond dissociated
                        mastermatrix[mrow][3], //initial snap
                        mastermatrix[mrow][4]); //final snap

                        mastermatrix[mrow][1]=mastermatrix[mrow][2]; //copy over to coln 1
                        mastermatrix[mrow][3]=mastermatrix[mrow][5]; //copy initial snap
                        mastermatrix[mrow][4]=snap; //update snap
                        mastermatrix[mrow][7]=1;
                        mastermatrix[mrow][8]=1;

                        mastermatrix[mrow][2]=0; //remove partner O
                        mastermatrix[mrow][5]=0;
                        mastermatrix[mrow][6]=0;
		     }
		  }//end of Zundel to Eigen cases
		  else if(mastermatrix[mrow][7]==1 && currentmatrix[crow][3]==2){
//Case 4a: When Eigen/Water forms Zundel
		     if(mastermatrix[mrow][1]==currentmatrix[crow][1]){
//When O matches O1 in mastermatrix
                        printf("4a: Eigen/water to Zundel! %d\n",snap);
                        printf("%d %d %d %d\n",
                        currentmatrix[crow][0],//H index
                        currentmatrix[crow][1],//O1
                        currentmatrix[crow][2], //initial snap
                        currentmatrix[crow][3]); //flag

                        mastermatrix[mrow][2]=currentmatrix[crow][2];
                        mastermatrix[mrow][4]=snap; //update final snap, bond didn't break
                        mastermatrix[mrow][5]=snap; //start new counter for new bond
                        mastermatrix[mrow][6]=snap;

                        mastermatrix[mrow][7]=2;
                        mastermatrix[mrow][8]=1;
		     }
		     else if(mastermatrix[mrow][1]==currentmatrix[crow][2]){
//When O matches O2 in mastermatrix
                        printf("4b: Eigen/water to Zundel! %d\n",snap);
                        printf("%d %d %d %d\n",
                        currentmatrix[crow][0],//H index
                        currentmatrix[crow][1],//O1
                        currentmatrix[crow][2], //initial snap
                        currentmatrix[crow][3]); //flag

			mastermatrix[mrow][1]=currentmatrix[crow][2];
                        mastermatrix[mrow][2]=currentmatrix[crow][1];
                        mastermatrix[mrow][4]=snap; //update final snap, bond didn't break
                        mastermatrix[mrow][5]=snap; //start new counter for new bond
                        mastermatrix[mrow][6]=snap;

                        mastermatrix[mrow][7]=2;
                        mastermatrix[mrow][8]=1;
		     }
		  }//end of case 4 situations
		  else if(mastermatrix[mrow][7]==0){
                     if(currentmatrix[crow][3]==1){
                        printf("Reappearing bond! %d %d\n",currentmatrix[crow][0], snap);
                        mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
                        mastermatrix[mrow][2]=0; //O2, not Zundel
                        mastermatrix[mrow][3]=snap; //initial snap
                        mastermatrix[mrow][4]=snap;
                        mastermatrix[mrow][5]=0;
                        mastermatrix[mrow][6]=0;
                        mastermatrix[mrow][7]=currentmatrix[crow][3];
                        mastermatrix[mrow][8]=1;
                     }
                     else if(currentmatrix[crow][3]==2){
                        printf("Reappearing bond! %d %d\n",currentmatrix[crow][0], snap);
                        mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
                        mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
                        mastermatrix[mrow][3]=snap; //initial snap
                        mastermatrix[mrow][4]=snap;
                        mastermatrix[mrow][5]=snap;
                        mastermatrix[mrow][6]=snap;
                        mastermatrix[mrow][7]=currentmatrix[crow][3]; //flag
                        mastermatrix[mrow][8]=1;
		     }
		  }
	       }
	    }//end mrow-loop
	 }//end crow-loop
	 
//***************************************************************************
//Section to update mastermatrix when bond disappears from a graph. 
//***************************************************************************
        for(mrow=1; mrow<213; mrow++){
           if(mastermatrix[mrow][8]==0){ //for non-updated connectivities
              if(mastermatrix[mrow][7]==1){
                 printf("Connectivity lost %d: %d %d %d %d\n",
                 snap,
                 mastermatrix[mrow][0],
                 mastermatrix[mrow][1],
                 mastermatrix[mrow][3],
                 mastermatrix[mrow][4]);

                 mastermatrix[mrow][1]=0; //bond disappeared
                 mastermatrix[mrow][3]=0;
                 mastermatrix[mrow][4]=0;
                 mastermatrix[mrow][7]=0;
                 mastermatrix[mrow][8]=1;
              }
              else if(mastermatrix[mrow][5]==2){
                 printf("Connectivity lost %d: %d %d %d %d %d %d %d %d\n",
                 snap,
                 mastermatrix[mrow][0],
                 mastermatrix[mrow][1],
                 mastermatrix[mrow][2],
                 mastermatrix[mrow][3],
                 mastermatrix[mrow][4],
                 mastermatrix[mrow][5],
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

//Print out mastermatrix at end of trajectory
         if(snap==maxsnap){
            for(mrow=1; mrow<213; mrow++){
               if(mastermatrix[mrow][7]==1){
                  fprintf(output,"%d %d %d %d\n",
                  mastermatrix[mrow][0], //H index
                  mastermatrix[mrow][1], //O1 
                  mastermatrix[mrow][3], //initial snap
                  mastermatrix[mrow][4]);
               }
               else if(mastermatrix[mrow][7]==2){
                  fprintf(output,"%d %d %d %d\n",
                  mastermatrix[mrow][0], //H index
                  mastermatrix[mrow][1], //O1 
                  mastermatrix[mrow][3], //initial snap
                  mastermatrix[mrow][4]);

                  fprintf(output,"%d %d %d %d\n",
                  mastermatrix[mrow][0], //H index
                  mastermatrix[mrow][2], //O1 
                  mastermatrix[mrow][5], //initial snap
                  mastermatrix[mrow][6]);
               }
            }
         }

//***************************************************************************
//Section to intialize currentmatrix for the proceeding snapshot.
//***************************************************************************

      for(y=0; y<clines; y++){
	 for(z=0; z<4; z++){
	    currentmatrix[y][z]=0;
	    currentweight[y][z]=0.000;
	 }
      }//end y-loop
      clines=0;
      countHindex=1;
      
      for(mrow=1; mrow<213; mrow++){
	 mastermatrix[mrow][8]=0;
      }

      }//end snap>1 loop
   }//end snap-loop
   fclose(output);

}//end main-l(op
