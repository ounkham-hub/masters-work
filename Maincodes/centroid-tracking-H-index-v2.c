#include<stdio.h>
#include<stdlib.h>

//*****************************************************************************
//Program created to read in GraphGeod files and track the coordination
//number of each H atom with respective O indices. The final output
//should yield information on peristence.
//
//Input: R-HCl.input.solB#.xyz.solC#.xyz.GraphGeod
//Output:
//
//Created by Lelee Ounkham
//Last modified on January 11, 2019
//
//*****************************************************************************

int main(){

   int snap,Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5; 
   int countHindex,y,z,Oindex,Hhold1,clines;
   int crow,mrow,indexno,maxsnap;
   float dist,angle; 
   FILE *input,*input1,*ZEoutput,*EZoutput;
   char ifile[100],ifile1[100];
  
   int currentsnap[500][7];
   int mastermatrix[500][7];
 
//Intialize matrices and counters
   for(y=0; y<500; y++){
      for(z=0; z<7; z++){
         currentsnap[y][z]=0;
	 mastermatrix[y][z]=0;
      }
   }//end z-loop
   countHindex=1;
   clines=0;
   indexno=103;
   maxsnap=61555;
  
   EZoutput=fopen("centroid-Eigen-to-Zundel.txt","a");
   ZEoutput=fopen("centroid-Zundel-to-Eigen.txt","a"); 
   for(snap=1; snap<=maxsnap; snap++){

/*      if(snap==maxsnap){ 

         for(y=1; y<213; y++){
            if(mastermatrix[y][0]!=0){
	       printf("%d %d %d %d %d %d\n",
	       mastermatrix[y][0],
	       mastermatrix[y][1],
	       mastermatrix[y][2],
	       mastermatrix[y][3],
	       mastermatrix[y][4],
               mastermatrix[y][5]);
           }
         }//end y-loop

         exit(0);
      }
*/
//*****************************************************************************
//Section I: Create mastermatrix which will encompass connectivity information
//per H atom..
//*****************************************************************************
       if(snap==1){
          sprintf(ifile,"H-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
          input=fopen(ifile,"r");

	  for(mrow=1; mrow<213; mrow++){
	     mastermatrix[mrow][0]=indexno;
	     indexno++;
	  }//end mrow-loop
 
/*          for(y=0; y<indexno; y++){
	     if(mastermatrix[y][0]!=0){
	        printf("Main: %d %d %d %d %d %d\n",
	        mastermatrix[y][0],
	        mastermatrix[y][1],
	        mastermatrix[y][2],
	        mastermatrix[y][3],
	        mastermatrix[y][4],
		mastermatrix[y][5]);
	     }
          }
*/
//*****************************************************************************
//Section II: Read in first snapshot and identify the coordination number for
//per H atom. 1 O partner, water or Eigen (indistinguishable) or 2 O for 
//Zundels.
//*****************************************************************************
          while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
             if(countHindex==1){
                Oindex=Oxy;
                Hhold1=Hyd;

                currentsnap[clines][0]=Hyd; //H index
                currentsnap[clines][1]=1; //CN
                currentsnap[clines][2]=Oxy;
                currentsnap[clines][3]=0; //Only has 1 O partner, it's a water
                currentsnap[clines][4]=snap;
                clines++;

             }
             else if(countHindex==2){
                if(Hyd==Hhold1){ //when H index matches
                   currentsnap[clines-1][0]=0;
                   currentsnap[clines-1][1]=0;
                   currentsnap[clines-1][2]=0;
                   currentsnap[clines-1][3]=0;
                   currentsnap[clines-1][4]=0;

                   currentsnap[clines][0]=Hyd; //indicative of Zundel
                   currentsnap[clines][1]=2; //CN
                   currentsnap[clines][2]=Oindex; //1st O partner
                   currentsnap[clines][3]=Oxy; //2nd O partner
                   currentsnap[clines][4]=snap;
                   clines++;
                   countHindex=1; //restart loop
                }
                else{
                   Oindex=Oxy;
                   Hhold1=Hyd;

                   currentsnap[clines][0]=Hyd; //H index
                   currentsnap[clines][1]=1; //CN
                   currentsnap[clines][2]=Oxy;
                   currentsnap[clines][3]=0; //Only has 1 O partner, it's a water
                   currentsnap[clines][4]=snap;
                   clines++;
                   countHindex=1; //restart loop
                }
             }
             countHindex++;
          }//end while-loop
          fclose(input);
/*
	  for(y=0; y<clines; y++){
	     if(currentsnap[y][0]!=0){
	        printf("%d %d %d %d %d\n",
	        currentsnap[y][0],
	        currentsnap[y][1],
	        currentsnap[y][2],
	        currentsnap[y][3],
	        currentsnap[y][4]);
	     }
	  }
*/
//*****************************************************************************
//Section III: Identify CN of H atoms in the first snapshot which will be used 
//to determine if there was a change in connectivity. Transfer info to 
//mastermatrix.
//*****************************************************************************
          for(crow=0; crow<clines; crow++){
	     for(mrow=1; mrow<213; mrow++){
		if(currentsnap[crow][0]!=0){
	           if(mastermatrix[mrow][0]==currentsnap[crow][0]){
//If H indices matches
		      if(currentsnap[crow][1]==1){
//If Eigen or water is identified
		         mastermatrix[mrow][1]=currentsnap[crow][1];//CN
		         mastermatrix[mrow][2]=currentsnap[crow][2]; //O1
		         mastermatrix[mrow][4]=snap; //initial snap
		         mastermatrix[mrow][5]=snap; //final snap

			 currentsnap[crow][0]=0;
		      }
		      else if(currentsnap[crow][1]==2){
//If Zundel is identified
		         mastermatrix[mrow][1]=currentsnap[crow][1];//CN
		         mastermatrix[mrow][2]=currentsnap[crow][2]; //O1
			 mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		         mastermatrix[mrow][4]=snap; //initial snap
		         mastermatrix[mrow][5]=snap; //final snap
			 
			 currentsnap[crow][0]=0;
		      }	      
		   }//end if H-index matches
		}
	     }//end crow-loop
	  }//end mrow-loop

//Delete currentsnap matrix for next snapshot
        for(y=0; y<clines; y++){
           for(z=0; z<7; z++){
              currentsnap[y][z]=0;
           }
        }//end z-loop
	clines=0;
	countHindex=1;
/*
          for(y=1; y<indexno; y++){
	     if(mastermatrix[y][0]!=0){
	        printf("%d %d %d %d %d %d\n",
	        mastermatrix[y][0],
	        mastermatrix[y][1],
	        mastermatrix[y][2],
	        mastermatrix[y][3],
	        mastermatrix[y][4],
		mastermatrix[y][5]);
	     }
          }//end y-loop
*/
       }//end if snap==1
//*****************************************************************************
//Section II: Read in proceeding GraphGeod file to identify the coordination # 
//for each H atom.
//*****************************************************************************

       else if(snap>1){

          sprintf(ifile1,"H-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
          input1=fopen(ifile1,"r");

          while(fscanf(input1,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
             if(countHindex==1){
                Oindex=Oxy;
                Hhold1=Hyd;

                currentsnap[clines][0]=Hyd; //H index
                currentsnap[clines][1]=1; //CN
                currentsnap[clines][2]=Oxy;
                currentsnap[clines][3]=0; //Only has 1 O partner, it's a water
                currentsnap[clines][4]=snap;
                clines++;

             }
             else if(countHindex==2){
                if(Hyd==Hhold1){ //when H index matches
                   currentsnap[clines-1][0]=0;
                   currentsnap[clines-1][1]=0;
                   currentsnap[clines-1][2]=0;
                   currentsnap[clines-1][3]=0;
                   currentsnap[clines-1][4]=0;

                   currentsnap[clines][0]=Hyd; //indicative of Zundel
                   currentsnap[clines][1]=2; //CN
                   currentsnap[clines][2]=Oindex; //1st O partner
                   currentsnap[clines][3]=Oxy; //2nd O partner
                   currentsnap[clines][4]=snap;
                   clines++;
                   countHindex=1; //restart loop
                }
                else{
                   Oindex=Oxy;
                   Hhold1=Hyd;

                   currentsnap[clines][0]=Hyd; //H index
                   currentsnap[clines][1]=1; //CN
                   currentsnap[clines][2]=Oxy;
                   currentsnap[clines][3]=0; //Only has 1 O partner, it's a water
                   currentsnap[clines][4]=snap;
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
	       currentsnap[y][0],
	       currentsnap[y][1],
	       currentsnap[y][2],
	       currentsnap[y][3],
	       currentsnap[y][4]);
	    }
         }
*/
//*****************************************************************************
//Section IV: Section compare connectivities in initial snapshots to proceeding
//ones. 
//*****************************************************************************
	for(crow=0; crow<clines; crow++){
	   for(mrow=1; mrow<213; mrow++){
	      if(currentsnap[crow][0]!=0){
	         if(mastermatrix[mrow][0]==currentsnap[crow][0]){
//If H index matches
		    if(mastermatrix[mrow][1]==currentsnap[crow][1]){
//Case 1: If structure does not change,same O partners
		       if(mastermatrix[mrow][2]==currentsnap[crow][2] && mastermatrix[mrow][3]==currentsnap[crow][3]){
/*			  printf("H index is unchanged\n");

 			  printf("%d %d %d %d %d %d\n",
			  mastermatrix[mrow][0], //H index
			  mastermatrix[mrow][1], //CN
			  mastermatrix[mrow][2], //O1	
			  mastermatrix[mrow][3], //O2
			  mastermatrix[mrow][4], //initial snap
			  mastermatrix[mrow][5]); //final snap
*/
			  mastermatrix[mrow][5]=snap;
			  mastermatrix[mrow][6]=1; //turn flag on 

			  currentsnap[crow][0]=0;
		        }
		        else if(mastermatrix[mrow][3]!=currentsnap[crow][3] && mastermatrix[mrow][3]==currentsnap[crow][2]){
//If O1 is different, update mastermatrix
			   printf("1a: Case when O1 is different, O2 is same!\n");

 			   printf("%d %d %d %d %d %d %d\n",
			   snap,
			   mastermatrix[mrow][0], //H index
			   mastermatrix[mrow][1], //CN
			   mastermatrix[mrow][2], //O1	
			   mastermatrix[mrow][3], //O2
			   mastermatrix[mrow][4], //initial snap
			   mastermatrix[mrow][5]); //final snap

			   mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
			   mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
			   mastermatrix[mrow][4]=snap; //update initial snap
			   mastermatrix[mrow][5]=snap; //update final snap
			   mastermatrix[mrow][6]=1; //turn flag on
			   currentsnap[crow][0]=0;
		        }
		        else if(mastermatrix[mrow][3]!=currentsnap[crow][3]){
//If O2 is different, update mastermatrix
			   printf("1b. Case when O2 is different, O1 is same!\n");

 			   printf("%d %d %d %d %d %d %d\n",
			   snap,
			   mastermatrix[mrow][0], //H index
			   mastermatrix[mrow][1], //flag
			   mastermatrix[mrow][2], //O1	
			   mastermatrix[mrow][3], //O2
			   mastermatrix[mrow][4], //initial snap
			   mastermatrix[mrow][5]); //final snap

			   mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
			   mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
			   mastermatrix[mrow][4]=snap; //update initial snap
			   mastermatrix[mrow][5]=snap; //update final snap
			   mastermatrix[mrow][6]=1; //turn flag on
			   currentsnap[crow][0]=0;
		       }
		    }//end of Case 1 situations
		    else if(mastermatrix[mrow][1]==1 && currentsnap[crow][1]==2){
//Case 2: If structure changes from Eigen/water to Zundel, update mastermatrix
//		       printf("2: Eigen/water transformed to Zundel\n");

 		       fprintf(EZoutput,"%d %d %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][1], //CN
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][3], //O2
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
		      
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;

		    }//end Case 2 situation
		    else if(mastermatrix[mrow][1]==2 && currentsnap[crow][1]==1){
//Case 3: If structure changes from Zundel to Eigen/water, update mastermatrix
//		       printf("3: Zundel transformed to Eigen/water\n");

 		       fprintf(ZEoutput,"%d %d %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][1], //CN
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][3], //O2
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
	  
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;
		    }//end Case 3 situation
		    else if(mastermatrix[mrow][1]==0 && currentsnap[crow][1]!=0){
//Case 4: Occurs when connectivity reappears in proceeding snapshot, update matrix	

		       mastermatrix[mrow][1]=currentsnap[crow][1];
		       mastermatrix[mrow][2]=currentsnap[crow][2];
		       mastermatrix[mrow][3]=currentsnap[crow][3];
		       mastermatrix[mrow][4]=snap;
		       mastermatrix[mrow][5]=snap;

		       printf("5: Connectivity reappeared!\n");
 		       printf("%d %d %d %d %d %d %d\n",
		       snap,
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][1], //CN
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][3], //O2
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap

		       mastermatrix[mrow][6]=1;
		       currentsnap[crow][0]=0;
		    }//end Case 4 situation
		 }//if non-zero element
	      }   
	   }//end mrow-loop
	}//end crow-loop

//*****************************************************************************
//Section V: Connectivities that "disappear" in a snapshot should also
//be accounted for. These would be considered structural changes.
//*****************************************************************************

	    for(mrow=1; mrow<213; mrow++){
	       if(mastermatrix[mrow][6]==0 && mastermatrix[mrow][0]!=0){
//If H index matches, update matrix

	          printf("6: Connectivity disappeared!\n");
 		  printf("%d %d %d %d %d %d %d\n",
		  snap,
		  mastermatrix[mrow][0], //Hindex
		  mastermatrix[mrow][1], //CN
		  mastermatrix[mrow][2], //O1	
		  mastermatrix[mrow][3], //O2
		  mastermatrix[mrow][4], //initial snap
		  mastermatrix[mrow][5]); //final snap

//Reset elements for H index 
		  mastermatrix[mrow][1]=0;
		  mastermatrix[mrow][2]=0;	
		  mastermatrix[mrow][3]=0;
		  mastermatrix[mrow][4]=0;
		  mastermatrix[mrow][5]=0;
	       }//end crow-loop
	    }//end mrow-loop

//Initialize currentsnap matrix for the snapshot
        for(y=0; y<clines; y++){
           for(z=0; z<7; z++){
              currentsnap[y][z]=0;
           }
        }//end z-loop
	clines=0;

//Initialize mastermatrix "switches" for next set of snapshots
        for(y=0; y<213; y++){
	   mastermatrix[y][6]=0;
        }//end z-loop
     }//end else-if snap>1

     if(snap==maxsnap){
        for(y=1; y<213; y++){
           if(mastermatrix[y][0]!=0){
	      printf("%d %d %d %d %d %d\n",
	      mastermatrix[y][0],
	      mastermatrix[y][1],
	      mastermatrix[y][2],
	      mastermatrix[y][3],
	      mastermatrix[y][4],
              mastermatrix[y][5]);
	   }
        }//end y-loop
     }//end snap=-maxsnap

   }//end snap-loop 
   fclose(EZoutput);
   fclose(ZEoutput);
}//end main-loop 
