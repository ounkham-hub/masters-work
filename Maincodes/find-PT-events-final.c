#include<stdio.h>
#include<stdlib.h>

//*******************************************************************
//Program created to read O-H (0.0-1.3A) GraphGeod files, obtain 
//the list of O-partners and identify when partner switching occurs
//in a given snapshot.There are two outputs:
//   1) zundel-list, contains the list of all the Zundels formed 
//      in a simulation
//   2) PT-list, contains a list of paths of OCO with respect 
//   	to time 	
//
// *Note: All GraphGeod,master-list, and snapshot-i files must be 
// organized prior to analyis. Use sort -k1 -n command that organizes 
// files according to first(or which ever) column of file*
//
// *IMPORTANT! Change the number of total # snapshots or maxsnap
// to appropriately set up matrixes. 
//
//Created by Lelee Ounkham
//Last update: July 31, 2017
//
//******************************************************************

int main()
{

//*********************************************************************************************
//Section I: Identify OCO from covalent bond GraphGeods files 
//*********************************************************************************************

//Parameters to read O-H GraphGeods
	int Oxy,Hyd,int1,int2,int3,int4,int5;
	float float1,float2;
	char ifile[100],ofile[100];
	FILE *input,*output;

//Parameters for arrays
	int mxline=300;
	int Oindex,Hhold1,Hhold2,maxsnap;
	int OCOarray[mxline][5];	 

//Parameters to read major file containing all OCO-lists
	int snapno,OCO,CBH1,CBH2,CBH3;

//Loop and counter parameters
	int i,j,z,lines;

	maxsnap=3; //!!CHANGE ACCORDING TO SYSTEM!! The total number of snapshots!
/*
//Read in GraphGeod files
	for(i=1; i<maxsnap+1; i++){
	   sprintf(ifile,"R-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",i,i);
	   sprintf(ofile,"OCO-list-%d.txt",i); 
	   input=fopen(ifile,"r");
	   output=fopen(ofile,"w");

//Zero out counters
	   j=1;
	   lines=0;

//Identifying OCO from O-H GraphGeod files
 
	   while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&int1,&int2,&int3,&int4,&int5,&float1,&float2)==9){
              if(j==1){
                 Oindex=Oxy;
                 Hhold1=Hyd;		
              }
              else if(j==2){
                 if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;
                 }
                 else{
                    j=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
                 }
              }
              else if(j==3){ //if there is a 3rd occurrence then send to memory/array
                 if(Oindex==Oxy){
		    OCOarray[lines][0]=i;
                    OCOarray[lines][1]=Oxy;
		    OCOarray[lines][2]=Hhold1;
		    OCOarray[lines][3]=Hhold2;
	 	    OCOarray[lines][4]=Hyd;
		    lines++;
		    j=0; //restarts Oindex to read the next line
		 }
	         else{
	            j=1;
	            Oindex=Oxy;
		    Hhold1=Hyd;
	         }
	      }	
	      j++;
	   } //end while

           for(z=0; z<lines; z++){
              fprintf(output,"%d %d %d %d %d\n",
	      OCOarray[z][0],
	      OCOarray[z][1],
	      OCOarray[z][2],
	      OCOarray[z][3],
	      OCOarray[z][4]);
           }//end z-loop

	   fclose(input);
	   fclose(output);
	 
	} //end i-loop 
*/
//*********************************************************************************************
//Section II: Find Zundels and save set of OCO in first snapshot 
//*********************************************************************************************

//Parameters to read in OCO-list
	char ifile2[100],ifile3[100];
	char ofile2[100],ofile3[100];
	FILE *input2,*input3;
	FILE *output2,*output3;

//Parameters for arrays
	int startlist[50][7];
	int nextlist[50][7];
	int pathlist[maxsnap][20][5];

//Loop and counter parameters
	int c,d,x,y,srow,scoln,nrow,ncoln;
	int snap,slines,nlines,countPT,countZtoE;
	int countEE,countEZ,countZdiss,countZZ;

        for(c=0; c<50; c++){//For startlist
           for(d=0; d<7; d++){
              startlist[c][d]=0;
           }
        }

        for(c=0; c<50; c++){//for nextlist
           for(d=0; d<7; d++){
              nextlist[c][d]=0;
           }
        }

        for(c=0; c<maxsnap+1; c++){//for pathlist
           for(d=0; d<20; d++){
              pathlist[c][d][0]=0;
              pathlist[c][d][1]=0;
              pathlist[c][d][2]=0;
              pathlist[c][d][3]=0;
	      pathlist[c][d][4]=0;
           }
        }

//Zero out counters
	slines=1;
	nlines=1;
	countPT=0;
	countZtoE=0;
	countZdiss=0;
	countEE=0;
	countEZ=0;
	countZZ=0;

	output2=fopen("zundel-complete-list.txt","a");
	output3=fopen("proton-history.txt","a");

	for(snap=1; snap<maxsnap; snap++){
  	   sprintf(ifile3,"OCO-list-%d.txt",snap+1);
	   input3=fopen(ifile3,"r");

//Identify Zundels in the first snapshot and flag them, Zundel=1 and Eigen=2

 	   if(snap==1){  
	      sprintf(ifile2,"OCO-list-%d.txt",snap);
    	      input2=fopen(ifile2,"r");

	      while(fscanf(input2,"%d %d %d %d %d\n",&snapno,&OCO,&CBH1,&CBH2,&CBH3)==5){
	         startlist[slines][0]=snapno;
	         startlist[slines][1]=OCO;
	         startlist[slines][2]=CBH1;
	         startlist[slines][3]=CBH2;
	         startlist[slines][4]=CBH3;
	         slines++;
	      }//end while

 	      fclose(input2);

	      for(srow=1; srow<slines; srow++){
		 for(scoln=2; scoln<5; scoln++){
		    for(c=2; c<slines; c++){
		       for(d=2; d<5; d++){
			  if(startlist[srow][5]==0){ //If the OCO hasn't been flagged
		             if(startlist[srow][scoln]==startlist[c][d]){ //if H matches
			        if(startlist[srow][1]!=startlist[c][1]){ //if O's don't match
			           fprintf(output2,"%d %d %d %d %d\n",
			           snap, 
			           snap,
			           startlist[srow][1], //O1 in Zundel
			           startlist[c][1], //O2 in Zundel
			           startlist[srow][scoln]); //shared proton

				   startlist[srow][5]=1;
				   startlist[c][5]=1; //flag Zundel with 1
				   startlist[srow][6]=startlist[c][1];
				   startlist[c][6]=startlist[srow][1]; //store in partner info for PT check
				}
			     }
		          }
		       }//end d-loop
		    }//end c-loop
		 }//end scoln-loop
	      }//end srow-loop

	      for(srow=1; srow<slines; srow++){
		 if(startlist[srow][5]!=1){ //if OCO has NOT been flagged, it's an Eigen (2)
		    startlist[srow][5]=2;
		 }

//for(z=1; z<slines; z++){
//printf("%d %d %d %d %d %d\n",startlist[z][1],startlist[z][2],startlist[z][3],startlist[z][4],startlist[z][5],startlist[z][6]);
//}
		 
//Fill in pathlist for the first snapshot. Note: pathlist contains H info for each OCO for tracking purposes
  	         pathlist[0][0][0]=snap;
	         pathlist[0][srow][0]=startlist[srow][1]; //OCO index
	         pathlist[0][srow][1]=startlist[srow][2]; //CBH1
	         pathlist[0][srow][2]=startlist[srow][3]; //CBH2
	         pathlist[0][srow][3]=startlist[srow][4]; //CBH3
	         pathlist[0][srow][4]=startlist[srow][5]; //flag
	      } //end srow-loop 

	      for(z=0; z<slines; z++){		
		 if(pathlist[0][z][0]!=0){ //if element is nonzero
		    fprintf(output3,"%d ", pathlist[0][z][0]);
		 }
	      }//end z-loop
	      fprintf(output3,"\n");

//Read in snapshot 2 and compare
              while(fscanf(input3,"%d %d %d %d %d\n",&snapno,&OCO,&CBH1,&CBH2,&CBH3)==5)
              {
                 nextlist[nlines][0]=snapno;
                 nextlist[nlines][1]=OCO;
                 nextlist[nlines][2]=CBH1;
                 nextlist[nlines][3]=CBH2;
                 nextlist[nlines][4]=CBH3;
                 nlines++;
              }//end while

              fclose(input3);	 

	      for(nrow=1; nrow<nlines; nrow++){
		 for(ncoln=2; ncoln<5; ncoln++){
		    for(c=2; c<nlines; c++){
		       for(d=2; d<5; d++){
			  if(nextlist[nrow][5]==0){ //If the OCO hasn't been flagged
		             if(nextlist[nrow][ncoln]==nextlist[c][d]){ //if H matches
			        if(nextlist[nrow][1]!=nextlist[c][1]){ //if O's don't match
				   fprintf(output2,"%d %d %d %d %d\n",
			           snap, 
			           snap+1,
			           nextlist[nrow][1], //O1 in Zundel
			           nextlist[c][1], //O2 in Zundel
			           nextlist[nrow][ncoln]); //shared proton

				   nextlist[nrow][5]=1;
				   nextlist[c][5]=1; //flag Zundel with 1
				   nextlist[nrow][6]=nextlist[c][1];
				   nextlist[c][6]=nextlist[nrow][1]; //store partner list for check
				}
			     }
		          }
		       }//end d-loop
		    }//end c-loop
		 }//end ncoln-loop
	      }//end nrow-loop

	      for(nrow=1; nrow<nlines; nrow++){
		 if(nextlist[nrow][5]!=1){ //if OCO has NOT been flagged, it's an Eigen (2)
		    nextlist[nrow][5]=2;
		 }
	      }//end nrow-loop

//for(z=1; z<nlines; z++){
//printf("%d %d %d %d %d %d\n",nextlist[z][1],nextlist[z][2],nextlist[z][3],nextlist[z][4],nextlist[z][5],nextlist[z][6]);
//}

//Check for changes and update pathlist accordingly

	      for(srow=1; srow<slines; srow++){
	         for(nrow=1; nrow<nlines; nrow++){
		    if(startlist[srow][1]==nextlist[nrow][1]){ //if OCO match
		       if(startlist[srow][5]==nextlist[nrow][5]){ //if flag did NOT change, no transfer
		       }
		       else if(startlist[srow][5]!=nextlist[nrow][5]){ //if flag DID change, transfer!
//No need to double count, when a Zundel becomes an Eigen - it will be counted as a flag change
//printf("%d %d %d %d\n",startlist[srow][1],nextlist[nrow][1],startlist[srow][5],nextlist[nrow][5]);
		       }
		    }
		    else if(startlist[srow][1]!=nextlist[nrow][1]){ //if OCO don't match
		       for(c=2; c<5; c++){ //loop through all H's
		          for(d=2; d<5; d++){
			     if(startlist[srow][c]==nextlist[nrow][d]){ //if a H's match, transfer!
			        if(startlist[srow][5]!=1){ //EE transfer OR Eigen to Zundel (not a Zundel initially)
				   if(nextlist[nrow][5]==1){ //Eigen to Zundel
				      countEZ++;
				   }
				   else if(nextlist[nrow][5]==2){ //Eigen to Eigen
				      countEE++;
				   }
//printf("Type 1 event:\n");
//printf("%d %d %d %d %d\n",snap,startlist[srow][1],nextlist[nrow][1],startlist[srow][5],nextlist[nrow][5]);
			        }
			        else if(startlist[srow][5]==1 && nextlist[nrow][5]!=1){
				   for(z=1; z<nlines; z++){ //does the OCO in Zundel appear in the next snap?
				      if(startlist[srow][1]==nextlist[z][1]){ //if OCO matches partner doest match ANY OCO in next snap
					 countZdiss++; //Occurs only when Zundel dissociates to water and Eigen
//printf("trigger 1: %d %d\n",startlist[srow][1],nextlist[z][1]);
				      }
			              else if(startlist[srow][6]==nextlist[z][1]){
					 countZdiss++; //Occurs only when Zundel dissociates to water and Eigen
//printf("trigger 2: %d %d\n",startlist[srow][1],nextlist[z][1]);
				      }
				      else if(startlist[srow][1]!=nextlist[z][1] && startlist[srow][6]!=nextlist[z][1]){
					 for(x=2; x<5; x++){
					    for(y=2; y<5; y++){
					       if(startlist[srow][x]==nextlist[z][y]){ //if any H's match
			                          countZtoE++; //Zundel transfer to water
					       }
					    }//end y-loop
				         }//end x-loop
				      }//end else-if
				   }//end z-loop
			        }
				else if(startlist[srow][5]==1 && nextlist[nrow][5]==1){ //If both structures initially and after are Zundel
				   if(startlist[srow][1]==nextlist[nrow][1]){ //if OCO matches
				      if(startlist[srow][6]!=nextlist[nrow][6]){//if second Zundel OCO pair does not match
				         for(x=2; x<5; x++){
					    for(y=2; y<5; y++){
					       if(startlist[srow][x]==nextlist[nrow][y]){ //if a set of H's is identified
						   countZZ++;
					       }
					    }//end y-loop
					 }//end x-loop
				      }
				   }
				   else if(startlist[srow][1]!=nextlist[nrow][1]){ //if OCO's don't match
				      if(startlist[srow][6]==nextlist[nrow][6]){ //second Z-OCO partner does match
					 for(x=2; x<5; x++){
					    for(y=2; y<5; y++){
					       if(startlist[srow][x]==nextlist[nrow][y]){
						  countZZ++;
					       }
					    }//end y-loop
					 }//end x-loop
				      }
				   }
				}
			     }
			  }//end d-loop
		       }//end c-loop
		    }//else if
	         }//end nrow-loop
	      }//end srow-loop

	      for(c=1; c<nlines; c++){
		 pathlist[1][0][0]=nextlist[c][0];
		 pathlist[1][c][0]=nextlist[c][1];
		 pathlist[1][c][1]=nextlist[c][2];
		 pathlist[1][c][2]=nextlist[c][3];
		 pathlist[1][c][3]=nextlist[c][4];
		 pathlist[1][c][4]=nextlist[c][5];
	      }//end c-loop

	      for(z=0; z<nlines; z++){		
		 if(pathlist[1][z][0]!=0){ //if element is nonzero
		    fprintf(output3,"%d ", pathlist[1][z][0]);
		 }
	      }//end z-loop
	      fprintf(output3,"\n");

//Intialize startlist
	     slines=nlines;

             for(c=0; c<50; c++){//For startlist
                for(d=0; d<7; d++){
                   startlist[c][d]=0;
                }
             }

//Copy nextlist into startlist
	       for(d=1; d<slines; d++){
	       startlist[d][0]=nextlist[d][0];
	       startlist[d][1]=nextlist[d][1];
	       startlist[d][2]=nextlist[d][2];
	       startlist[d][3]=nextlist[d][3];
	       startlist[d][4]=nextlist[d][4];
	       startlist[d][5]=nextlist[d][5];
	       startlist[d][6]=nextlist[d][6];
	    }

//Intialize nextlist	
             for(c=0; c<50; c++){//for nextlist
                for(d=0; d<7; d++){
                   nextlist[c][d]=0;
                }
            }
	 }//end if snap==1

//***********************************************************************************
//Only for the first snapshot - do the above loop!
//***********************************************************************************

	 else if(snap!=1){ //if not the last snapshot store memory to compare arrays
	     nlines=1;
             while(fscanf(input3,"%d %d %d %d %d\n",&snapno,&OCO,&CBH1,&CBH2,&CBH3)==5)
              {
                 nextlist[nlines][0]=snapno;
                 nextlist[nlines][1]=OCO;
                 nextlist[nlines][2]=CBH1;
                 nextlist[nlines][3]=CBH2;
                 nextlist[nlines][4]=CBH3;
                 nlines++;
              }//end while

              fclose(input3);
//	   }//end else
//*********************************************************************************************
//Section III: Count PT events,identify Zundels and update pathlist
//*********************************************************************************************

//Identify Zundels in nextlist then compare to startlist to track transfers.

	      for(nrow=1; nrow<nlines; nrow++){
		 for(ncoln=1; ncoln<5; ncoln++){
		    for(c=2; c<nlines; c++){
		       for(d=2; d<5; d++){
			  if(nextlist[nrow][5]==0){ //If the OCO hasn't been flagged
		             if(nextlist[nrow][ncoln]==nextlist[c][d]){ //if H matches
			        if(nextlist[nrow][1]!=nextlist[c][1]){ //if O's don't match
				   fprintf(output2,"%d %d %d %d %d\n",
			           snap, 
			           snap+1,
			           nextlist[nrow][1], //O1 in Zundel
			           nextlist[c][1], //O2 in Zundel
			           nextlist[nrow][ncoln]); //shared proton

				   nextlist[nrow][5]=1;
				   nextlist[c][5]=1; //flag Zundel with 1
				   nextlist[nrow][6]=nextlist[c][1];
				   nextlist[c][6]=nextlist[nrow][1]; //store partner info
				}
			     }
		          }
		       }//end d-loop
		    }//end c-loop
		 }//end ncoln-loop
	      }//end nrow-loop

	      for(nrow=1; nrow<nlines; nrow++){
		 if(nextlist[nrow][5]!=1){ //if OCO has NOT been flagged, it's an Eigen (2)
		    nextlist[nrow][5]=2;
		 }
	      }//end nrow-loop

//for(z=1; z<nlines; z++){
//printf("%d %d %d %d %d %d\n",nextlist[z][1],nextlist[z][2],nextlist[z][3],nextlist[z][4],nextlist[z][5],nextlist[z][6]);
//}

//Check for changes and update pathlist accordingly

	      for(srow=1; srow<slines; srow++){
	         for(nrow=1; nrow<nlines; nrow++){
		    if(startlist[srow][1]==nextlist[nrow][1]){ //if OCO match
		       if(startlist[srow][5]==nextlist[nrow][5]){ //if flag did NOT change, no transfer
		       }
		       else if(startlist[srow][5]!=nextlist[nrow][5]){ //if flag DID change, transfer!
//No need to double count, when a Zundel becomes an Eigen - it will be counted as a flag change
//printf("%d %d %d %d\n",startlist[srow][1],nextlist[nrow][1],startlist[srow][5],nextlist[nrow][5]);
		       }
		    }
		    else if(startlist[srow][1]!=nextlist[nrow][1]){//if OCO don't match
		       for(c=2; c<5; c++){ //loop through all H's
		          for(d=2; d<5; d++){
			     if(startlist[srow][c]==nextlist[nrow][d]){ //if a H's match, transfer!
			        if(startlist[srow][5]!=1 && startlist[srow][5]!=0){ //if not a Zundel
				   if(nextlist[nrow][5]==1){ //Eigen to Zundel
				      countEZ++;
				   }
				   else if(nextlist[nrow][5]==2){ //Eigen to Eigen
				      countEE++;
				   }
//printf("Type III: EZ or EE\n");
//printf("%d %d %d %d %d\n",snap,startlist[srow][1],nextlist[nrow][1],startlist[srow][5],nextlist[nrow][5]);
			        }
			        else if(startlist[srow][5]==1 && nextlist[nrow][5]!=1){
//printf("list: %d %d %d %d\n",snap,startlist[srow][1],startlist[srow][6],startlist[srow][5]);
				   for(z=1; z<nlines; z++){ //does the OCO in Zundel appear in the next snap?
				      if(startlist[srow][1]==nextlist[z][1]){ //if OCO matches partner doest match ANY OCO in next snap
					 countZdiss++; //Occurs only when Zundel dissociates to water and Eigen
//printf("trigger 1:\n");
				      }
			              else if(startlist[srow][6]==nextlist[z][1]){
					 countZdiss++; //Occurs only when Zundel dissociates to water and Eigen
//printf("trigger 2: %d %d\n",startlist[srow][6],nextlist[z][1]);				   	 
				      }
				      else if(startlist[srow][1]!=nextlist[z][1] && startlist[srow][6]!=nextlist[z][1]){;
					for(x=2; x<5; x++){
					   for(y=2; y<5; y++){
					      if(startlist[srow][x]==nextlist[z][y]){	
//printf("%d %d %d\n",startlist[srow][1],nextlist[z][1],startlist[srow][x]);
			                         countZtoE++; //Zundel transfer to water
				              }
					   }//end y-loop
				        }//end x-loop
				      }
				   }//end z-loop
				}
				else if(startlist[srow][5]==1 && nextlist[nrow][5]==1){ //If both structures initially and after are Zundel
				   if(startlist[srow][1]==nextlist[nrow][1]){ //if OCO matches
				      if(startlist[srow][6]!=nextlist[nrow][6]){//if second Zundel OCO pair does not match
				         for(x=2; x<5; x++){
					    for(y=2; y<5; y++){
					       if(startlist[srow][x]==nextlist[nrow][y]){ //if a set of H's is identified
						   countZZ++;
					       }
					    }//end y-loop
					 }//end x-loop
				      }
				   }
				   else if(startlist[srow][1]!=nextlist[nrow][1]){ //if OCO's don't match
				      if(startlist[srow][6]==nextlist[nrow][6]){ //second Z-OCO partner does match
					 for(x=2; x<5; x++){
					    for(y=2; y<5; y++){
					       if(startlist[srow][x]==nextlist[nrow][y]){
						  countZZ++;
					       }
					    }//end y-loop
					 }//end x-loop
				      }
				   }
				}
			     }
			  }//end d-loop
		       }//end c-loop
		    }//else if
	         }//end nrow-loop
	      }//end srow-loop
	
	      for(c=1; c<nlines; c++){
		 pathlist[snap+1][0][0]=nextlist[c][0];
		 pathlist[snap+1][c][0]=nextlist[c][1];
		 pathlist[snap+1][c][1]=nextlist[c][2];
		 pathlist[snap+1][c][2]=nextlist[c][3];
		 pathlist[snap+1][c][3]=nextlist[c][4];
		 pathlist[snap+1][c][4]=nextlist[c][5];	 
	      }//end c-loop
	       
	      for(z=0; z<nlines; z++){		
		 if(pathlist[snap+1][z][0]!=0){ //if element is nonzero
		    fprintf(output3,"%d ", pathlist[snap+1][z][0]);
		 }
	      }//end d-loop
	      if(pathlist[snap+1][1][0]!=0){
	         fprintf(output3,"\n");
	      }


//Intialize startlist
	     slines=nlines;
             for(c=0; c<50; c++){//For startlist
                for(d=0; d<7; d++){
                    startlist[c][d]=0;
                }
             }
//Copy nextlist into startlist
	    
	    for(d=1; d<nlines; d++){
	       startlist[d][0]=nextlist[d][0];
	       startlist[d][1]=nextlist[d][1];
	       startlist[d][2]=nextlist[d][2];
	       startlist[d][3]=nextlist[d][3];
	       startlist[d][4]=nextlist[d][4];
	       startlist[d][5]=nextlist[d][5];	
	       startlist[d][6]=nextlist[d][6];
	    }

//Intialize nextlist	
             for(c=0; c<50; c++){//for nextlist
                for(d=0; d<7; d++){
                   nextlist[c][d]=0;
                }
            }
	    nlines=1;
	   }//end else snap!=1
	}//end snap-loop
	printf("# of Z transfers: %d\n", countZtoE);
	printf("# of E to Z: %d\n",countEZ);
	printf("# of E to E: %d\n",countEE);
	printf("# Z dissociations: %d\n",countZdiss);
	printf("# of Z to Z transfers: %d\n",countZZ);

}//end main
