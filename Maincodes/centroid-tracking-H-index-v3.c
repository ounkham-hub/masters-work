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
//Legend for flags:
//  1, H's apart of Eigens
//  2, shared proton in Zundel
//  4, H's bridging waters in Zundels
//  55, H's a part of water molecules
//
//*****************************************************************************

int main(){

   int snap,Oxy,Hyd,PBC1,PBC2,PBC3,PBC4,PBC5; 
   int countHindex,y,z,Oindex,Hhold1,clines;
   int Hhold2,row,coln;
   int crow,mrow,indexno,maxsnap,countOindex,lines;
   float dist,angle; 
   FILE *input,*input1,*output,*hinput,*hinput1;
   char ifile[100],ifile1[100],hfile[100],hfile1[100];
  
   int currentsnap[500][7];
   int mastermatrix[500][7];
   int labelmatrix[500][7];
 
//Intialize matrices and counters
   for(y=0; y<500; y++){
      for(z=0; z<7; z++){
         currentsnap[y][z]=0;
	 mastermatrix[y][z]=0;
	 labelmatrix[y][z]=0;
      }
   }//end z-loop
   countHindex=1;
   countOindex=1;
   clines=0;
   lines=0;
   indexno=103;
   maxsnap=61555;
  
   output=fopen("centroid-total-event-persistences.txt","a");

   for(snap=1; snap<=maxsnap; snap++){
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
//Section II: Distinguish Eigens containing 3 H's from water molecules by
//examining another set of GraphGeod files.
//
//*****************************************************************************
         sprintf(hfile,"R-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
         hinput=fopen(hfile,"r");

      while(fscanf(hinput,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
         if(countOindex==1){
             Oindex=Oxy;
             Hhold1=Hyd;
          }
          else if(countOindex==2){
             if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;

                    labelmatrix[lines][0]=snap;
                    labelmatrix[lines][1]=Oxy;
                    labelmatrix[lines][2]=Hhold1;
                    labelmatrix[lines][3]=Hhold2;
                    labelmatrix[lines][5]=55; //flagged as water
                    lines++;
                 }
                 else{ //Should never trigger unless OH- exist       
                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
                 }
              }
              else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
                 if(Oindex==Oxy){
                    labelmatrix[lines][0]=snap;
                    labelmatrix[lines][1]=Oxy;
                    labelmatrix[lines][2]=Hhold1;
                    labelmatrix[lines][3]=Hhold2;
                    labelmatrix[lines][4]=Hyd;
                    labelmatrix[lines][5]=1; //flagged as eigen or H3O+

//Remove previous data on O atom, when H3O+ was no identified
                    labelmatrix[lines-1][0]=0;
                    labelmatrix[lines-1][1]=0;
                    labelmatrix[lines-1][2]=0;
                    labelmatrix[lines-1][3]=0;
                    labelmatrix[lines-1][4]=0;
                    lines++;

                    countOindex=0; //restarts Oindex to read the next line;
                 }
                 else{ //Water molecule identified and countOindex==2
                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
                 }
              }
              countOindex++;
         }//end while-loop
         fclose(hinput);

/*
        for(z=0; z<lines; z++){
           if(labelmatrix[z][0]!=0){
              printf("Sec 1: %d %d %d %d %d %d\n",
              labelmatrix[z][0],
              labelmatrix[z][1],
              labelmatrix[z][2],
              labelmatrix[z][3],
              labelmatrix[z][4],
              labelmatrix[z][5]);
           }
        }//end z-loop
*/

//*****************************************************************************
//Section IIa: Read in first snapshot and identify the coordination number for
//per H atom. 1 O partner, water or Eigen (indistinguishable) or 2 O for 
//Zundels.
//*****************************************************************************
          while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
             if(countHindex==1){
                Oindex=Oxy;
                Hhold1=Hyd;

                currentsnap[clines][0]=Hyd; //H index
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
	     if(currentsnap[y][1]==2){
	        printf("Sec. 2 %d %d %d %d %d\n",
	        currentsnap[y][0],
	        currentsnap[y][1],
	        currentsnap[y][2],
	        currentsnap[y][3],
	        currentsnap[y][4]);
	     }
	  }
*/
//*****************************************************************************
//Section IIb: Additional filter to determine  which H's and O's corresponds 
//to an Eigen, shared proton in Zundel, or bridging waters in Zundels or water.
//*****************************************************************************

	for(crow=0; crow<clines; crow++){
//If Eigen is identified in labelmatrix, flag==1
           if(currentsnap[crow][1]!=2){
	      for(row=0; row<lines; row++){
		 for(coln=2; coln<5; coln++){
		    if(currentsnap[crow][0]==labelmatrix[row][coln]){ 	
//If H indices matches, then change flag
		       if(labelmatrix[row][5]==1){
//If flags don't match
/*			   printf("%d %d %d %d %d %d\n",
			   labelmatrix[row][0],	
			   labelmatrix[row][1],
			   labelmatrix[row][2],
			   labelmatrix[row][3],
		           labelmatrix[row][4],
			   labelmatrix[row][5]);
*/			
			   currentsnap[crow][1]=1; //change CN to Eigen
		       }
	               else if(labelmatrix[row][5]==55){
//If water is identified in labelmatrix, flag==55
			   currentsnap[crow][1]=55; //change CN to water
		       }
		    }  
	         }//end coln-loop
	      }//end row-loop
           }
	}//end crow-loop

//Subsection to distinguish bridging H's in Zundels vs. Eigens

	for(crow=0; crow<clines; crow++){
	   if(currentsnap[crow][1]==2){ //If flagged Zundel
	      for(row=0; row<lines; row++){ 
		 if(labelmatrix[row][1]==currentsnap[crow][2]){
//If O1 indices matches, change flag
	            for(coln=2; coln<5; coln++){
		       for(z=0; z<clines; z++){
		          if(currentsnap[z][0]==labelmatrix[row][coln] && currentsnap[z][0]!=currentsnap[crow][0]){
//If H indices matches, but isn't the shared proton
			     currentsnap[z][1]=4;
			  }
		       }//end z-loop
	            }//end coln-loop
                 }
		 else if(labelmatrix[row][1]==currentsnap[crow][3]){
//If O2 indices matches, change flag
		    for(coln=2; coln<5; coln++){
		       for(z=0; z<clines; z++){
		          if(currentsnap[z][0]==labelmatrix[row][coln] && currentsnap[z][0]!=currentsnap[crow][0]){
//If H indices matches, but isn't the shared proton
			     currentsnap[z][1]=4;
			  }
		       }//end z-loop	
		    }		       
	         }//end coln-loop
	      }//end row-loop
           }
	}//end crow-loop

/* 
	for(y=0; y<clines; y++){
	   if(currentsnap[y][0]!=0){
	      printf("Relabel: %d %d %d %d %d\n",
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
		      if(currentsnap[crow][1]>=4){
//If Eigen or water is identified
		         mastermatrix[mrow][1]=currentsnap[crow][1];//CN
		         mastermatrix[mrow][2]=currentsnap[crow][2]; //O1
		         mastermatrix[mrow][4]=snap; //initial snap
		         mastermatrix[mrow][5]=snap; //final snap

			 currentsnap[crow][0]=0;
		      }
		      else if(currentsnap[crow][1]<4){
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
        for(y=0; y<500; y++){
           for(z=0; z<7; z++){
              currentsnap[y][z]=0;
	      labelmatrix[y][z]=0;
           }
        }//end z-loop
	clines=0;
	countHindex=1;
	lines=0;
	countOindex=1;

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
//Section IV: Distinguish Eigens containing 3 H's from water molecules by
//examining another set of GraphGeod files.
//
//*****************************************************************************

       else if(snap>1){
         sprintf(hfile1,"R-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
         hinput1=fopen(hfile1,"r");

      while(fscanf(hinput1,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
         if(countOindex==1){
             Oindex=Oxy;
             Hhold1=Hyd;
          }
          else if(countOindex==2){
             if(Oindex==Oxy){ //if O is the same on next line
                    Hhold2=Hyd;

                    labelmatrix[lines][0]=snap;
                    labelmatrix[lines][1]=Oxy;
                    labelmatrix[lines][2]=Hhold1;
                    labelmatrix[lines][3]=Hhold2;
                    labelmatrix[lines][5]=55; //flagged as water
                    lines++;
                 }
                 else{ //Should never trigger unless OH- exist       
                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
                 }
              }
              else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
                 if(Oindex==Oxy){
                    labelmatrix[lines][0]=snap;
                    labelmatrix[lines][1]=Oxy;
                    labelmatrix[lines][2]=Hhold1;
                    labelmatrix[lines][3]=Hhold2;
                    labelmatrix[lines][4]=Hyd;
                    labelmatrix[lines][5]=1; //flagged as eigen or H3O+

//Remove previous data on O atom, when H3O+ was no identified
                    labelmatrix[lines-1][0]=0;
                    labelmatrix[lines-1][1]=0;
                    labelmatrix[lines-1][2]=0;
                    labelmatrix[lines-1][3]=0;
                    labelmatrix[lines-1][4]=0;
                    lines++;

                    countOindex=0; //restarts Oindex to read the next line;
                 }
                 else{ //Water molecule identified and countOindex==2
                    countOindex=1;
                    Oindex=Oxy;
                    Hhold1=Hyd;
                 }
              }
              countOindex++;
         }//end while-loop
         fclose(hinput1);
/*
        for(z=0; z<lines; z++){
           if(labelmatrix[z][0]!=0){
              printf("%d %d %d %d %d %d\n",
              labelmatrix[z][0],
              labelmatrix[z][1],
              labelmatrix[z][2],
              labelmatrix[z][3],
              labelmatrix[z][4],
              labelmatrix[z][5]);
           }
        }//end z-loop
*/

//*****************************************************************************
//Section IVa: Read in proceeding GraphGeod file to identify the coordination # 
//for each H atom.
//*****************************************************************************

          sprintf(ifile1,"H-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap); //open current snap
          input1=fopen(ifile1,"r");

          while(fscanf(input1,"%d %d %d %d %d %d %d %f %f\n",&Oxy,&Hyd,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&angle)==9){
             if(countHindex==1){
                Oindex=Oxy;
                Hhold1=Hyd;

                currentsnap[clines][0]=Hyd; //H index
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
//Section IVb: Additional filter to determine  which H's and O's corresponds 
//to an Eigen, shared proton in Zundel, or bridging waters in Zundels or water.
//*****************************************************************************

	for(crow=0; crow<clines; crow++){
//If Eigen is identified in labelmatrix, flag==1
           if(currentsnap[crow][1]!=2){
	      for(row=0; row<lines; row++){
		 for(coln=2; coln<5; coln++){
		    if(currentsnap[crow][0]==labelmatrix[row][coln]){ 	
//If H indices matches, then change flag
		       if(labelmatrix[row][5]==1){
//If flags don't match
/*			   printf("%d %d %d %d %d %d\n",
			   labelmatrix[row][0],	
			   labelmatrix[row][1],
			   labelmatrix[row][2],
			   labelmatrix[row][3],
		           labelmatrix[row][4],
			   labelmatrix[row][5]);
*/			
			   currentsnap[crow][1]=1; //change CN to Eigen
		       }
	               else if(labelmatrix[row][5]==55){
//If water is identified in labelmatrix, flag==55
			   currentsnap[crow][1]=55; //change CN to water
		       }
		    }  
	         }//end coln-loop
	      }//end row-loop
           }
	}//end crow-loop

//Subsection to distinguish bridging H's in Zundels vs. Eigens

	for(crow=0; crow<clines; crow++){
	   if(currentsnap[crow][1]==2){ //If flagged Zundel
	      for(row=0; row<lines; row++){ 
		 if(labelmatrix[row][1]==currentsnap[crow][2]){
//If O1 indices matches, change flag
	            for(coln=2; coln<5; coln++){
		       for(z=0; z<clines; z++){
		          if(currentsnap[z][0]==labelmatrix[row][coln] && currentsnap[z][0]!=currentsnap[crow][0]){
//If H indices matches, but isn't the shared proton
			     currentsnap[z][1]=4;
			  }
		       }//end z-loop
	            }//end coln-loop
                 }
		 else if(labelmatrix[row][1]==currentsnap[crow][3]){
//If O2 indices matches, change flag
		    for(coln=2; coln<5; coln++){
		       for(z=0; z<clines; z++){
		          if(currentsnap[z][0]==labelmatrix[row][coln] && currentsnap[z][0]!=currentsnap[crow][0]){
//If H indices matches, but isn't the shared proton
			     currentsnap[z][1]=4;
			  }
		       }//end z-loop	
		    }		       
	         }//end coln-loop
	      }//end row-loop
           }
	}//end crow-loop

/*
	  for(y=0; y<clines; y++){
	     if(currentsnap[y][0]!=0){
	        printf("Relabel 2: %d %d %d %d %d\n",
	        currentsnap[y][0],
	        currentsnap[y][1],
	        currentsnap[y][2],
	        currentsnap[y][3],
	        currentsnap[y][4]);
	     }
	  }
*/
//*****************************************************************************
//Section V: Section compare connectivities in initial snapshots to proceeding
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
		        }//end Case 1b Situation
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
		       }//end Case 1c. situation
		    }//end of Case 1 situations
		    else if(mastermatrix[mrow][1]==1 && currentsnap[crow][1]==2){
//Case 2: If structure changes from Eigen to Zundel (shared proton), update mastermatrix
//		       printf("2: Eigen transformed to Zundel (shared proton)\n");

 		       fprintf(output,"1 2 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
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
		    else if(mastermatrix[mrow][1]==1 && currentsnap[crow][1]==4){
//Case 3: If structure changes from Eigen to Zundel (bridging water), update mastermatrix
//		       printf("3.: Eigen transformed to Zundel (bridging water)\n");

 		       fprintf(output,"1 4 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
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
		    else if(mastermatrix[mrow][1]==1 && currentsnap[crow][1]==55){
//Case 4: If structure changes from Eigen to Water, update mastermatrix
//		       printf("3.: Eigen transformed to Zundel (bridging water)\n");

 		       fprintf(output,"1 55 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
		      
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;

		    }//end Case 4 situation
		    else if(mastermatrix[mrow][1]==2 && currentsnap[crow][1]==1){
//Case 5: If structure changes from Zundel (shared proton) to Eigen, update mastermatrix
//		       printf("5: Zundel (shared proton) transformed to Eigen\n");

 		       fprintf(output,"2 1 %d %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
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
		    }//end Case 5 situation
		    else if(mastermatrix[mrow][1]==2 && currentsnap[crow][1]==4){
//Case 6: If structure changes from Zundel (shared proton) to Zundel (bridging water), update mastermatrix
//		       printf("6: Zundel (shared proton) transformed to Zundel bridging water\n");

 		       fprintf(output,"2 4 %d %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
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
		    }//end Case 6 situation
		    else if(mastermatrix[mrow][1]==2 && currentsnap[crow][1]==55){
//Case 7: If structure changes from Zundel (shared proton) to Water, update mastermatrix
//		       printf("7: Zundel (shared proton) transformed to Water\n");

 		       fprintf(output,"2 55 %d %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
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
		    }//end Case 7 situation
		    else if(mastermatrix[mrow][1]==4 && currentsnap[crow][1]==1){
//Case 8: If structure changes from Zundel (bridging water) to Eigen, update mastermatrix
//		       printf("8: Zundel (bridging water) transformed to Eigen\n");

 		       fprintf(output,"4 1 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
  
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;
		    }//end Case 8 situation
		    else if(mastermatrix[mrow][1]==4 && currentsnap[crow][1]==2){
//Case 9: If structure changes from Zundel (bridging water) to Zundel (shared proton), update mastermatrix
//		       printf("9: Zundel (bridging water) transformed to Zundel (shared proton)\n");

 		       fprintf(output,"4 2 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
  
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;
		    }//end Case 9 situation
		    else if(mastermatrix[mrow][1]==4 && currentsnap[crow][1]==55){
//Case 10: If structure changes from Zundel (bridging water) to Water, update mastermatrix
//		       printf("9: Zundel (bridging water) transformed to Water\n");

 		       fprintf(output,"4 55 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
  
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;
		    }//end Case 10 situation
		    else if(mastermatrix[mrow][1]==55 && currentsnap[crow][1]==1){
//Case 11: If structure changes from water to Eigen, update mastermatrix
//		       printf("3: Water transformed to Eigen\n");

 		       fprintf(output,"55 1 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
  
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;
		    }//end Case 11 situation
		    else if(mastermatrix[mrow][1]==55 && currentsnap[crow][1]==2){
//Case 12: If structure changes from water to Zundel (shared proton), update mastermatrix
//		       printf("3: Water transformed to Zundel (shared proton)\n");

 		       fprintf(output,"55 2 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
  
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;
		    }//end Case 12 situation
		    else if(mastermatrix[mrow][1]==55 && currentsnap[crow][1]==4){
//Case 13: If structure changes from water to Zundel (bridging water), update mastermatrix
//		       printf("3: Water transformed to Zundel (bridging water)\n");

 		       fprintf(output,"55 4 %d %d %d %d\n",
		       mastermatrix[mrow][0], //Hindex
		       mastermatrix[mrow][2], //O1	
		       mastermatrix[mrow][4], //initial snap
		       mastermatrix[mrow][5]); //final snap
  
		       mastermatrix[mrow][1]=currentsnap[crow][1]; //Change structure label
		       mastermatrix[mrow][2]=currentsnap[crow][2]; //change O1
		       mastermatrix[mrow][3]=currentsnap[crow][3]; //O2
		       mastermatrix[mrow][4]=snap; //update initial snap
		       mastermatrix[mrow][5]=snap; //update final snap
		       mastermatrix[mrow][6]=1; //turn flag on
		       currentsnap[crow][0]=0;
		    }//end Case 12 situation

		    else if(mastermatrix[mrow][1]==0 && currentsnap[crow][1]!=0){
//Case : Occurs when connectivity reappears in proceeding snapshot, update matrix	

		       mastermatrix[mrow][1]=currentsnap[crow][1];
		       mastermatrix[mrow][2]=currentsnap[crow][2];
		       mastermatrix[mrow][3]=currentsnap[crow][3];
		       mastermatrix[mrow][4]=snap;
		       mastermatrix[mrow][5]=snap;

		       printf("5: Connectivity reappeared!\n");
 		       printf("%d %d\n",
		       snap,
		       mastermatrix[mrow][0]); //Hindex

		       mastermatrix[mrow][6]=1;
		       currentsnap[crow][0]=0;
		    }//end Case 4 situation
		 }//if non-zero element
	      }   
	   }//end mrow-loop
	}//end crow-loop

//*****************************************************************************
//Section VI: Connectivities that "disappear" in a snapshot should also
//be accounted for. These would be considered structural changes.
//*****************************************************************************

	    for(mrow=1; mrow<213; mrow++){
	       if(mastermatrix[mrow][6]==0 && mastermatrix[mrow][0]!=0){
//If H index matches, update matrix

	          printf("6: Connectivity disappeared!\n");
 		  printf("%d %d\n",
		  snap,
		  mastermatrix[mrow][0]); //Hindex

//Reset elements for H index 
		  mastermatrix[mrow][1]=0;
		  mastermatrix[mrow][2]=0;	
		  mastermatrix[mrow][3]=0;
		  mastermatrix[mrow][4]=0;
		  mastermatrix[mrow][5]=0;
	       }//end crow-loop
	    }//end mrow-loop

//Initialize currentsnap matrix for the snapshot
        for(y=0; y<500; y++){
           for(z=0; z<7; z++){
              currentsnap[y][z]=0;
	      labelmatrix[y][z]=0;
           }
        }//end z-loop
	clines=0;
	lines=0;
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
   fclose(output);
}//end main-loop 
