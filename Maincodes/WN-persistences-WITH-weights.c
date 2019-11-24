#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//***************************************************************************
//Program created to determine the persistence of CB weights based on a set of 
//predetermined weights.
//
//Created by Lelee Ounkham
//Last modified on February 17, 2019
//
//***************************************************************************

int main(){
	
    int Oxy,Hyd,snap,crow,Oindex,clines,y,z,countHindex,maxsnap;
    int indexno,mrow,lines,O,H,countOindex,Hhold1,Hhold2,row,coln;
    float dist,weight,whold1,w1,w2,prec,diff;
    FILE *input,*hinput,*input2,*hinput2,*output;
    char ifile[100],ifile2[100],hfile[100],hfile2[100];

    int currentmatrix[300][6];
    int labelmatrix[300][6];
    int mastermatrix[300][9];
    float currentweight[300][6];
    float masterweight[300][6];
	
 //Initialize matrix and counters
    for(y=0; y<300; y++){
       for(z=0; z<7; z++){
 	 	 currentmatrix[y][z]=0;
		 labelmatrix[y][z]=0;
 	 	 currentweight[y][z]=0.000;
       }
	   
       for(z=0; z<9; z++){
 		   mastermatrix[y][z]=0;
 		   masterweight[y][z]=0.000;
       }//end z-loop
    }//end y-loop

    clines=0;
	lines=0;
    countHindex=1;
	countOindex=1;
    indexno=103;
    maxsnap=61555;
	prec=0.000001;


	
//*****************************************************************************
//Section I: Create mastermatrix which will encompass connectivity information
//per H atom..
//*****************************************************************************
         for(mrow=1; mrow<213; mrow++){
            mastermatrix[mrow][0]=indexno;
            indexno++;
         }//end mrow-loop

//        for(y=1; y<213; y++){
//            printf("%d\n",mastermatrix[y][0]);
//         }//end y-loop
	
	output=fopen("updated-CB-weight-persistence.txt","a");
    for(snap=1; snap<=61555; snap++){
    	if(snap==1){		
//*****************************************************************************
//Section IIa: Distinguish Eigens containing 3 H's from water molecules by
//examining another set of GraphGeod files.
//
//*****************************************************************************
			sprintf(hfile,"HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap,snap); //open current snap
			hinput=fopen(hfile,"r");
			
	        while(fscanf(hinput,"%d %d %f %f\n",&O,&H,&w1,&w2)==4){
	           if(countOindex==1){
	               Oindex=O;
	               Hhold1=H;
	            }
	            else if(countOindex==2){
	               if(Oindex==O){ //if O is the same on next line
	                      Hhold2=H;

	                      labelmatrix[lines][0]=snap;
	                      labelmatrix[lines][1]=O;
	                      labelmatrix[lines][2]=Hhold1;
	                      labelmatrix[lines][3]=Hhold2;
	                      labelmatrix[lines][5]=55; //flagged as water
	                      lines++;
	                   }
	                   else{ //Should never trigger unless OH- exist       
	                      countOindex=1;
	                      Oindex=O;
	                      Hhold1=H;
	                   }
	                }
	                else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
	                   if(Oindex==O){
	                      labelmatrix[lines][0]=snap;
	                      labelmatrix[lines][1]=O;
	                      labelmatrix[lines][2]=Hhold1;
	                      labelmatrix[lines][3]=Hhold2;
	                      labelmatrix[lines][4]=H;
	                      labelmatrix[lines][5]=1; //flagged as eigen or H3O+

	  //Remove previous data on O atom, when H3O+ was no identified
	                      labelmatrix[lines-1][0]=0;
	                      labelmatrix[lines-1][1]=0;
	                      labelmatrix[lines-1][2]=0;
	                      labelmatrix[lines-1][3]=0;
	                      labelmatrix[lines-1][4]=0;
			      labelmatrix[lines-1][5]=0;
	                      lines++;

	                      countOindex=0; //restarts Oindex to read the next line;
	                   }
	                   else{ //Water molecule identified and countOindex==2
	                      countOindex=1;
	                      Oindex=O;
	                      Hhold1=H;
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
//Section IIb: Read in first snapshot and identify the coordination number for
//per H atom. 1 O partner, water or Eigen (indistinguishable) or 2 O for 
//Zundels.
//*****************************************************************************
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

 			currentweight[clines][0]=weight;	
                 	clines++;

              	}
              	else if(countHindex==2){
                	if(Hyd==Hhold1){ //when H index matches
                    	currentmatrix[clines-1][0]=0;
                    	currentmatrix[clines-1][1]=0;
                    	currentmatrix[clines-1][2]=0;
                    	currentmatrix[clines-1][3]=0;
 		   	currentweight[clines-1][0]=0.000;

                    	currentmatrix[clines][0]=Hyd; //indicative of Zundel
                    	currentmatrix[clines][1]=Oindex; //O1
                   	currentmatrix[clines][2]=Oxy; //O2
                    	currentmatrix[clines][3]=2; //CN

 		   	currentweight[clines][0]=whold1;
 		   	currentweight[clines][1]=weight;
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
 		   	currentweight[clines][0]=weight;
                   	clines++;
 		        countHindex=1;
                 	}
              	}
             	countHindex++;
           	}//end while-loop
           	fclose(input);
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
//*****************************************************************************
//Section IIb: Additional filter to determine  which H's and O's corresponds 
//to an Eigen, shared proton in Zundel, or bridging waters in Zundels or water.
//*****************************************************************************

	  	 for(crow=0; crow<clines; crow++){
//If Eigen is identified in labelmatrix, flag==1
	        if(currentmatrix[crow][3]!=2){
	  	        for(row=0; row<lines; row++){
	  		        for(coln=2; coln<5; coln++){
	  		           if(currentmatrix[crow][0]==labelmatrix[row][coln]){ 	
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
	  			             currentmatrix[crow][3]=1; //change CN to Eigen
	  		             }
	  	                 else if(labelmatrix[row][5]==55){
//If water is identified in labelmatrix, flag==55
	  			            currentmatrix[crow][3]=55; //change CN to water
	  		             }
	  		          }  
	  	           }//end coln-loop
	  	        }//end row-loop
	         }
	  	  }//end crow-loop

//Subsection to distinguish bridging H's in Zundels vs. Eigens

		  for(crow=0; crow<clines; crow++){
	  	     if(currentmatrix[crow][3]==2){ //If flagged Zundel
	  	        for(row=0; row<lines; row++){ 
	  		   if(labelmatrix[row][1]==currentmatrix[crow][1]){
//If O1 indices matches, change flag
	  	              for(coln=2; coln<5; coln++){
	  		           for(z=0; z<clines; z++){
	  		                if(currentmatrix[z][0]==labelmatrix[row][coln] && currentmatrix[z][0]!=currentmatrix[crow][0]){
//If H indices matches, but isn't the shared proton and not a double zundel
	  			               currentmatrix[z][3]=4;

	  			            }
	  		             }//end z-loop
	  	              }//end coln-loop
	               }
	  	        else if(labelmatrix[row][1]==currentmatrix[crow][2]){
//If O2 indices matches, change flag
	  		         for(coln=2; coln<5; coln++){
	  		            for(z=0; z<clines; z++){
	  		               if(currentmatrix[z][0]==labelmatrix[row][coln] && currentmatrix[z][0]!=currentmatrix[crow][0]){
//If H indices matches, but isn't the shared proton or a double zundel
	  			              currentmatrix[z][3]=4;

	  			           }
	  		            }//end z-loop	
	  		         }		       
	  	           }//end coln-loop
	  	       }//end row-loop
	            }
	  	}//end crow-loop
		
/*
	  	for(y=0; y<clines; y++){
	  	   if(currentmatrix[y][0]!=0){
	  	      printf("Relabel 1: %d %d %d %d %0.3f %0.3f\n",
	  	      currentmatrix[y][0],
	  	      currentmatrix[y][1],
	  	      currentmatrix[y][2],
	  	      currentmatrix[y][3],
			  currentweight[y][0],
			  currentweight[y][1]);
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
	   		if(currentmatrix[crow][0]!=0){
	   	           if(mastermatrix[mrow][0]==currentmatrix[crow][0]){
//If H indices matches
	   		          if(currentmatrix[crow][3]>=4){
//Water or outer flag is identified
	   		             mastermatrix[mrow][3]=currentmatrix[crow][3];//CN
						 
	   		             mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
				     mastermatrix[mrow][2]=0;
	   		             mastermatrix[mrow][4]=snap; //initial snap
	   		             mastermatrix[mrow][5]=snap; //final snap
						 
				     masterweight[mrow][0]=currentweight[crow][0];
				     masterweight[mrow][1]=currentweight[crow][1];

	   			         currentmatrix[crow][0]=0;
	   		          }
	   		          else if(currentmatrix[crow][3]<4){
//If Zundel or Eigen is identified
	   		             mastermatrix[mrow][3]=currentmatrix[crow][3];//CN
						 
	   		             mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
	   			     mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2, should be 0 for Eigen
	   		             mastermatrix[mrow][4]=snap; //initial snap
	   		             mastermatrix[mrow][5]=snap; //final snap
				     mastermatrix[mrow][6]=snap; //initial snap
				     mastermatrix[mrow][7]=snap;
						 
				     masterweight[mrow][0]=currentweight[crow][0];
				     masterweight[mrow][1]=currentweight[crow][1];
			 
	   			     currentmatrix[crow][0]=0;
	   		          }	      
	   		       }//end if H-index matches
	   		   }
	   	     }//end crow-loop
	   	  }//end mrow-loop

 
/*
	      for(y=1; y<213; y++){
	   	     if(mastermatrix[y][0]!=0){
	   	        printf("%d %d %d %d %d %d %0.3f %0.3f\n",
	   	        mastermatrix[y][0],
	   	        mastermatrix[y][1],
	   	        mastermatrix[y][2],
	   	        mastermatrix[y][3],
	   	        mastermatrix[y][4],
	   		mastermatrix[y][5],
			masterweight[y][0],
			masterweight[y][1]);			
	   	     }
	     }//end y-loop
*/		 		 
//Delete currentmatrix matrix for next snapshot
	      for(y=0; y<clines; y++){
	 	    for(z=0; z<7; z++){
	 	       currentmatrix[y][z]=0;
	 	       currentweight[y][z]=0.000;
	 	     }
	      }//end y-loop
	      for(y=0; y<lines; y++){
	 	    for(z=0; z<7; z++){
	 	       labelmatrix[y][z]=0;
	 	     }
	      }//end y-loop

	      clines=0;
		  lines=0;
	      countHindex=1;
		  countOindex=1; 

      }//end snap==1 loop
//***************************************************************************
//Section IV: Transfer connectivity and weights information from proceeding 
//snapshot to currentmatrix.
//***************************************************************************
	  else if(snap>1){
		sprintf(hfile2,"HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap,snap); //open current snap
		hinput2=fopen(hfile2,"r");
			
	        while(fscanf(hinput2,"%d %d %f %f\n",&O,&H,&w1,&w2)==4){
	           if(countOindex==1){
	               Oindex=O;
	               Hhold1=H;
	            }
	            else if(countOindex==2){
	               if(Oindex==O){ //if O is the same on next line
	                      Hhold2=H;

	                      labelmatrix[lines][0]=snap;
	                      labelmatrix[lines][1]=O;
	                      labelmatrix[lines][2]=Hhold1;
	                      labelmatrix[lines][3]=Hhold2;
	                      labelmatrix[lines][5]=55; //flagged as water
	                      lines++;
	                   }
	                   else{ //Should never trigger unless OH- exist       
	                      countOindex=1;
	                      Oindex=O;
	                      Hhold1=H;
	                   }
	                }
	                else if(countOindex==3){ //if there is a 3rd occurrence then send to memory/array
	                   if(Oindex==O){
	                      labelmatrix[lines][0]=snap;
	                      labelmatrix[lines][1]=O;
	                      labelmatrix[lines][2]=Hhold1;
	                      labelmatrix[lines][3]=Hhold2;
	                      labelmatrix[lines][4]=H;
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
	                      Oindex=O;
	                      Hhold1=H;
	                   }
	                }
	                countOindex++;
	           }//end while-loop
	           fclose(hinput2);

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
//Section Va: Read in first snapshot and identify the coordination number for
//per H atom. 1 O partner, water or Eigen (indistinguishable) or 2 O for 
//Zundels.
//*****************************************************************************
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


 				currentweight[clines][0]=weight;	
                 		clines++;

              		}
              		else if(countHindex==2){
                		if(Hyd==Hhold1){ //when H index matches
                    		currentmatrix[clines-1][0]=0;
                    		currentmatrix[clines-1][1]=0;
                    		currentmatrix[clines-1][2]=0;
                    		currentmatrix[clines-1][3]=0;
 		   		currentweight[clines-1][0]=0.000;

                    		currentmatrix[clines][0]=Hyd; //indicative of Zundel
                    		currentmatrix[clines][1]=Oindex; //O1
                   	 	currentmatrix[clines][2]=Oxy; //O2
                    		currentmatrix[clines][3]=2; //CN

 		   		currentweight[clines][0]=whold1;
 		   		currentweight[clines][1]=weight;
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


 		   		   	currentweight[clines][0]=weight;
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
 	       	   	  currentweight[crow][0],
 	       	   	  currentweight[crow][1],
 	       	   	  currentmatrix[crow][3]);
 	    	  }
 		  }//end crow-loop
*/	 
//*****************************************************************************
//Section Vb: Additional filter to determine  which H's and O's corresponds 
//to an Eigen, shared proton in Zundel, or bridging waters in Zundels or water.
//*****************************************************************************

	  	 for(crow=0; crow<clines; crow++){
//If Eigen is identified in labelmatrix, flag==1
	            if(currentmatrix[crow][3]!=2){
	  	        for(row=0; row<lines; row++){
	  		    for(coln=2; coln<5; coln++){
	  		        if(currentmatrix[crow][0]==labelmatrix[row][coln]){ 	
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
	  			           currentmatrix[crow][3]=1; //change CN to Eigen
	  		             }
	  	                     else if(labelmatrix[row][5]==55){
//If water is identified in labelmatrix, flag==55
	  			           currentmatrix[crow][3]=55; //change CN to water
	  		             }
	  		          }  
	  	             }//end coln-loop
	  	        }//end row-loop
	             }
	  	  }//end crow-loop

//Subsection to distinguish bridging H's in Zundels vs. Eigens

		  for(crow=0; crow<clines; crow++){
	  	     if(currentmatrix[crow][3]==2){ //If flagged Zundel
	  	        for(row=0; row<lines; row++){ 
	  		    if(labelmatrix[row][1]==currentmatrix[crow][1]){
//If O1 indices matches, change flag
	  	                for(coln=2; coln<5; coln++){
	  		             for(z=0; z<clines; z++){
	  		                if(currentmatrix[z][0]==labelmatrix[row][coln] && currentmatrix[z][0]!=currentmatrix[crow][0]){
//If H indices matches, but isn't the shared proton
	  			               currentmatrix[z][3]=4;

	  			         }
	  		             }//end z-loop
	  	                 }//end coln-loop
	                      }
	  		      else if(labelmatrix[row][1]==currentmatrix[crow][2]){
//If O2 indices matches, change flag
	  		         for(coln=2; coln<5; coln++){
	  		            for(z=0; z<clines; z++){
	  		               if(currentmatrix[z][0]==labelmatrix[row][coln] && currentmatrix[z][0]!=currentmatrix[crow][0]){
//If H indices matches, but isn't the shared proton
	  			              currentmatrix[z][3]=4;
						  
	  			           }
	  		             }//end z-loop	
	  		         }		       
	  	           }//end coln-loop
	  	        }//end row-loop
	            }
	  	}//end crow-loop
	
	  	for(y=0; y<clines; y++){
	  	   if(currentmatrix[y][0]!=0){
if(currentmatrix[y][0]==155 && currentmatrix[y][3]==2){
	  	      printf("Relabel 2: %d %d %d %d %d %0.3f %0.3f\n",
		      snap,
	  	      currentmatrix[y][0],
	  	      currentmatrix[y][1],
	  	      currentmatrix[y][2],
	  	      currentmatrix[y][3],
		      currentweight[y][0],
		      currentweight[y][1]);
}
	  	   }
	  	}
				
//*****************************************************************************
//Section VI: Determine the weight persistence of CB for a given molecular
//species.
//*****************************************************************************
		
		for(crow=0; crow<clines; crow++){
			for(mrow=1; mrow<213; mrow++){
			    if(currentmatrix[crow][3]!=0){ // if currentmatrix is nonzero
				if(currentmatrix[crow][0]==mastermatrix[mrow][0]){ //If H index matches
						if(currentmatrix[crow][3]==1 && mastermatrix[mrow][3]==1){ //If Eigen flag is identified	
							if(currentmatrix[crow][1]==mastermatrix[mrow][1]){ //O's matches
								if(fabs(currentweight[crow][0] - masterweight[mrow][0]) <= prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) >= prec){ //O1H matches, update info
									mastermatrix[mrow][5]=snap; //update final snap
									mastermatrix[mrow][3]=currentmatrix[crow][3]; //update flag
									
									currentmatrix[crow][3]=0; //initialize currentmatrix
									mastermatrix[mrow][8]=1; //switch on
								}
								else if(fabs(currentweight[crow][0]-masterweight[mrow][0]) > prec && fabs(currentweight[crow][0]+masterweight[mrow][0]) > prec){ 
//if indices are the same, but weight changed - output persistence
								
									//Weight in eigen changed
									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3],
									currentmatrix[crow][3],
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][1], //O1 index
									mastermatrix[mrow][4], //initial snap
									mastermatrix[mrow][5], //final snap
									masterweight[mrow][0]); 	
									
									
									mastermatrix[mrow][3]=currentmatrix[crow][3]; //update flag
									
									masterweight[mrow][0]=currentweight[crow][0]; //update weight
									mastermatrix[mrow][4]=snap; //restart time
									mastermatrix[mrow][5]=snap;
								
									mastermatrix[mrow][8]=1;//switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
								}
							}			
						}
						else if(currentmatrix[crow][3]==55 && mastermatrix[mrow][3]==55){ //If Water flag is identified	
							if(currentmatrix[crow][1]==mastermatrix[mrow][1]){ //O's matches
								if(fabs(currentweight[crow][0] - masterweight[mrow][0]) <= prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) >= prec){ //O1H matches, update info
					
									mastermatrix[mrow][5]=snap; //update final snap
									
									mastermatrix[mrow][3]=currentmatrix[crow][3]; //update flag
									mastermatrix[mrow][8]=1; //switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
								}
								else if(fabs(currentweight[crow][0]-masterweight[mrow][0]) > prec && fabs(currentweight[crow][0]+masterweight[mrow][0]) > prec){ 
//if indices are the same, but weight changed - output persistence
				/*				printf("1w diff: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],
								currentmatrix[crow][3]);
*/

									
									//Bond in water changed
									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3],
									currentmatrix[crow][3],
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][1], //O1 index
									mastermatrix[mrow][4], //initial snap
									mastermatrix[mrow][5], //final snap
									masterweight[mrow][0]); 	
				
									masterweight[mrow][0]=currentweight[crow][0]; //update weight
									mastermatrix[mrow][3]=currentmatrix[crow][3]; //update flag
									mastermatrix[mrow][4]=snap; //restart time
									mastermatrix[mrow][5]=snap;
									
									mastermatrix[mrow][8]=1;//switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
								}
							}	
						}	
						else if(currentmatrix[crow][3]==4 && mastermatrix[mrow][3]==4){ //If outer H is identified								
							if(currentmatrix[crow][1]==mastermatrix[mrow][1]){ //O's matches							
								if(fabs(currentweight[crow][0] - masterweight[mrow][0]) <= prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) >= prec){ //O1H matches, update info
						
									mastermatrix[mrow][5]=snap; //update final snap
									mastermatrix[mrow][3]=currentmatrix[crow][3]; //update flag
					
									mastermatrix[mrow][8]=1; //switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
								}
								else if(fabs(currentweight[crow][0]-masterweight[mrow][0]) > prec && fabs(currentweight[crow][0]+masterweight[mrow][0]) > prec){ 
//if indices are the same, but weight changed - output persistence
							/*	printf("1z diff: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);
*/

									//Outer H weight in zundel changed
									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3],
									currentmatrix[crow][3],
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][1], //O1 index
									mastermatrix[mrow][4], //initial snap
									mastermatrix[mrow][5], //final snap
									masterweight[mrow][0]); 	
				
									masterweight[mrow][0]=currentweight[crow][0]; //update weight
									mastermatrix[mrow][3]=currentmatrix[crow][3]; //update flag
									mastermatrix[mrow][4]=snap; //restart time
									mastermatrix[mrow][5]=snap;
				
									mastermatrix[mrow][8]=1;//switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
								}
							}					
						}
						else if(currentmatrix[crow][3]==2 && mastermatrix[mrow][3]==2){ //If Zundel H is identified									
							if(currentmatrix[crow][1]==mastermatrix[mrow][1] && currentmatrix[crow][2]==mastermatrix[mrow][2]){ 
//O1 and O2 matches in zundel
								if(fabs(currentweight[crow][0] - prec) <= masterweight[mrow][0] && fabs(currentweight[crow][0] + prec) >= masterweight[mrow][0]){
									if(fabs(currentweight[crow][1] - prec) <= masterweight[mrow][1] && fabs(currentweight[crow][1] + prec) >= masterweight[mrow][1]){ 	
//O1 and O2 weight matches, update snap.
							/*										
										printf("Same Zundel %0.3f %0.3f to %0.3f %0.3f: %d %d %d %d %d\n",
										masterweight[mrow][0],
										currentweight[crow][0],
										masterweight[mrow][1],
										currentweight[crow][1],
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][2], //O index
										mastermatrix[mrow][4], //initial snap
										mastermatrix[mrow][5], //final snap
										mastermatrix[mrow][3]); //CN											
*/										
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //update flag
										mastermatrix[mrow][5]=snap; //update final snap
										mastermatrix[mrow][7]=snap; //update final snap									

										mastermatrix[mrow][8]=1; //switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix
									}								
//Same O indices, O1H has same weight and O2H is diff
									else if(fabs(currentweight[crow][1] - masterweight[mrow][1]) > prec && fabs(currentweight[crow][1] + masterweight[mrow][1]) > prec){ 
							/*			printf("1sp same O1 O2 diff: %d %d %d %d %0.3f %0.3f %d %d\n",
										snap,
										mastermatrix[mrow][0],
										mastermatrix[mrow][1],
										mastermatrix[mrow][2],
										masterweight[mrow][0],
										masterweight[mrow][1],
										mastermatrix[mrow][3],currentmatrix[crow][3]);
							*/	
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3],
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][2], //O2 index
										mastermatrix[mrow][6], //O2H time start
										mastermatrix[mrow][7], //O2H time fin
										masterweight[mrow][1]); //O2H weight			
									
										mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
										mastermatrix[mrow][2]=currentmatrix[crow][2]; //initialize O index
							    			mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][5]=snap; //update time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart time
										masterweight[mrow][0]=currentweight[crow][0]; //O1H weight
										masterweight[mrow][1]=currentweight[crow][1]; //O2H weight
									
										mastermatrix[mrow][8]=1;//switch on	
										currentmatrix[crow][3]=0; //initialize currentmatrix
									}
	
								}
								else if(fabs(currentweight[crow][1] - prec) <= masterweight[mrow][1]  && fabs(currentweight[crow][1] + prec) >= masterweight[mrow][1]){
//O2H weight is the same, but O1H is not
									if(fabs(currentweight[crow][0] - masterweight[mrow][0]) > prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) > prec){

							/*			printf("1sp same O2 O1 diff: %d %d %d %d %0.3f %0.3f %d %d\n",
										snap,
										mastermatrix[mrow][0],
										mastermatrix[mrow][1],
										mastermatrix[mrow][2],
										masterweight[mrow][0],
										masterweight[mrow][1],
										mastermatrix[mrow][3],currentmatrix[crow][3]);

							*/		 	fprintf(output,"%d %d %d %d %d %d %0.3f\n",
                                                                                mastermatrix[mrow][3],
                                                                                currentmatrix[crow][3], //new flag
                                                                                mastermatrix[mrow][0], //H index
                                                                                mastermatrix[mrow][1], //O1 index
                                                                                mastermatrix[mrow][4], //O1H time start
                                                                                mastermatrix[mrow][5], //O1H time fin
                                                                                masterweight[mrow][1]); //O2H weight

										mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
                                                                                mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
                                                                                mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
                                                                                mastermatrix[mrow][4]=snap; //update time
                                                                                mastermatrix[mrow][5]=snap; //restart time
                                                                                mastermatrix[mrow][7]=snap; //restart final snap
                                                                                masterweight[mrow][0]=currentweight[crow][0]; //O1H
                                                                                masterweight[mrow][1]=currentweight[crow][1]; //O2H

                                                                                mastermatrix[mrow][8]=1;//switch on
                                                                                currentmatrix[crow][3]=0; //initialize currentmatrix
									}

								}
								else if(fabs(currentweight[crow][0] - masterweight[mrow][0]) > prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) > prec){
									if(fabs(currentweight[crow][1] - masterweight[mrow][1]) > prec && fabs(currentweight[crow][1] + masterweight[mrow][1]) > prec){
										printf("1sp both diff: %d %d %d %d %0.3f %0.3f %d %d\n",
										snap,
										mastermatrix[mrow][0],
										mastermatrix[mrow][1],
										mastermatrix[mrow][2],
										masterweight[mrow][0],
										masterweight[mrow][1],
										mastermatrix[mrow][3],currentmatrix[crow][3]);
								
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3],
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O2 index
										mastermatrix[mrow][4], //O2H time start
										mastermatrix[mrow][5], //O2H time fin
										masterweight[mrow][0]); //O2H weight
								
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3],
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][2], //O2 index
										mastermatrix[mrow][6], //O2H time start
										mastermatrix[mrow][7], //O2H time fin
										masterweight[mrow][1]); //O2H weight			
									
										mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
										mastermatrix[mrow][2]=currentmatrix[crow][2]; //initialize O index
							    			mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap; //update time
										mastermatrix[mrow][5]=snap; //update time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart time
										masterweight[mrow][0]=currentweight[crow][0]; //O1H weight
										masterweight[mrow][1]=currentweight[crow][1]; //O2H weight
									
										mastermatrix[mrow][8]=1;//switch on	
										currentmatrix[crow][3]=0; //initialize currentmatrix
		
									}
								}

							}//if O1 and O2 matches in Zundel
					}//End Case 1: When flags are unchanged
					else if(mastermatrix[mrow][1]==0 && mastermatrix[mrow][3]==0){ //Bond reappeared!
						mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
						mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
						mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
						mastermatrix[mrow][4]=snap; //restart time
						mastermatrix[mrow][5]=snap; //restart time
						mastermatrix[mrow][6]=snap; //restart time
						mastermatrix[mrow][7]=snap; //restart final snap
						masterweight[mrow][0]=currentweight[crow][0]; //O1H
						masterweight[mrow][1]=currentweight[crow][1]; //O2H
						
						mastermatrix[mrow][8]=1;
					}
				}//end if H-index matches
}
			}//end mrow-loop
		}//end crow-loop
		
		
//*****************************************************************************
//Section VII: Output previous persistences when Eigen or Zundel changes 
//structures.
//*****************************************************************************
		
		for(crow=0; crow<clines; crow++){
			for(mrow=1; mrow<213; mrow++){
			    if(currentmatrix[crow][3]!=0){
				if(mastermatrix[mrow][8]==0){ //for unchecked bonds
					if(currentmatrix[crow][0]==mastermatrix[mrow][0]){ //If H index matches
							if(mastermatrix[mrow][3]==1 && currentmatrix[crow][3]!=2){ //If Eigen CB breaks, output persistence	
								//Bond in Eigen breaks
						/*		printf("2e change: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);

*/
								fprintf(output,"%d %d %d %d %d %d %0.3f\n",
								mastermatrix[mrow][3],
								currentmatrix[crow][3], //new flag
								mastermatrix[mrow][0], //H index
								mastermatrix[mrow][1], //O index
								mastermatrix[mrow][4], //time start
								mastermatrix[mrow][5], //final start
								masterweight[mrow][0]);
							
								mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
							        mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
								mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
								mastermatrix[mrow][4]=snap; //restart time
								mastermatrix[mrow][5]=snap; //restart time
								masterweight[mrow][0]=currentweight[crow][0]; //update weight
							
								mastermatrix[mrow][8]=1;//switch on
								currentmatrix[crow][3]=0; //initialize currentmatrix
							
							}
							else if(mastermatrix[mrow][3]==1 && currentmatrix[crow][3]==2){ //If Eigen forms Zundel
								if(mastermatrix[mrow][1]==currentmatrix[crow][1]){//If Eigen O matches in Zundel O1, examine weights
									if(fabs(currentweight[crow][0] - prec) <= masterweight[mrow][0] && fabs(currentweight[crow][0] + prec) >= masterweight[mrow][0]){
																																							//Case 1a. Eigen formed Zundel!
/*										printf("%d %d %d %d %d %d %d %d %d %0.3f %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										currentmatrix[crow][1], //O2
										mastermatrix[mrow][4], //time start
										snap,//final snap
										mastermatrix[mrow][6], //time start
										snap,//final snap									
										masterweight[mrow][0], //O1H weight
										currentweight[crow][1]);
*/ 								
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3],
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O index
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5], //final start
										masterweight[mrow][0]);

										
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap;
										mastermatrix[mrow][5]=snap; //update time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][1]=currentweight[crow][1]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix
									}
									else if(fabs(currentweight[crow][0] - masterweight[mrow][0]) > prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) > prec){
									//CB in Eigen changed and formed Zundel!
		/*								printf("Case 1b. Eigen formed Zundel! %d %d %d %d %d %d\n",
										snap,
										mastermatrix[mrow][0],
										mastermatrix[mrow][1],
										currentmatrix[crow][2],
										mastermatrix[3],currentmatrix[crow][3]);
		
*/
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight
										
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap; //restart time
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][0]=currentweight[crow][0]; //update weight
										masterweight[mrow][1]=currentweight[crow][1]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix

									}
//If Eigen O mathches, but O1H don't
								}
								else if(mastermatrix[mrow][1]==currentmatrix[crow][2]){//If Eigen O matches in Zundel O2, examine weights
									if(fabs(currentweight[crow][1] - prec) <= masterweight[mrow][0] && fabs(currentweight[crow][1] + prec) >= masterweight[mrow][0]){
//If Eigen O matches in Zundel O2, examine weights
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3],
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O index
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5], //final start
										masterweight[mrow][0]);

										//Case 1b. Eigen formed Zundel!
/*										printf("%d %d %d %d %d %d %d %d %d %0.3f %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										currentmatrix[crow][1], //O2
										mastermatrix[mrow][4], //time start
										snap,//final snap
										mastermatrix[mrow][6], //time start
										snap,//final snap									
										masterweight[mrow][0], //O1H weight
										currentweight[crow][0]); 
*/										
							    			mastermatrix[mrow][2]=currentmatrix[crow][1]; //store new O index
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap;
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][1]=currentweight[crow][0]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix
									}
									else if(fabs(currentweight[crow][0] - masterweight[mrow][1]) > prec && fabs(currentweight[crow][0] + masterweight[mrow][1]) > prec){

						/*		printf("2e to Zundel, diff weight: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);
*/

									//CB in Eigen changed and formed Zundel with O2!
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight
										
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap; //restart time
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][0]=currentweight[crow][0]; //update weight
										masterweight[mrow][1]=currentweight[crow][1]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix

									}

								}	
							}
							else if(mastermatrix[mrow][3]==2 && currentmatrix[crow][3]!=2){ //If SP CB breaks, output persistence
								if(mastermatrix[mrow][1]==currentmatrix[crow][1]){ //If O1 matches, examine weight
									if(fabs(currentweight[crow][0] - prec) <= masterweight[mrow][0] && fabs(currentweight[crow][0] + prec) >= masterweight[mrow][0]){
//Zundel dissociated and weight of O1H didn't change, update snap for bond
								/*printf("2z diss O1 same: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);

							*/

									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3],
									currentmatrix[crow][3], //new flag
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][1], //O index
									mastermatrix[mrow][4], //time start
									mastermatrix[mrow][5], //final start
									masterweight[mrow][0]);

									//O1H bond in Zundel broke
									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3], //previous flag
									currentmatrix[crow][3], //new flag
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][2], //O2 broken bond
									mastermatrix[mrow][6], //time start
									mastermatrix[mrow][7],//final snap
									masterweight[mrow][1]); //O2H weight
								
									mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
									mastermatrix[mrow][2]=0; //should be 0
							    		mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
									mastermatrix[mrow][4]=snap;
									mastermatrix[mrow][5]=snap; //update time
									mastermatrix[mrow][6]=snap; //restart time
									mastermatrix[mrow][7]=snap; //restart time
									masterweight[mrow][0]=currentweight[crow][0]; //O1H
									masterweight[mrow][1]=0.000; //initialize weight
							
									mastermatrix[mrow][8]=1;//switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
									
									}
									else if(fabs(currentweight[crow][0] - masterweight[mrow][0]) > prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) > prec){
//Zundel dissociated and weight of O1H DID change, print persistence
							/*	fprintf(output1,"2z diss O1 diff: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);
*/

									//O1H bond in Zundel broke
									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3], //previous flag
									currentmatrix[crow][3], //new flag
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][1], //O2 broken bond
									mastermatrix[mrow][4], //time start
									mastermatrix[mrow][5],//final snap
									masterweight[mrow][0]); //O2H weight


									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3], //previous flag
									currentmatrix[crow][3], //new flag
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][2], //O2 broken bond
									mastermatrix[mrow][6], //time start
									mastermatrix[mrow][7],//final snap
									masterweight[mrow][1]); //O2H weight
								
									mastermatrix[mrow][2]=0; //update O2, should be 0
							    		mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
									masterweight[mrow][4]=snap; //restart time
									mastermatrix[mrow][5]=snap; //restart time
									mastermatrix[mrow][6]=snap; //restart time
									mastermatrix[mrow][7]=snap; //restart time
									masterweight[mrow][0]=currentweight[crow][0]; //O1H
									masterweight[mrow][1]=0.000; //initialize weight
							
									mastermatrix[mrow][8]=1;//switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
									
									}

								}
								else if(mastermatrix[mrow][2]==currentmatrix[crow][1]){ //If O2 matches, examine weight
									if(fabs(currentweight[crow][0] - prec) <= masterweight[mrow][1] && fabs(currentweight[crow][0] + prec) >= masterweight[mrow][1]){
//Zundel dissociated and weight of O2H didn't change, update snap for bond
										
							/*	printf("2z diss O2 same: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);

*/
										//O2H bond in Zundel broke
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 broken bond
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight

										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][2], //O2 broken bond
										mastermatrix[mrow][6], //time start
										mastermatrix[mrow][7],//final snap
										masterweight[mrow][1]); //O2H weight


										mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
										mastermatrix[mrow][2]=0; //initialize O index
								    		mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap; //copy over time
										mastermatrix[mrow][5]=snap; //update time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart time
										masterweight[mrow][0]=currentweight[crow][0]; //copy over weight
										masterweight[mrow][1]=0.000; //O1H
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix
								
									}
									else if(fabs(currentweight[crow][0] - masterweight[mrow][1]) > prec && fabs(currentweight[crow][0] + masterweight[mrow][1]) > prec){
//Zundel dissociated and weight of O1H DID change, print persistence
							/*	printf("2z diss, O2 diff: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);

*/
									
									//O2H bond in Zundel broke
									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3], //previous flag
									currentmatrix[crow][3], //new flag
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][1], //O2 broken bond
									mastermatrix[mrow][4], //time start
									mastermatrix[mrow][5],//final snap
									masterweight[mrow][0]); //O2H weight


									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3], //previous flag
									currentmatrix[crow][3], //new flag
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][2], //O2 broken bond
									mastermatrix[mrow][6], //time start
									mastermatrix[mrow][7],//final snap
									masterweight[mrow][1]); //O2H weight

									
									mastermatrix[mrow][1]=currentmatrix[crow][1]; //copy over O index
									mastermatrix[mrow][2]=0;	
							    		mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
									masterweight[mrow][4]=snap; //restart time
									mastermatrix[mrow][5]=snap; //restart time
									mastermatrix[mrow][6]=snap; //restart time
									mastermatrix[mrow][7]=snap; //restart time
									masterweight[mrow][0]=currentweight[crow][0]; //O1H
									masterweight[mrow][1]=0.000; //initialize weight
							
									mastermatrix[mrow][8]=1;//switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
									
									}

								}
							}
							else if(mastermatrix[mrow][3]==4 && currentmatrix[crow][3]!=2){ //If Zundel CB breaks, output persistence
/*
								printf("3o change: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);


*/								//Outer bond in Zundel changed
								fprintf(output,"%d %d %d %d %d %d %0.3f\n",
								mastermatrix[mrow][3],
								currentmatrix[crow][3], //new flag
								mastermatrix[mrow][0], //H index
								mastermatrix[mrow][1], //O index
								mastermatrix[mrow][4], //time start
								mastermatrix[mrow][5], //final start
								masterweight[mrow][0]);
	
							
								mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
							    	mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
								mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
								mastermatrix[mrow][4]=snap; //restart time
								mastermatrix[mrow][5]=snap; //restart time
								mastermatrix[mrow][6]=snap; //restart time
								mastermatrix[mrow][7]=snap; //restart final snap
								masterweight[mrow][0]=currentweight[crow][0]; //O1H
								masterweight[mrow][1]=currentweight[crow][1]; //O2H
							
								mastermatrix[mrow][8]=1;//switch on
								currentmatrix[crow][3]=0; //initialize currentmatrix
							}
							else if(mastermatrix[mrow][3]==4 && currentmatrix[crow][3]==2){ //If outer H becomes shared bond
								if(mastermatrix[mrow][1]==currentmatrix[crow][1]){//If water O matches in Zundel O1, examine weights
									if(fabs(currentweight[crow][0] - prec) <= masterweight[mrow][0] && fabs(currentweight[crow][0] + prec) >= masterweight[mrow][0]){
										
						/*		printf("3o to Zundel, same O1: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1]),
								mastermatrix[mrow][3],currentmatrix[crow][3];
*/
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3],
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O index
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5], //final start
										masterweight[mrow][0]);
									
										//Bond in Zundel changed to SP O1, weight is the same
										printf("%d %d %d %d %d %d %d %d %d %0.3f %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										currentmatrix[crow][1], //O2
										mastermatrix[mrow][4], //time start
										snap,//final snap
										mastermatrix[mrow][6], //time start
										snap,//final snap									
										masterweight[mrow][0], //O1H weight
										currentweight[crow][1]); 
										
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap;
										mastermatrix[mrow][5]=snap; //update time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][0]=currentweight[crow][0];
										masterweight[mrow][1]=currentweight[crow][1]; //update weight
										
										mastermatrix[mrow][8]=1;//switch on		
										currentmatrix[crow][3]=0; //initialize currentmatrix						
									}
									else if(fabs(currentweight[crow][0] - masterweight[mrow][0]) > prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) > prec){
									//CB in Zundel changed and formed SP!

							/*	printf("3o to Zundel, diff O1: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);
*/

										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight
										
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap; //restart time
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][0]=currentweight[crow][0]; //update weight
										masterweight[mrow][1]=currentweight[crow][1]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix


									}
								}
								else if(mastermatrix[mrow][1]==currentmatrix[crow][2]){//If Zundel O matches in SP  O2, examine weights
									if(fabs(currentweight[crow][1] - prec) <= masterweight[mrow][0] && fabs(currentweight[crow][1] + prec) >= masterweight[mrow][0]){
//If Eigen O matches in Zundel O2, examine weights
							/*	printf("3o to Zundel, same O2: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],
								mastermatrix[mrow][3],currentmatrix[crow][3]);
*/

										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight


										//Zundel formed to SP O2
										printf("%d %d %d %d %d %d %d %d %d %0.3f %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										currentmatrix[crow][1], //O2
										mastermatrix[mrow][4], //time start
										snap,//final snap
										mastermatrix[mrow][6], //time start
										snap,//final snap									
										masterweight[mrow][0], //O1H weight
										currentweight[crow][0]); 
										
							    			mastermatrix[mrow][2]=currentmatrix[crow][1]; //store new O index
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap;
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][0]=currentweight[crow][1];
										masterweight[mrow][1]=currentweight[crow][0]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix
									}
									else if(fabs(currentweight[crow][1] - masterweight[mrow][0]) > prec && fabs(currentweight[crow][1] + masterweight[mrow][0]) > prec){
/*								printf("3o to Zundel, diff O2: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],mastermatrix[mrow][3],currentmatrix[crow][3]);

										
*/									//CB in Water changed and formed Zundel with O2!
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight
									
												
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap; //restart time
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][0]=currentweight[crow][0]; //update weight
										masterweight[mrow][1]=currentweight[crow][1]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix

									}

								}
							}	
							else if(mastermatrix[mrow][3]==55 && currentmatrix[crow][3]!=2){ //If wateroutput CB breaks, output persistence
							/*	printf("4w change: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],mastermatrix[mrow][3],currentmatrix[crow][3]);
*/
									//Bond in water breaks
									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
									mastermatrix[mrow][3],
									currentmatrix[crow][3], //new flag
									mastermatrix[mrow][0], //H index
									mastermatrix[mrow][1], //O index
									mastermatrix[mrow][4], //time start
									mastermatrix[mrow][5], //final start
									masterweight[mrow][0]);
							
									mastermatrix[mrow][1]=currentmatrix[crow][1]; //O1
							    		mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
									mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
									mastermatrix[mrow][4]=snap; //restart time
									mastermatrix[mrow][5]=snap; //restart time
									masterweight[mrow][0]=currentweight[crow][0]; //update weight
							
									mastermatrix[mrow][8]=1;//switch on
									currentmatrix[crow][3]=0; //initialize currentmatrix
							}
							else if(mastermatrix[mrow][3]==55 && currentmatrix[crow][3]==2){ //If water CB breaks 
								if(mastermatrix[mrow][1]==currentmatrix[crow][1]){//If water O matches in Zundel O1, examine weights
									if(fabs(currentweight[crow][0] - prec) <= masterweight[mrow][0] && fabs(currentweight[crow][0] + prec) >= masterweight[mrow][0]){
							/*	printf("4w to Zundel, same O1: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],mastermatrix[mrow][3],currentmatrix[crow][3]);
	*/
										//Bond in water changed to Zundel O1
										printf("%d %d %d %d %d %d %d %d %d %0.3f %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										currentmatrix[crow][1], //O2
										mastermatrix[mrow][4], //time start
										snap,//final snap
										mastermatrix[mrow][6], //time start
										snap,//final snap									
										masterweight[mrow][0], //O1H weight
										currentweight[crow][1]);

 										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight
	
										
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap;
										mastermatrix[mrow][5]=snap; //update time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap	
										masterweight[mrow][0]=currentweight[crow][0];
										masterweight[mrow][1]=currentweight[crow][1]; //update weight
										
										mastermatrix[mrow][8]=1;//switch on		
										currentmatrix[crow][3]=0; //initialize currentmatrix						
									}
									else if(fabs(currentweight[crow][0] - masterweight[mrow][0]) > prec && fabs(currentweight[crow][0] + masterweight[mrow][0]) > prec){
									//CB in Water changed and formed Zundel!

										
							/*	fprintf(output1,"4w to Zundel, diff O1: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],mastermatrix[mrow][3],currentmatrix[crow][3]);
*/
									fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight
										
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap; //restart time
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][0]=currentweight[crow][0]; //update weight
										masterweight[mrow][1]=currentweight[crow][1]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix


									}
								}
								else if(mastermatrix[mrow][1]==currentmatrix[crow][2]){//If water O matches in Zundel O1, examine weights
									if(fabs(currentweight[crow][1] - prec) <= masterweight[mrow][0] && fabs(currentweight[crow][1] + prec) >= masterweight[mrow][0]){
//If Eigen O matches in Zundel O2, examine weights
							/*	printf("4w to Zundel, same O2: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],mastermatrix[mrow][3],currentmatrix[crow][3]);

							*/			
										// Water formed to Zundel O2
										printf("%d %d %d %d %d %d %d %d %d %0.3f %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										currentmatrix[crow][1], //O2
										mastermatrix[mrow][4], //time start
										snap,//final snap
										mastermatrix[mrow][6], //time start
										snap,//final snap									
										masterweight[mrow][0], //O1H weight
										currentweight[crow][0]); 
									
										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight
	
	
							    			mastermatrix[mrow][2]=currentmatrix[crow][1]; //store new O index
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap;
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][1]=currentweight[crow][0]; //update weight
										masterweight[mrow][0]=currentweight[crow][1];

										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix
									}
									else if(fabs(currentweight[crow][1] - masterweight[mrow][0]) > prec && fabs(currentweight[crow][1] + masterweight[mrow][0]) > prec){
									//CB in Water changed and formed Zundel with O2!
							/*	printf("4w to Zundel, diff O2: %d %d %d %d %0.3f %0.3f %d %d\n",
								snap,
								mastermatrix[mrow][0],
								mastermatrix[mrow][1],
								mastermatrix[mrow][2],
								masterweight[mrow][0],
								masterweight[mrow][1],mastermatrix[mrow][3],currentmatrix[crow][3]);
*/

										fprintf(output,"%d %d %d %d %d %d %0.3f\n",
										mastermatrix[mrow][3], //previous flag
										currentmatrix[crow][3], //new flag
										mastermatrix[mrow][0], //H index
										mastermatrix[mrow][1], //O1 
										mastermatrix[mrow][4], //time start
										mastermatrix[mrow][5],//final snap
										masterweight[mrow][0]); //O1H weight
									
												
							    			mastermatrix[mrow][2]=currentmatrix[crow][2]; //O2
										mastermatrix[mrow][3]=currentmatrix[crow][3]; //change flag
										mastermatrix[mrow][4]=snap; //restart time
										mastermatrix[mrow][5]=snap; //restart time
										mastermatrix[mrow][6]=snap; //restart time
										mastermatrix[mrow][7]=snap; //restart final snap
										masterweight[mrow][0]=currentweight[crow][1]; //update weight
										masterweight[mrow][1]=currentweight[crow][0]; //update weight
							
										mastermatrix[mrow][8]=1;//switch on
										currentmatrix[crow][3]=0; //initialize currentmatrix

									}

								}
						}	
					} //end if H index matches
}
				}//end unchecked-bonds
			}//end mrow-loop
		}//end crow-loop

//*****************************************************************************
//Section VII: This will only occur if a bond is missing
//*****************************************************************************
		
			for(mrow=1; mrow<213; mrow++){
				if(mastermatrix[mrow][8]==0 && currentmatrix[crow][0]!=0){ //for unchecked bonds
					if(mastermatrix[mrow][3]!=2){
						printf("Bond disappeared!\n");
						
						fprintf(output,"%d 0 %d %d %d %d %0.3f\n",
						mastermatrix[mrow][3],
						mastermatrix[mrow][0],
						mastermatrix[mrow][1],
						mastermatrix[mrow][4],
						mastermatrix[mrow][5],
						masterweight[mrow][0]);
					
						mastermatrix[mrow][1]=0; //O1
				 		mastermatrix[mrow][2]=0; //O2
				 		mastermatrix[mrow][3]=0; //change flag
						mastermatrix[mrow][4]=snap; //restart time
						mastermatrix[mrow][5]=snap; //restart time
						mastermatrix[mrow][6]=snap; //restart time
						mastermatrix[mrow][7]=snap; //restart final snap
						masterweight[mrow][0]=0.000; //O1H
						masterweight[mrow][1]=0.000; //O2H
						
						mastermatrix[mrow][8]=1;//switch on	

					}
					
					else if(mastermatrix[mrow][3]==2){
						printf("Zundel bond disappeared!\n");
						fprintf(output,"%d 0 %d %d %d %d %0.3f\n",
						mastermatrix[mrow][3],
						mastermatrix[mrow][0],
						mastermatrix[mrow][1],
						mastermatrix[mrow][4],
						mastermatrix[mrow][5],
						masterweight[mrow][0]);
						
						fprintf(output,"%d 0 %d %d %d %d %0.3f\n",
						mastermatrix[mrow][3],
						mastermatrix[mrow][0],
						mastermatrix[mrow][2],
						mastermatrix[mrow][6],
						mastermatrix[mrow][7],
						masterweight[mrow][1]);
									
						mastermatrix[mrow][1]=0; //O1
				 		mastermatrix[mrow][2]=0; //O2
				 		mastermatrix[mrow][3]=0; //change flag
						mastermatrix[mrow][4]=snap; //restart time
						mastermatrix[mrow][5]=snap; //restart time
						mastermatrix[mrow][6]=snap; //restart time
						mastermatrix[mrow][7]=snap; //restart final snap
						masterweight[mrow][0]=0.000; //O1H
						masterweight[mrow][1]=0.000; //O2H
						
						mastermatrix[mrow][8]=1;//switch on	
					}	
				}
			}
	
			if(snap==maxsnap){		
				for(y=1; y<213; y++){
					if(mastermatrix[y][3]==2){ //zundels				
						fprintf(output,"%d %d %d %d %d %d %0.3f\n",
						mastermatrix[y][3],
						mastermatrix[y][3],
						mastermatrix[y][0],
						mastermatrix[y][1],
						mastermatrix[y][4],
						mastermatrix[y][5],
						masterweight[y][0]);
						
						fprintf(output,"%d %d %d %d %d %d %0.3f\n",
						mastermatrix[y][3],
						mastermatrix[y][3],
						mastermatrix[y][0],
						mastermatrix[y][2],
						mastermatrix[y][6],
						mastermatrix[y][7],
						masterweight[y][1]);
					}
					else if(mastermatrix[y][3]!=2){ //for nonzundels
						fprintf(output,"%d %d %d %d %d %d %0.3f\n",
						mastermatrix[y][3],
						mastermatrix[y][3],
						mastermatrix[y][0],
						mastermatrix[y][1],
						mastermatrix[y][4],
						mastermatrix[y][5],
						masterweight[y][0]);
					}
				}
			}

//Delete currentmatrix matrix for next snapshot
	      for(y=0; y<clines; y++){
	 	    for(z=0; z<7; z++){
	 	       currentmatrix[y][z]=0;
	 	       currentweight[y][z]=0.000;
	 	     }
	      }//end y-loop
		  
	      for(y=0; y<lines; y++){
	          for(z=0; z<6; z++){
			labelmatrix[y][z]=0;
		   }
	      }

		  for(y=1; y<213; y++){
			  mastermatrix[y][8]=0;
		  }
	      clines=0;
		  lines=0;
	      countHindex=1;
		  countOindex=1; 
		  
	  }//end if snap<1
   }//end snap-loop
   fclose(output);

	
}//end main-loop
