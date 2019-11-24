#include<stdio.h>
#include<stdlib.h>

//*****************************************************************
//Program written to sort updated GraphGeod files in a single 
//snapshot and bin them according to angle (increments of 
//2 degrees) and delta (increments of 0.1A).
//
//To obtain deltaOH, bin-diff-updated-GG.sh is required!
//
//Output: 
//1st column: angle range
//1st row: delta value
//2nd - 12th column: count of structures with specific angle and
//delta value	  
//
//Created by Lelee Ounkham
//Last Updated: September 6, 2017
//*****************************************************************
int main ()
{
//Setting up parameters to read in GG file
	int O1,O2,H,lines;
	float OOdist,O1H,O2H,ang,diff;
	float angle,newangle;
	char ifile[100],ofile1[100];
	FILE *input,*output;

//Parameters for diffarray
	int maxlines=175000;
	float diffarray[maxlines][3];
	int count[13];
	int printerarray[16][13];

//Counter parameters
	int i,j,p,y,z,row;

	for(y=0; y<16; y++){
	   for(p=0; p<13; p++){
	      count[p]=0;
 	      printerarray[y][p]=0.0;
	   }
	}

	z=1;
	for(angle=150; angle <= 180; (angle = angle+5)){
	   printerarray[z][0]=angle;//first column of array corresponds to angle
	   z++;
	}

//for(j=1; j<16; j++){
//printf("%f\n",printerarray[j][0]);
//}
//
//Zero out parameters

	   diffarray[lines][0]=0.0;
	   diffarray[lines][1]=0.0;
	   diffarray[lines][2]=0.0;
	   lines=0;
	   j=0;	

//	for(i=1; i<=4; i++){
//         sprintf(ifile,"OH-diff%d.GG",i);
	   input=fopen("total-2m-anglevsdelta.txt","r");
//           output=fopen("total-angle-vs-delta.txt","a");
	 
//Storing OO distances and O-H distances into array name diffarray	   
	  while(fscanf(input,"%d %d %d %f %f %f %f %f\n",&O1,&O2,&H,&OOdist,&O1H,&O2H,&ang,&diff)==8){	     
	     diffarray[lines][0]=OOdist;
	     diffarray[lines][1]=diff;
	     diffarray[lines][2]=ang;
	     lines++;
	  }//end while
	  fclose(input);

	        for(angle=150.0; angle < 180.0; (angle = angle+5)){
	           newangle=angle+5.0; //increments of 2.0 degrees
		   j++; //increment j for printerarray every angle range
		   for(row=0; row<lines; row++){
                      if(diffarray[row][2]!=0.0){
		         if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50){ //if delta is within a range
			    if(diffarray[row][2] > angle && diffarray[row][2] < newangle){
			             count[1]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][1]=count[1];
                            }
			 }
                         else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40){
		            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[2]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][2]=count[2];
                            }//end if
		         }
                         else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30){
        	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[3]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][3]=count[3];
		            }//end if
		         }
                         else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20){
        	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[4]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][4]=count[4];
                            }//end if
		         }
                         else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10){
         	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[5]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][5]=count[5];
		            }//end if
		         }
                         else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00){
         	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[6]++;
 				     diffarray[row][2]=0.0;
				     printerarray[j][6]=count[6];
                            }//end if
		         }
                         else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10){
         	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[7]++;
 				     diffarray[row][2]=0.0;
				     printerarray[j][7]=count[7];
                            }//end if
		         }
                         else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20){
		            if(diffarray[row][2] > angle && diffarray[row][2] < newangle){
                                     count[8]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][8]=count[8];
		      	    }
                         }//end if
                         else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30){
        	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[9]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][9]=count[9];
                            }//end if
		         }
                         else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40){
        	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[10]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][10]=count[10];
                            }//end if
		         }
                         else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50){
        	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[11]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][11]=count[11];
                            }//end if
		         }
                         else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60){
        	            if(diffarray[row][2] > angle  && diffarray[row][2] < newangle){
                                     count[12]++;
				     diffarray[row][2]=0.0;
				     printerarray[j][12]=count[12];
                            }//end if 
		         }

	           }//end diffarray not 0
    	        }//end row loop
	        for(y=0; y<13; y++){
		count[y]=0;
	        }
             }//end angle-loop
//       }//end i-loop
 
       for(p=0; p<=j; p++){
          printf("%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
	  printerarray[p][0],
	  printerarray[p][1],
	  printerarray[p][2],
	  printerarray[p][3],
	  printerarray[p][4],
	  printerarray[p][5],
	  printerarray[p][6],
	  printerarray[p][7],	
	  printerarray[p][8],
	  printerarray[p][9],
	  printerarray[p][10],
	  printerarray[p][11],
	  printerarray[p][12]);
       }//end p-loop

  //     fclose(output);

}//end main-loop
