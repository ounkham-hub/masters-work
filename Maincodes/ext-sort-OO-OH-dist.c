#include<stdio.h>
#include<stdlib.h>

//*****************************************************************
//Program written to sort updated GraphGeod files in a single 
//snapshot and bin them according to a range of values and 
//count the frequency per bin.
//
//To obtain updated GG, post-select-angle.c must be run prior to.
//
//4 outputs are created:
//map1: O-O distance from 2.3-2.4A, angle from 150 - 180
//map2: O-O distance from 2.4-2.5A, angle from 150 - 180
//map3: O-O distance from 2.5-2.6A, angle from 150 - 180
//map4: O-O distance from 2.6-2.7A, angle from 150 - 180

//Created by Lelee Ounkham
//Last Updated: March 06, 2017
//*****************************************************************
int main ()
{
//Setting up parameters to read in GG file
	int O1,O2,H;
	float OOdist,O1H,O2H,ang,diff;
	char ifile[100],ofile1[100],ofile2[100],ofile3[100],ofile4[100];
	FILE *input,*output1,*output2,*output3,*output4;

//Parameters for diffarray
	int l=500;
	int i,n,x;
	float diffarray[l][3];
//Counter parameters
	int row,coln;
    int count1,count2,count3,count4,count5,count6;
	int count7,count8,count9,count10,count11,count12;
	int count13,count14,count15,count16,count17,count18;

	for(i=1; i<=61556; i++)
   	{
  	   sprintf(ifile,"HCl.input.solB%d.xyz.solB%d.xyz.GraphGeod-angle",i,i);
 	   sprintf(ofile1,"angle1.txt",i);
	   sprintf(ofile2,"angle2.txt",i);
	   sprintf(ofile3,"angle3.txt",i);
	   sprintf(ofile4,"angle4.txt",i);
        
	   input=fopen(ifile,"r");
        
	   output1=fopen(ofile1,"a");
	   output2=fopen(ofile2,"a");
	   output3=fopen(ofile3,"a");
	   output4=fopen(ofile4,"a");
	   
//Zero out parameters
	   n=0;
	   count1=0;
	   count2=0;
	   count3=0;
	   count4=0;
	   count5=0;
	   count6=0;
	   count7=0;
	   count8=0;
	   count9=0;
	   count10=0;
	   count11=0;
	   count12=0;
	   count13=0;
	   count14=0;
	   count15=0;
	   count16=0;
	   count17=0;
	   count18=0;
	   diffarray[n][0]=0;
	   diffarray[n][1]=0;
	   diffarray[n][2]=0;
	   
	   //Storing OO distances and O-H distances into array name diffarray
	   
	  while(fscanf(input,"%d %d %d %f %f %f %f\n",&O1,&O2,&H,&OOdist,&O1H,&O2H,&ang)==7)
	  {	     
	     diffarray[n][0]=OOdist;
	     diffarray[n][1]=O1H;
	     diffarray[n][2]=ang;
	     n++;
	  }//end while

	  fclose(input);
	   
 	  for(row=0; row<n; row++)
 	  {	
 	        if(diffarray[row][0] >= 2.3 && diffarray[row][0] <= 2.4)
 	        {
 	          if(diffarray[row][2] >= 150.0 && diffarray[row][2] <= 160.0)
 	          {
 				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
 				  {
 					  count1++;
 				  }
 				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
 				  {
 					  count2++;
 				  }
 				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
 				  {
 					  count3++;
 				  }
 				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
 				  {
 					  count4++;
 				  }
 				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
 				  {
 					  count5++;
 				  }
 				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
 				  {
 					  count6++;
 				  }
 	          }//end ang-if
 	          else if(diffarray[row][2] >= 160.0 && diffarray[row][2] <= 170.0)
 	          {
 				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
 				  {
 					  count7++;
 				  }
 				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
 				  {
 					  count8++;
 				  }
 				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
 				  {
 					  count9++;
 				  }
 				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
 				  {
 					  count10++;
 				  }
 				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
 				  {
 					  count11++;
 				  }
 				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
 				  {
 					  count12++;
 				  }
 	          }//end ang-if
 	          else if(diffarray[row][2] >= 170.0 && diffarray[row][2] <= 180.0)
 	          {
 				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
 				  {
 					  count13++;
 				  }
 				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
 				  {
 					  count14++;
 				  }
 				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
 				  {
 					  count15++;
 				  }
 				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
 				  {
 					  count16++;
 				  }
 				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
 				  {
 					  count17++;
 				  }
 				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
 				  {
 					  count18++;
 				  }
 	          }//end if
 	        }//end 2.3-2.4-if-loop
 	   }//end row-loop
 	    fprintf(output1,"%d %d %d %d %d %d   %d %d %d %d %d %d   %d %d %d %d %d %d\n",
 		count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,
 		count11,count12,count13,count14,count15,count16,count17,count18);
 		fclose(output1);
	         
//zero out counters here
	   count1=0;
	   count2=0;
	   count3=0;
	   count4=0;
	   count5=0;
	   count6=0;
	   count7=0;
	   count8=0;
	   count9=0;
	   count10=0;
	   count11=0;
	   count12=0;
	   count13=0;
	   count14=0;
	   count15=0;
	   count16-0;
	   count17=0;
	   count18=0; 
	   
  	  for(row=0; row<n; row++)
  	  {	
  	        if(diffarray[row][0] >= 2.4 && diffarray[row][0] <= 2.5)
  	        {
  	          if(diffarray[row][2] >= 150.0 && diffarray[row][2] <= 160.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count1++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count2++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count3++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count4++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count5++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count6++;
  				  }
  	          }//end ang-if
  	          else if(diffarray[row][2] >= 160.0 && diffarray[row][2] <= 170.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count7++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count8++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count9++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count10++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count11++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count12++;
  				  }
  	          }//end ang-if
  	          else if(diffarray[row][2] >= 170.0 && diffarray[row][2] <= 180.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count13++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count14++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count15++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count16++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count17++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count18++;
  				  }
  	          }//end if
  	        }//end 2.4-2.5-if-loop
  	   }//end row-loop
  	    fprintf(output2,"%d %d %d %d %d %d   %d %d %d %d %d %d   %d %d %d %d %d %d\n",
  		count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,
  		count11,count12,count13,count14,count15,count16,count17,count18);
  		fclose(output2);
		
 ////zero out counters here
 	   count1=0;
 	   count2=0;
 	   count3=0;
 	   count4=0;
 	   count5=0;
 	   count6=0;
 	   count7=0;
 	   count8=0;
 	   count9=0;
 	   count10=0;
 	   count11=0;
 	   count12=0;
	   count13=0;
	   count14=0;
	   count15=0;
	   count16-0;
	   count17=0;
	   count18=0; 
	   
  	  for(row=0; row<n; row++)
  	  {	
  	        if(diffarray[row][0] >= 2.5 && diffarray[row][0] <= 2.6)
  	        {
  	          if(diffarray[row][2] >= 150.0 && diffarray[row][2] <= 160.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count1++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count2++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count3++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count4++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count5++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count6++;
  				  }
  	          }//end ang-if
  	          else if(diffarray[row][2] >= 160.0 && diffarray[row][2] <= 170.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count7++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count8++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count9++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count10++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count11++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count12++;
  				  }
  	          }//end ang-if
  	          else if(diffarray[row][2] >= 170.0 && diffarray[row][2] <= 180.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count13++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count14++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count15++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count16++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count17++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count18++;
  				  }
  	          }//end if
  	        }//end 2.5-2.6-if-loop
  	   }//end row-loop
  	    fprintf(output3,"%d %d %d %d %d %d   %d %d %d %d %d %d   %d %d %d %d %d %d\n",
  		count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,
  		count11,count12,count13,count14,count15,count16,count17,count18);
  		fclose(output3);
		
 ////zero out counters here
 	   count1=0;
 	   count2=0;
 	   count3=0;
 	   count4=0;
 	   count5=0;
 	   count6=0;
 	   count7=0;
 	   count8=0;
 	   count9=0;
 	   count10=0;
 	   count11=0;
 	   count12=0;
	   count13=0;
	   count14=0;
	   count15=0;
	   count16-0;
	   count17=0;
	   count18=0; 
	   
  	  for(row=0; row<n; row++)
  	  {	
  	        if(diffarray[row][0] >= 2.6 && diffarray[row][0] <= 2.7)
  	        {
  	          if(diffarray[row][2] >= 150.0 && diffarray[row][2] <= 160.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count1++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count2++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count3++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count4++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count5++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count6++;
  				  }
  	          }//end ang-if
  	          else if(diffarray[row][2] >= 160.0 && diffarray[row][2] <= 170.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count7++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count8++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count9++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count10++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count11++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count12++;
  				  }
  	          }//end ang-if
  	          else if(diffarray[row][2] >= 170.0 && diffarray[row][2] <= 180.0)
  	          {
  				  if(diffarray[row][1] >= 0.5 && diffarray[row][1] <= 0.75)
  				  {
  					  count13++;
  				  }
  				  if(diffarray[row][1] > 0.75 && diffarray[row][1] <= 1.0)
  				  {
  					  count14++;
  				  }
  				  if(diffarray[row][1] > 1.0 && diffarray[row][1] <= 1.25)
  				  {
  					  count15++;
  				  }
  				  if(diffarray[row][1] > 1.25 && diffarray[row][1] <= 1.50)
  				  {
  					  count16++;
  				  }
  				  if(diffarray[row][1] > 1.50 && diffarray[row][1] <= 1.75)
  				  {
  					  count17++;
  				  }
  				  if(diffarray[row][1] > 1.75 && diffarray[row][1] <= 2.0)
  				  {
  					  count18++;
  				  }
  	          }//end if
  	        }//end 2.6-2.7-if-loop
  	   }//end row-loop
  	    fprintf(output4,"%d %d %d %d %d %d   %d %d %d %d %d %d   %d %d %d %d %d %d\n",
  		count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,
  		count11,count12,count13,count14,count15,count16,count17,count18);
  		fclose(output4);
		
 ////zero out counters here
 	   count1=0;
 	   count2=0;
 	   count3=0;
 	   count4=0;
 	   count5=0;
 	   count6=0;
 	   count7=0;
 	   count8=0;
 	   count9=0;
 	   count10=0;
 	   count11=0;
 	   count12=0;
	   count13=0;
	   count14=0;
	   count15=0;
	   count16-0;
	   count17=0;
	   count18=0; 

	}//end i-loop
    
} 
