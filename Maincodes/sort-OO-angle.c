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
//map1: O-O distance from 2.3-2.35A, angle from 150 - 180
//map2: O-O distance from 2.35-2.4A, angle from 150 - 180
//map3: O-O distance from 2.4-2.45A, angle from 150 - 180
//map4: O-O distance from 2.45-2.5A, angle from 150 - 180
//map5: O-O distance from 2.5-2.55A, angle from 150 - 180
//map6: O-O distance from 2.55-2.6A, angle from 150 - 180
//map7: O-O distance from 2.6-2.65A, angle from 150 - 180
//map8: O-O distance from 2.65-2.7A, angle from 150 - 180
//
//Created by Lelee Ounkham
//Last Updated: January 22, 2017
//*****************************************************************
int main ()
{
//Setting up parameters to read in GG file
	int O1,O2,H;
	float OOdist,O1H,O2H,ang,diff;
	char ifile[100],ofile1[100],ofile2[100],ofile3[100],ofile4[100];
	char ofile5[100],ofile6[100],ofile7[100],ofile8[100],ofile9[100];
    	char ofile10[100];
	FILE *input,*output1,*output2,*output3,*output4;
	FILE *output5,*output6,*output7,*output8,*output9,*output10;

//Parameters for diffarray
	int l=500;
	int i,n,x;
	float diffarray[l][2];
//Counter parameters
	int row,coln;
    int count1,count2,count3;

	for(i=1; i<=34004; i++)
   	{
  	   sprintf(ifile,"OH-diff%d.GG",i);
	   sprintf(ofile1,"map2.0-2.1",i);
 	   sprintf(ofile2,"map2.1-2.2",i);
	   sprintf(ofile3,"map2.2-2.3",i);
           sprintf(ofile4,"map2.3-2.4",i);
	   sprintf(ofile5,"map2.4-2.5",i);
	   sprintf(ofile6,"map2.5-2.6",i);
	   sprintf(ofile7,"map2.6-2.7",i);
           sprintf(ofile8,"map2.7-2.8",i);
           sprintf(ofile9,"map2.8-2.9",i);
           sprintf(ofile10,"map2.9-3.0",i);
        
	   input=fopen(ifile,"r");
        
       	   output1=fopen(ofile1,"a");
	   output2=fopen(ofile2,"a");
	   output3=fopen(ofile3,"a");
	   output4=fopen(ofile4,"a");
	   output5=fopen(ofile5,"a");
	   output6=fopen(ofile6,"a");
	   output7=fopen(ofile7,"a");
	   output8=fopen(ofile8,"a");
           output9=fopen(ofile9,"a");
	   output10=fopen(ofile10,"a");
	   
//Zero out parameters
	   n=0;
	   count1=0;
	   count2=0;
	   count3=0;
	   diffarray[n][0]=0;
	   diffarray[n][1]=0;
	   //Storing OO distances and O-H distances into array name diffarray
	   
	  while(fscanf(input,"%d %d %d %f %f %f %f %f\n",&O1,&O2,&H,&OOdist,&O1H,&O2H,&ang,&diff)==8)
	  {	     
	     diffarray[n][0]=OOdist;
	     diffarray[n][1]=ang;
	     n++;
	  }//end while

	  fclose(input);

	  for(row=0; row<n; row++)
	  {	
	        if(diffarray[row][0] >= 2.0 && diffarray[row][0] <= 2.1)
	        {
	          if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
	          {
	            count1++;
	          }//end if
	          else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
	          {
	            count2++;
	          }//end if
	          else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
	          {
	            count3++;
	          }//end if
	        }//end 2.0-2.1-for-loop
	   }//end for-loop
	        fprintf(output1,"%d %d %d\n",count1,count2,count3);
		fclose(output1);
	         
//zero out counters here
        count1=0;
        count2=0;
        count3=0;
        
        for(row=0; row<n; row++)
        {
            if(diffarray[row][0] >= 2.1 && diffarray[row][0] <= 2.2)
            {
                if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
                {
                    count1++;
                }//end if
                else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
                {
                    count2++;
                }//end if
                else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
                {
                    count3++;
                }//end if
            }//end 2.1-2.2-for-loop
        }//end for-loop
        fprintf(output2,"%d %d %d\n",count1,count2,count3);
        
        fclose(output2);
      
//zero out counters here
        count1=0;
        count2=0;
        count3=0;
        
        for(row=0; row<n; row++)
        {
            if(diffarray[row][0] >= 2.2 && diffarray[row][0] <= 2.3)
            {
                if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
                {
                    count1++;
                }//end if
                else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
                {
                    count2++;
                }//end if
                else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
                {
                    count3++;
                }//end if
            }//end 2.2-2.3-for-loop
        }//end for-loop
        fprintf(output3,"%d %d %d\n",count1,count2,count3);
        
        fclose(output3);

//zero out counters here
        count1=0;
        count2=0;
        count3=0;

        for(row=0; row<n; row++)
        {
            if(diffarray[row][0] >= 2.3 && diffarray[row][0] <= 2.4)
            {
                if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
                {
                    count1++;
                }//end if
                else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
                {
                    count2++;
                }//end if
                else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
                {
                    count3++;
                }//end if
            }//end 2.2-2.3-for-loop
        }//end for-loop
        fprintf(output4,"%d %d %d\n",count1,count2,count3);        
        fclose(output4);
        
//zero out counters here
        count1=0;
        count2=0;
        count3=0;
        
        for(row=0; row<n; row++)
        {
            if(diffarray[row][0] >= 2.4 && diffarray[row][0] <= 2.5)
            {
                if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
                {
                    count1++;
                }//end if
                else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
                {
                    count2++;
                }//end if
                else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
                {
                    count3++;
                }//end if
            }//end 2.3-2.4-for-loop
        }//end for-loop
        fprintf(output5,"%d %d %d\n",count1,count2,count3);    
        fclose(output5);
        
//zero out counters here
        count1=0;
        count2=0;
        count3=0;
  
        for(row=0; row<n; row++)
        {
            if(diffarray[row][0] >= 2.5 && diffarray[row][0] <= 2.6)
            {
                if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
                {
                    count1++;
                }//end if
                else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
                {
                    count2++;
                }//end if
                else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
                {
                    count3++;
                }//end if
            }//end 2.4-2.5-for-loop
        }//end for-loop
        fprintf(output6,"%d %d %d\n",count1,count2,count3);        
        fclose(output6);
       
//zero out counters here
        count1=0;
        count2=0;
        count3=0;
        
        for(row=0; row<n; row++)
        {
            if(diffarray[row][0] >= 2.6 && diffarray[row][0] <= 2.7)
            {
                if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
                {
                    count1++;
                }//end if
                else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
                {
                    count2++;
                }//end if
                else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
                {
                    count3++;
                }//end if
            }//end 2.6-2.7-for-loop
        }//end for-loop
        fprintf(output7,"%d %d %d\n",count1,count2,count3);      
        fclose(output7);

//zero out counters here
        count1=0;
        count2=0;
        count3=0;
        
        for(row=0; row<n; row++)
        {
            if(diffarray[row][0] >= 2.7 && diffarray[row][0] <= 2.8)
            {
                if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
                {
                    count1++;
                }//end if
                else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
                {
                    count2++;
                }//end if
                else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
                {
                    count3++;
                }//end if
            }//end 2.7-2.8-for-loop
        }//end for-loop
        fprintf(output8,"%d %d %d\n",count1,count2,count3);        
        fclose(output8);
        
//zero out counters here
        count1=0;
        count2=0;
        count3=0;
        for(row=0; row<n; row++)
        {
	    if(diffarray[row][0] >= 2.7 && diffarray[row][0] <= 2.8)
            {
            	if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
            	{
                	count1++;
            	}//end if
            	else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
            	{
                	count2++;
            	}//end if
            	else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
            	{
                	count3++;
            	}//end if
           }//end 2.8-2.9-for-loop
        }//end for-loop
        fprintf(output9,"%d %d %d\n",count1,count2,count3);        
        fclose(output9);
        
//zero out counters here
        count1=0;
        count2=0;
        count3=0;
        for(row=0; row<n; row++)
        {
            if(diffarray[row][0] >= 2.9 && diffarray[row][0] <= 3.0)
            {
                if(diffarray[row][1] >= 150.0 && diffarray[row][1] <= 160.0)
                {
                    count1++;
                }//end if
                else if(diffarray[row][1] >= 160.0 && diffarray[row][1] <= 170.0)
                {
                    count2++;
                }//end if
                else if(diffarray[row][1] >= 170.0 && diffarray[row][1] <= 180.0)
                {
                    count3++;
                }//end if
            }//end 2.9-3.0-for-loop
        }//end for-loop
        fprintf(output10,"%d %d %d\n",count1,count2,count3);        
        fclose(output10);
	
	}//end i-loop
    
} 
