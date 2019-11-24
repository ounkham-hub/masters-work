#include<stdio.h>
#include<stdlib.h>

//*****************************************************************
//Program written to sort updated GraphGeod files in a single 
//snapshot and bin them according to a range of values and 
//count the frequency per bin.
//
//To obtain deltaOH, bin-diff-updated-GG.sh is required!
//
//4 outputs are created:
//Map 1: O-O distance from 2.3-2.4A,deltaOH from -2 to +2
//Map 2: O-O distance from 2.4-2.5A,deltaOH from -2 to +2
//Map 3: O-O distance from 2.5-2.6A,deltaOH from -2 to +2
//Map 4: O-O distance from 2.6-2.7A,deltaOH from -2 to +2	
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
	FILE *input,*output1,*output2,*output3,*output4;
//Parameters for diffarray
	int l=300;
	int i,n,x;
	float diffarray[l][2];
//Counter parameters
	int row,coln;
	int countmin,count1,count2,count3,count4,zundel1,zundel2,countmax;

	for(i=33; i<=33; i++)
   	{
           sprintf(ifile,"OH-diff%d.GG",i);
	   sprintf(ofile1,"map1-%d",i);
 	   sprintf(ofile2,"map2-%d",i);
	   sprintf(ofile3,"map3-%d",i);
           sprintf(ofile4,"map4-%d",i);
	   input=fopen(ifile,"r");
           output1=fopen(ofile1,"w");
	   output2=fopen(ofile2,"w");
	   output3=fopen(ofile3,"w");
	   output4=fopen(ofile4,"w");
		   
//Zero out parameters
	   n=0;
	   countmin=0;
	   count1=0;
	   count2=0;
	   count3=0;
	   count4=0;
	   zundel1=0;
	   zundel2=0;
	   countmax=0;
	   diffarray[n][0]=0;
	   diffarray[n][1]=0;
	   //Storing OO distances and O-H distances into array name diffarray
	   
	  while(fscanf(input,"%d %d %d %f %f %f %f %f\n",&O1,&O2,&H,&OOdist,&O1H,&O2H,&ang,&diff)==8)
	  {	     
	     diffarray[n][0]=OOdist;
	     diffarray[n][1]=diff;
	     n++;
	  }//end while

	  fclose(input);

	  for(row=0; row<n; row++)
	  {	
	        if(diffarray[row][0] >= 2.3 && diffarray[row][0] <= 2.4)
	        {
	          if(diffarray[row][1] >= -2.0 && diffarray[row][1] <= -0.60)
	          {
	            countmin++;
	          }//end if
	          else if(diffarray[row][1] >= -0.60 && diffarray[row][1] <= -0.40)
	          {
	            count1++;
	          }//end if
	          else if(diffarray[row][1] >= -0.40 && diffarray[row][1] <= -0.20)
	          {
	            count2++;
	          }//end if
	          else if(diffarray[row][1] >= -0.20 && diffarray[row][1] <= 0.00)
	          {
	            zundel1++;
	          }//end if
	          else if(diffarray[row][1] >= 0.00 && diffarray[row][1] <= 0.20)
	          {
	            zundel2++;
	          }//end if
	          else if(diffarray[row][1] >= 0.20 && diffarray[row][1] <= 0.40)
	          {
	            count3++;
	          }//end if
	          else if(diffarray[row][1] >= 0.40 && diffarray[row][1] <= 0.60)
	          {
	            count4++;
	          }//end if
                  else if(diffarray[row][1] >= 0.60 && diffarray[row][1] <= 2.0)
	          {    
	            countmax++;
	          }//end if
	        }//end 2.3-2.4-for-loop
	      }//end for-loop
	        fprintf(output1,"%d %d %d %d %d %d %d %d\n",
	        countmin,count1,count2,zundel1,zundel2,count3,count4,countmax);
	
		fclose(output1);
//zero out counters here
		countmin=0;
		count1=0;
	 	count2=0;
		count3=0;
		count4=0;
		countmax=0;
		zundel1=0;
		zundel2=0;

	     for(row=0; row<n; row++)
	     {
	      if(diffarray[row][0] >= 2.4 && diffarray[row][0] <= 2.5)
	      {
	        if(diffarray[row][1] >= -2.0 && diffarray[row][1] <= -0.60)
	        {
	          countmin++;
	        }//end if
	        else if(diffarray[row][1] >= -0.60 && diffarray[row][1] <= -0.40)
	        {
	          count1++;
	        }//end if
	        else if(diffarray[row][1] >= -0.40 && diffarray[row][1] <= -0.20)
	        {
	          count2++;
	        }//end if
	        else if(diffarray[row][1] >= -0.20 && diffarray[row][1] <= 0.00)
	        {
	          zundel1++;
	        }//end if
	        else if(diffarray[row][1] >= 0.00 && diffarray[row][1] <= 0.20)
	        {
	          zundel2++;
	        }//end if
	        else if(diffarray[row][1] >= 0.20 && diffarray[row][1] <= 0.40)
	        {
	          count3++;
	        }//end if
	        else if(diffarray[row][1] >= 0.40 && diffarray[row][1] <= 0.60)
	        {
	          count4++;
	        }//end if
                else if(diffarray[row][1] >= 0.60 && diffarray[row][1] <= 2.0)
	        {    
	          countmax++;
		}//end if
	       }//end 2.4-2.5-for-loop
	      }//end for-loop	   
        fprintf(output2,"%d %d %d %d %d %d %d %d\n",
        countmin,count1,count2,zundel1,zundel2,count3,count4,countmax);
	
	fclose(output2);	
//zero out counters here
   		  countmin=0;
   		  count1=0;
   	 	  count2=0;
   		  count3=0;
   		  count4=0;
   		  countmax=0;
   	   	  zundel1=0;
   		  zundel2=0;
	 
	     for(row=0; row<n; row++)
	     {
	      if(diffarray[row][0] >= 2.5 && diffarray[row][0] <= 2.6)
	      {
	        if(diffarray[row][1] >= -2.0 && diffarray[row][1] <= -0.60)
	        {
	          countmin++;
	        }//end if
	        else if(diffarray[row][1] >= -0.60 && diffarray[row][1] <= -0.40)
	        {
	          count1++;
	        }//end if
	        else if(diffarray[row][1] >= -0.40 && diffarray[row][1] <= -0.20)
	        {
	          count2++;
	        }//end if
	        else if(diffarray[row][1] >= -0.20 && diffarray[row][1] <= 0.00)
	        {
	          zundel1++;
	        }//end if
	        else if(diffarray[row][1] >= 0.00 && diffarray[row][1] <= 0.20)
	        {
	          zundel2++;
	        }//end if
	        else if(diffarray[row][1] >= 0.20 && diffarray[row][1] <= 0.40)
	        {
	          count3++;
	        }//end if
	        else if(diffarray[row][1] >= 0.40 && diffarray[row][1] <= 0.60)
	        {
	          count4++;
	        }//end if
                else if(diffarray[row][1] >= 0.60 && diffarray[row][1] <= 2.0)
	        {    
	          countmax++;
		    }//end if
	       }//end 2.5-2.6-for-loop
	      }//end for-loop	   
           fprintf(output3,"%d %d %d %d %d %d %d %d\n",
           countmin,count1,count2,zundel1,zundel2,count3,count4,countmax);
		
	   fclose(output3);
  //zero out counters here
      		  countmin=0;
      		  count1=0;
      	 	  count2=0;
      		  count3=0;
      		  count4=0;
      		  countmax=0;
      	   	  zundel1=0;
      		  zundel2=0;

	      for(row=0; row<n; row++)
 	      {
	       if(diffarray[row][0] >= 2.6 && diffarray[row][0] <= 2.7)
	       {
	         if(diffarray[row][1] >= -2.0 && diffarray[row][1] <= -0.60)
	         {
		   countmin++;
	         }//end if
	         else if(diffarray[row][1] >= -0.60 && diffarray[row][1] <= -0.40)
	         {
		   count1++;
	         }//end if
	         else if(diffarray[row][1] >= -0.40 && diffarray[row][1] <= -0.20)
	         {
	           count2++;
                 }//end if
	         else if(diffarray[row][1] >= -0.20 && diffarray[row][1] <= 0.00)
	         {
	           zundel1++;
	         }//end if
	         else if(diffarray[row][1] >= 0.00 && diffarray[row][1] <= 0.20)
	         {
                   zundel2++;
	         }//end if
	         else if(diffarray[row][1] >= 0.20 && diffarray[row][1] <= 0.40)
	         {
		   count3++;
	         }//end if
	         else if(diffarray[row][1] >= 0.40 && diffarray[row][1] <= 0.60)
	         {
		   count4++;
	         }//end if
	         else if(diffarray[row][1] >= 0.60 && diffarray[row][1] <= 2.0)
	         {      
		   countmax++;
	         }//end if
	        }//end 2.6-2.7-for-loop
	       }//end for	   
	       fprintf(output4,"%d %d %d %d %d %d %d %d\n",
	       countmin,count1,count2,zundel1,zundel2,count3,count4,countmax);
		   
	       fclose(output4);
	    }//end row-loop
	
/*	  for(x=0; x<n; x++)
          {
             printf("%f\n",diffarray[x]);
          }//end for
*/	  

//	}//end main for

}
