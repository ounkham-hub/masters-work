#include<stdio.h>
#include<stdlib.h>

//*****************************************************************
//Program written to sort updated GraphGeod files in a single 
//snapshot and bin them according to a range of deltaOH, -2 to +2
//
//Requires using bin-diff-GG.sh to obtain the OH-diff and add
//to updated GraphGeod files!
//
//Created by Lelee Ounkham
//Last Updated: January 22, 2017
//*****************************************************************
int main ()
{
//Setting up parameters to read in GG file
	int O1,O2,H;
	float OOdist,O1H,O2H,ang,diff;
	char ifile[100],ofile[100];
	FILE *input,*output;
//Parameters for diffarray
	int l=150;
	int i,n,x;
	float diffarray[l];
//Counter parameters
	int row,countmin,count1,count2,count3,count4,zundel1,zundel2,countmax;

	for(i=1; i<=34004; i++)
   	{
           sprintf(ifile,"OH-diff%d.GG",i);
           sprintf(ofile,"diff-count%d.txt",i);
           input=fopen(ifile,"r");
           output=fopen(ofile,"w");
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
	   diffarray[n]=0;
 
	  while(fscanf(input,"%d %d %d %f %f %f %f %f\n",&O1,&O2,&H,&OOdist,&O1H,&O2H,&ang,&diff)==8)
	  {	     
	     diffarray[n]=diff;
	     n++;
	  }//end while

	  fclose(input);

	  for(row=0; row<n; row++)
	  {	
	  
	    if(diffarray[row] >= -2.0 && diffarray[row] <= -0.60)
	    {
	      countmin++;
	    }//end if
	
	    else if(diffarray[row] >= -0.60 && diffarray[row] <= -0.40)
	    {
	      count1++;
	    }//end if

	    else if(diffarray[row] >= -0.40 && diffarray[row] <= -0.20)
	    {
	      count2++;
	    }//end if

	    else if(diffarray[row] >= -0.20 && diffarray[row] <= 0.00)
	    {
	      zundel1++;
	    }//end if

	    else if(diffarray[row] >= 0.00 && diffarray[row] <= 0.20)
	    {
	      zundel2++;
	    }//end if

	    else if(diffarray[row] >= 0.20 && diffarray[row] <= 0.40)
	    {
	      count3++;
	    }//end if
	
	    else if(diffarray[row] >= 0.40 && diffarray[row] <= 0.60)
	    {
	      count4++;
	    }//end if
 
            else if(diffarray[row] >= 0.60 && diffarray[row] <= 2.0)
	    {
	      countmax++;
	    }//end if


	  }//end row-loop
	
	  fprintf(output,"%d %d %d %d %d %d %d %d\n",
	  countmin,
	  count1,
	  count2,
          zundel1,
	  zundel2,
	  count3,
	  count4,
	  countmax);

/*	  for(x=0; x<n; x++)
          {
             printf("%f\n",diffarray[x]);
          }//end for
*/	  
	  fclose(output);

	}//end main for

}
