#include<stdio.h>
#include<stdlib.h>

//*****************************************************************
//Program written to sort updated GraphGeod files in a single 
//snapshot and bin them according to a range of values and 
//count the frequency per bin.
//
//To obtain deltaOH, bin-diff-updated-GG.sh is required!
//
//8 outputs are created:
//map1: O-O distance from 2.3-2.35A,deltaOH from -2 to +2
//map2: O-O distance from 2.35-2.4A,deltaOH from -2 to +2
//map3: O-O distance from 2.4-2.45A,deltaOH from -2 to +2
//map4: O-O distance from 2.45-2.5A,deltaOH from -2 to +2	
//map5: O-O distance from 2.5-2.55A,deltaOH from -2 to +2
//map6: O-O distance from 2.55-2.6A,deltaOH from -2 to +2
//map7: O-O distance from 2.6-2.65A,deltaOH from -2 to +2
//map8: O-O distance from 2.65-2.7A,deltaOH from -2 to +2	
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
	char ofile5[100],ofile6[100],ofile7[100],ofile8[100];
	FILE *input,*output1,*output2,*output3,*output4;
	FILE *output5,*output6,*output7,*output8;

//Parameters for diffarray
	int l=300;
	int i,n,x;
	float diffarray[l][2];
//Counter parameters
	int row,coln;
    int count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12,count13;

	for(i=1; i<=61555; i++)
   	{
           sprintf(ifile,"OH-diff%d.GG",i);
	   sprintf(ofile1,"map1.txt",i);
 	   sprintf(ofile2,"map2.txt",i);
	   sprintf(ofile3,"map3.txt",i);
           sprintf(ofile4,"map4.txt",i);
	   sprintf(ofile5,"map5.txt",i);
	   sprintf(ofile6,"map6.txt",i);
	   sprintf(ofile7,"map7.txt",i);
           sprintf(ofile8,"map8.txt",i);
	   input=fopen(ifile,"r");
           output1=fopen(ofile1,"a");
	   output2=fopen(ofile2,"a");
	   output3=fopen(ofile3,"a");
	   output4=fopen(ofile4,"a");
	   output5=fopen(ofile5,"a");
	   output6=fopen(ofile6,"a");
	   output7=fopen(ofile7,"a");
	   output8=fopen(ofile8,"a");
	   
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
	        if(diffarray[row][0] >= 2.3 && diffarray[row][0] <= 2.35)
	        {
	          if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50)
	          {
	            count1++;
	          }//end if
	          else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40)
	          {
	            count2++;
	          }//end if
	          else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30)
	          {
	            count3++;
	          }//end if
	          else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20)
	          {
	            count4++;
	          }//end if
              else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10)
              {
                  count5++;
              }//end if
              else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00)
              {
                  count6++;
              }//end if
	          else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10)
	          {
	            count7++;
	          }//end if
	          else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20)
	          {
	            count8++;
	          }//end if
	          else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30)
	          {
	            count9++;
	          }//end if
              else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40)
	          {    
	            count10++;
	          }//end if
              else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50)
              {
                count11++;
              }//end if
              else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60)
              {
                count12++;
              }//end if
	        }//end 2.3-2.35-for-loop
	      }//end for-loop
	        fprintf(output1,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
	        count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12);
	
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
        
	     for(row=0; row<n; row++)
	     {
	      if(diffarray[row][0] >= 2.35 && diffarray[row][0] <= 2.4)
	      {
              if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50)
              {
                  count1++;
              }//end if
              else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40)
              {
                  count2++;
              }//end if
              else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30)
              {
                  count3++;
              }//end if
              else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20)
              {
                  count4++;
              }//end if
              else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10)
              {
                  count5++;
              }//end if
              else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00)
              {
                  count6++;
              }//end if
              else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10)
              {
                  count7++;
              }//end if
              else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20)
              {
                  count8++;
              }//end if
              else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30)
              {
                  count9++;
              }//end if
              else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40)
              {
                  count10++;
              }//end if
              else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50)
              {
                  count11++;
              }//end if
              else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60)
              {
                  count12++;
              }//end if
	       }//end 2.35-2.4-for-loop
	      }//end for-loop	   
        fprintf(output2,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
        count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12);

	fclose(output2);	
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
       
        
	     for(row=0; row<n; row++)
	     {
	      if(diffarray[row][0] >= 2.4 && diffarray[row][0] <= 2.45)
	      {
              if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50)
              {
                  count1++;
              }//end if
              else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40)
              {
                  count2++;
              }//end if
              else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30)
              {
                  count3++;
              }//end if
              else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20)
              {
                  count4++;
              }//end if
              else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10)
              {
                  count5++;
              }//end if
              else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00)
              {
                  count6++;
              }//end if
              else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10)
              {
                  count7++;
              }//end if
              else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20)
              {
                  count8++;
              }//end if
              else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30)
              {
                  count9++;
              }//end if
              else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40)
              {
                  count10++;
              }//end if
              else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50)
              {
                  count11++;
              }//end if
              else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60)
              {
                  count12++;
              }//end if
	       }//end 2.4-2.45-for-loop
	      }//end for-loop
        
        fprintf(output3,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
        count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12);

	   fclose(output3);
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
       
        
	      for(row=0; row<n; row++)
 	      {
	       if(diffarray[row][0] >= 2.45 && diffarray[row][0] <= 2.5)
	       {
               if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50)
               {
                   count1++;
               }//end if
               else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40)
               {
                   count2++;
               }//end if
               else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30)
               {
                   count3++;
               }//end if
               else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20)
               {
                   count4++;
               }//end if
               else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10)
               {
                   count5++;
               }//end if
               else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00)
               {
                   count6++;
               }//end if
               else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10)
               {
                   count7++;
               }//end if
               else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20)
               {
                   count8++;
               }//end if
               else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30)
               {
                   count9++;
               }//end if
               else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40)
               {
                   count10++;
               }//end if
               else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50)
               {
                   count11++;
               }//end if
               else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60)
               {
                   count12++;
               }//end if
	        }//end 2.45-2.5-for-loop
	       }//end for
        
        fprintf(output4,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
        count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12);
        
        fclose(output4);

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
        
	    for(row=0; row<n; row++)
 	      {
	       if(diffarray[row][0] >= 2.5 && diffarray[row][0] <= 2.55)
	       {
               if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50)
               {
                   count1++;
               }//end if
               else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40)
               {
                   count2++;
               }//end if
               else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30)
               {
                   count3++;
               }//end if
               else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20)
               {
                   count4++;
               }//end if
               else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10)
               {
                   count5++;
               }//end if
               else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00)
               {
                   count6++;
               }//end if
               else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10)
               {
                   count7++;
               }//end if
               else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20)
               {
                   count8++;
               }//end if
               else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30)
               {
                   count9++;
               }//end if
               else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40)
               {
                   count10++;
               }//end if
               else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50)
               {
                   count11++;
               }//end if
               else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60)
               {
                   count12++;
               }//end if
	        }//end 2.5-2.55-for-loop
	       }//end for
        
        fprintf(output5,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
        count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12);
        
	       fclose(output5);

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
        
	     for(row=0; row<n; row++)
 	      {
	       if(diffarray[row][0] >= 2.55 && diffarray[row][0] <= 2.60)
	       {
               if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50)
               {
                   count1++;
               }//end if
               else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40)
               {
                   count2++;
               }//end if
               else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30)
               {
                   count3++;
               }//end if
               else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20)
               {
                   count4++;
               }//end if
               else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10)
               {
                   count5++;
               }//end if
               else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00)
               {
                   count6++;
               }//end if
               else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10)
               {
                   count7++;
               }//end if
               else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20)
               {
                   count8++;
               }//end if
               else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30)
               {
                   count9++;
               }//end if
               else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40)
               {
                   count10++;
               }//end if
               else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50)
               {
                   count11++;
               }//end if
               else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60)
               {
                   count12++;
               }//end if
	        }//end 2.55-2.6-for-loop
	       }//end for	   

        fprintf(output6,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
        count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12);
        
	       fclose(output6);

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
        
	     for(row=0; row<n; row++)
 	      {
	       if(diffarray[row][0] >= 2.60 && diffarray[row][0] <= 2.65)
	       {
               if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50)
               {
                   count1++;
               }//end if
               else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40)
               {
                   count2++;
               }//end if
               else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30)
               {
                   count3++;
               }//end if
               else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20)
               {
                   count4++;
               }//end if
               else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10)
               {
                   count5++;
               }//end if
               else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00)
               {
                   count6++;
               }//end if
               else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10)
               {
                   count7++;
               }//end if
               else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20)
               {
                   count8++;
               }//end if
               else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30)
               {
                   count9++;
               }//end if
               else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40)
               {
                   count10++;
               }//end if
               else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50)
               {
                   count11++;
               }//end if
               else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60)
               {
                   count12++;
               }//end if
	        }//end 2.6-2.65-for-loop
	       }//end for
        
        fprintf(output7,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
        count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12);
        

	       fclose(output7);

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
        
	      for(row=0; row<n; row++)
 	      {
	        if(diffarray[row][0] >= 2.65 && diffarray[row][0] <= 2.70)
	       {
               if(diffarray[row][1] > -0.60 && diffarray[row][1] <= -0.50)
               {
                   count1++;
               }//end if
               else if(diffarray[row][1] > -0.50 && diffarray[row][1] <= -0.40)
               {
                   count2++;
               }//end if
               else if(diffarray[row][1] > -0.40 && diffarray[row][1] <= -0.30)
               {
                   count3++;
               }//end if
               else if(diffarray[row][1] > -0.30 && diffarray[row][1] <= -0.20)
               {
                   count4++;
               }//end if
               else if(diffarray[row][1] > -0.20 && diffarray[row][1] <= -0.10)
               {
                   count5++;
               }//end if
               else if(diffarray[row][1] > -0.10 && diffarray[row][1] < 0.00)
               {
                   count6++;
               }//end if
               else if(diffarray[row][1] > 0.00 && diffarray[row][1] <= 0.10)
               {
                   count7++;
               }//end if
               else if(diffarray[row][1] > 0.10 && diffarray[row][1] <= 0.20)
               {
                   count8++;
               }//end if
               else if(diffarray[row][1] > 0.20 && diffarray[row][1] <= 0.30)
               {
                   count9++;
               }//end if
               else if(diffarray[row][1] > 0.30 && diffarray[row][1] <= 0.40)
               {
                   count10++;
               }//end if
               else if(diffarray[row][1] > 0.40 && diffarray[row][1] <= 0.50)
               {
                   count11++;
               }//end if
               else if(diffarray[row][1] > 0.50 && diffarray[row][1] <= 0.60)
               {
                   count12++;
               }//end if
	        }//end 2.65-2.7-for-loop
	       }//end for	   

        fprintf(output8,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
        count1,count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12);
        

	       fclose(output8);
	
	}//end i-loop

}
