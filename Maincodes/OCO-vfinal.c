#include<stdio.h>
#include<stdlib.h>

//******************************************************
//!!!The previous program was reading in the wrong
//!!!variables - that was modifed in this version.
//
//Program created to obtain the edge distribution about
//an H atom and surrounding over-coordinated oxygens.
//
//Input: OCO-list.txt files from identify-EZ.c program  
//	 Covalent bond (0.0-1.3A) GraphGeod files
//    
//Output: Main file containing the counts of OCO to an H
//
//!!GraphGeod files need to be organized by the second 
//column (H indices) prior to analysis!!
//
//Created by: Lelee Ounkham
//Modified: August 8, 2017
//
//******************************************************

int main()
{

//Parameters to read inputs
     int Oindex,Hindex,Hhold1,Hhold2,O,H;
     int holdOxy,int1,int2,int3,int4,int5;
     int snap,OCO,maxsnap;
     int OCOrow,Hrow,countOCO;
     float float1,float2;
     char ifile[100],ifile2[100],ofile[100];
     FILE *input,*input2,*output;  
  
//Parameters for the arrays
     int mxline=300;
     int Harray[mxline][2];
     int OCOarray[mxline][2]; 
     int holdarray[3][2];

//Parameters for loops and counters
     int i,j,k,x,y,z;
     int lines,lines2;
     int OCOwater,zundel,unknown;

//Intialize counters

    maxsnap=55001; //Change accordingly!

//Identifying OCO's 

	 OCOwater=0;
	 zundel=0;
	 unknown=0;

        for(i=1; i<maxsnap+1; i++){
           sprintf(ifile,"R-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",i,i);
           input=fopen(ifile,"r");
   
	   j=1;
	   lines=0;
	   lines2=0;

          while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oindex,&Hindex,&int1,&int2,&int3,&int4,&int5,&float1,&float2)==9){
              if(j==1){
                 holdOxy=Oindex;
                 Hhold1=Hindex;
              }
              else if(j==2){
                 if(Oindex==holdOxy){ //if O is the same on next line
                    Hhold2=Hindex;
                 }
                 else{
                    j=1;
                    holdOxy=Oindex;
                    Hhold1=Hindex;
                 }
              }
              else if(j==3){ //if there is a 3rd occurrence then print into output
                 if(Oindex==holdOxy){
		    OCOarray[lines][0]=holdOxy;
		    OCOarray[lines][1]=Hhold1;
		    lines++;

		    OCOarray[lines][0]=holdOxy;
		    OCOarray[lines][1]=Hhold2;
		    lines++;

		    OCOarray[lines][0]=holdOxy;
		    OCOarray[lines][1]=Hindex;
		    lines++;

                    j=0; //restarts Oindex to read the next line
                 }
                 else{
                    j=1;
                    holdOxy=Oindex;
                    Hhold1=Hindex;
                 }
              }
              j++;
           } //end while

           fclose(input);

/*	   for(z=0; z<lines; z++){
	      printf("%d %d\n",OCOarray[z][0],OCOarray[z][1]);
	   }
*/
//Read in OCO-list from previous section and GraphGeod files
	   sprintf(ifile2,"S-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",i,i);
           input2=fopen(ifile2,"r");

           while(fscanf(input2,"%d %d %d %d %d %d %d %f %f\n",&O,&H,&int1,&int2,&int3,&int4,&int5,&float1,&float2)==9){
	      Harray[lines2][0]=O;
	      Harray[lines2][1]=H;
	      lines2++;

           }//end while-loop
   
	   fclose(input2);

// 	   for(z=0; z<lines2; z++)
//	   {
//	      printf("%d %d\n",Harray[z][0],Harray[z][1]);
//	   }

//******************************************************
//Section to delete any duplicates in Harray
//******************************************************

	   for(x=0; x<lines2; x++){
	      for(y=1; y<lines2; y++){
	         if(Harray[x][1]!=0){
		    if(Harray[x][1]==Harray[y][1]){ //if first H matches the next - duplicate
		       if(Harray[x][0]!=Harray[y][0]){ //If O's don't match
	    	          Harray[y][1]=0;
		       }
		    }
	         }
	      }//end y-loop
	   }//end x-loop

/*	   for(z=0; z<lines2; z++){
	      if(Harray[z][1]!=0){
		 printf("%d %d\n",Harray[z][0],Harray[z][1]);
	      }
	   }
*/
//******************************************************
//Section to count the number of OCO per H in a 
//GraphGeod file by comparing Harray to OCOarray
//******************************************************

//Intialize counters and arrays	
	   countOCO=0;

	   for(z=0; z<4; z++){
	      holdarray[z][0]=0;
	      holdarray[z][1]=0;
	   }
  
	   for(Hrow=0; Hrow<lines2; Hrow++){
	      for(OCOrow=0; OCOrow<lines; OCOrow++){ 
		 if(Harray[Hrow][1]!=0 && OCOarray[OCOrow][1]!=0){
	            if(OCOarray[OCOrow][1]==Harray[Hrow][1]){ //If H's matches
		       holdarray[0][0]=OCOarray[OCOrow][0]; //hold OCO
		       holdarray[0][1]=Harray[Hrow][1];  //hold H index
		       countOCO++;
//Loop over OCOarray to ensure there aren't more OCO's bonded to that H index
		       for(z=0; z<lines; z++){ //go through OCOarray again
			  if(OCOarray[OCOrow][0]!=OCOarray[z][0]){//If O's are different
			     if(holdarray[0][1]==OCOarray[z][1]){ //and H's are the same 
				if(holdarray[1][0]==0){ //for H with edge of 2 - Zundel identified
				   holdarray[1][0]=OCOarray[z][0]; //hold info 
				   holdarray[1][1]=OCOarray[z][1];
			           countOCO++;
				   OCOarray[z][0]=0;
				   OCOarray[z][1]=0; //Zero out array
//printf("%d %d %d\n",OCOarray[z][0],OCOarray[z][1],countOCO);
				}
				else if(holdarray[1][0]!=0){ //for H with edge of 3
				   holdarray[2][0]=OCOarray[z][0];
				   holdarray[2][1]=OCOarray[z][1];
				   countOCO++;
			        } 
			     }
		          }
		       }//end z-loop
		       Harray[Hrow][1]=0; //remove H from list
		    }
		 }
		 if(holdarray[0][1]!=0){
		    if(countOCO==1){
		       OCOwater++;
		    }
		    else if(countOCO==2){
		       zundel++;
		    }
		    else if(countOCO==3){
		       unknown++;
		    }
		 }
		 countOCO=0;

		    for(z=0; z<3; z++){
		       holdarray[z][0]=0;
		       holdarray[z][1]=0;
		    }
	      }//end OCOrow-loop
	   }//end Hrow-loop

	   for(z=0; z<lines; z++){
	      OCOarray[z][0]=0;
	      OCOarray[z][1]=0;
	   }

	   for(z=0; z<lines2; z++){
	      Harray[z][0]=0;
	      Harray[z][1]=0;
	   }

	}//end i-loop
	
     //   output=fopen("OCO-edge-distribution.txt","w");
      
	printf("1 edge: %d\n2 edge: %d\n3 edge: %d\n", OCOwater,zundel,unknown);
       // fclose(output);

} //end of main-loop
