#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
	
	int x,row1,row2,wlines,freq1,freq2,num;
	float O1Hweight,O2Hweight,threshold;
	FILE *input,*output;
	
	int countmatrix[500][2];
	float weightmatrix[500][2];
	
//Initialize matrices
	for(x=0; x<500; x++){
		countmatrix[x][0]=0;
		countmatrix[x][1]=0;
		weightmatrix[x][0]=0.0;
		weightmatrix[x][1]=0.0;
	}	
	wlines=0;
	threshold=0.0000001;
	
//Allocate data into matrix
   output=fopen("z-0.3-total-counts.txt","a");	
   input=fopen("z-0.3.txt","r");
   while(fscanf(input,"%f %f %d\n",&O1Hweight,&O2Hweight,&num)==3){
	   weightmatrix[wlines][0]=O1Hweight;
	   weightmatrix[wlines][1]=O2Hweight;
	   countmatrix[wlines][0]=num;
	   wlines++;
//printf("%0.3f %0.3f %d\n",O1Hweight,O2Hweight,freq);
   }
   fclose(input);
   
//   for(x=0; x<wlines; x++){
//     printf("%0.3f %0.3f %d\n",weightmatrix[x][0],weightmatrix[x][1],countmatrix[x][0]);
//   } 

   
//Identify structures with identical weights, but in different positions
   for(row1=0; row1<wlines; row1++){
	  for(row2=1; row2<wlines; row2++){
		 if(countmatrix[row2][0]!=2000 && countmatrix[row1][0]!=2000 && countmatrix[row1][0]!=countmatrix[row2][0]){ //if element is nonzero
			if(weightmatrix[row1][0]-weightmatrix[row2][1] <= threshold && weightmatrix[row1][0]+threshold >= weightmatrix[row2][1]
			&& weightmatrix[row1][1]-weightmatrix[row2][0] <= threshold && weightmatrix[row1][1]+threshold >= weightmatrix[row2][0]){
//If boths weights are different position, then difference and sum will be less than 0.0000001
				freq2=countmatrix[row2][0]; //store count
				freq1=countmatrix[row1][0]; //store count
//printf("trigger: %d %d %0.3f %0.3f %0.3f %0.3f\n",freq2,countmatrix[row1][0],weightmatrix[row1][0],weightmatrix[row2][0],weightmatrix[row1][1],weightmatrix[row2][1]);
				  countmatrix[row1][1]=freq1+freq2; //add to overall sum
//printf("trigger2: %d %0.3f %0.3f %0.3f %0.3f\n",countmatrix[row1][1],weightmatrix[row1][0],weightmatrix[row1][1],weightmatrix[row1][0],weightmatrix[row1][1]);

				countmatrix[row2][0]=2000;
				countmatrix[row1][0]=2000;
		    }
   		    else if(weightmatrix[row1][0]-weightmatrix[row2][0] <= threshold && weightmatrix[row1][0]+threshold >= weightmatrix[row2][0]
   			&& weightmatrix[row1][1]-weightmatrix[row2][1] <= threshold && weightmatrix[row1][1]+threshold >= weightmatrix[row2][1]){
 //If boths weights are different position, then difference and sum will be less than 0.0000001
   				freq2=countmatrix[row2][0]; //store count
   				freq1=countmatrix[row1][0]; //store count
//printf("trigger2: %d %d %0.3f %0.3f %0.3f %0.3f\n",freq2,countmatrix[row1][0],weightmatrix[row1][0],weightmatrix[row2][0],weightmatrix[row1][1],weightmatrix[row2][1]);
   				countmatrix[row1][1]=freq1+freq2; //add to overall sum
   //printf("trigger2: %d %0.3f %0.3f %0.3f %0.3f\n",countmatrix[row1][1],weightmatrix[row1][0],weightmatrix[row1][1],weightmatrix[row1][0],weightmatrix[row1][1]);

   				countmatrix[row2][0]=2000;
   				countmatrix[row1][0]=2000; 
		    } 
	     }
      }//end row2-loop
   }//end row1-loop
 
   for(x=0; x<wlines; x++){
      if(countmatrix[x][0]!=2000){
	 fprintf(output,"%0.3f %0.3f %d\n",weightmatrix[x][0],weightmatrix[x][1],countmatrix[x][0]);
      }
      else if(countmatrix[x][1]!=0){
         fprintf(output,"%0.3f %0.3f %d\n",
	 weightmatrix[x][0],
	 weightmatrix[x][1],
	 countmatrix[x][1]);	   
      }
   }//end x-loop
   fclose(output);
}//end main-loop
