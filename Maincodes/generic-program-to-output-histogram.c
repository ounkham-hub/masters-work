#include<stdio.h>
#include<stdlib.h>

int main(){
   int indexnum,num,plines,x,y,count,persist,clines;
   int crow,prow,sum,numbins,numlines;
   FILE *input,*output;

   numbins=28821; //# of desired bins
   numlines=755; //# of lines in file

   int countmatrix[numbins][3]; //countmatrix[# of predefined bins+1][# colns]
   int persistmatrix[numlines][2]; //persistmatrix[# of lines in inputfile][# colns]

//Initialize matrix
   for(x=0; x<numbins; x++){
      for(y=0; y<3; y++){
	 countmatrix[x][y]=0;
      }
   }
   for(x=0; x<numlines; x++){
      persistmatrix[x][0]=0;
      persistmatrix[x][1]=0;
   }

   clines=0;
   plines=1;
   sum=0;

   output=fopen("lifetimes-histogram.txt","w");
   input=fopen("lifetimes.txt","r");
//Create matrix with proper x and y markers
    for(indexnum=1; indexnum<=numbins; indexnum++){
       if(indexnum==1){
          countmatrix[clines][0]=indexnum;
          countmatrix[clines][1]=indexnum+9;//should be (bin value subtracted by 1, ex. 20-1=19)
//          printf("%d %d\n",indexnum,indexnum+19);
	  clines++;
	} 
	else if(indexnum!=1){
	  countmatrix[clines][0]=countmatrix[clines-1][0]+10; //change accordingly
	  countmatrix[clines][1]=countmatrix[clines][0]+9; //change
	  clines++;
	}
      }//end index-loop
   
//      for(x=1; x<clines; x++){ 
//         printf("%d %d %d\n",x,countmatrix[x][0],countmatrix[x][1]);
//      }
//Read in lifetime vs. count input file
      while(fscanf(input,"%d %d\n",&persist,&count)==2){
	 persistmatrix[plines][0]=persist;
	 persistmatrix[plines][1]=count;
         plines++;
      }//end while-loop
      fclose(input);

//      for(x=1; x<plines; x++){ 
//         printf("%d %d %d\n",x,persistmatrix[x][0],persistmatrix[x][1]);
//      }


      for(prow=1; prow<plines; prow++){
         for(crow=0; crow<clines; crow++){
	    if(persistmatrix[prow][0]!=0){
 	       if(persistmatrix[prow][0]>=countmatrix[crow][0] && persistmatrix[prow][0]<=countmatrix[crow][1]){
//If persistence value is between index range
	         sum=persistmatrix[prow][0];  //change back to 1!
//printf("%d %d %d\n",countmatrix[crow][0],countmatrix[crow][1],persistmatrix[prow][1]);
	         countmatrix[crow][2]=countmatrix[crow][2]+sum; //update sum

//	         persistmatrix[prow][0]=0;
//	         persistmatrix[prow][1]=0;
//	         persistmatrix[prow][2]=0;
	      } 
	   }//end crow-loop
	}
      }//end row-loop  
    
      for(x=0; x<clines; x++){
         fprintf(output,"%d %d %d\n",countmatrix[x][0],countmatrix[x][1],countmatrix[x][2]);
//	   printf("%d %d\n",persistmatrix[x][0],persistmatrix[x][1]);
      }
      fclose(output);
}//end main-loop
