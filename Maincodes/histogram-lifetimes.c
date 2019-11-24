#include<stdio.h>
#include<stdlib.h>

int main(){
   int indexnum,num,plines,x,y,count,persist,clines;
   int crow,prow,sum;
   FILE *input;

   int countmatrix[1442][3];
   int persistmatrix[750][2];

//Initialize matrix
   for(x=0; x<1442; x++){
      for(y=0; y<3; y++){
	 countmatrix[x][y]=0;
      }
   }
   for(x=0; x<750; x++){
      persistmatrix[x][0]=0;
      persistmatrix[x][1]=0;
   }

   clines=0;
   plines=1;
   sum=0;

   input=fopen("lifetimes.txt","r");
//Create matrix with proper x and y markers
    for(indexnum=1; indexnum<=1442; indexnum++){
       if(indexnum==1){
          countmatrix[clines][0]=indexnum;
          countmatrix[clines][1]=indexnum+19;//should be 20
//          printf("%d %d\n",indexnum,indexnum+19);
	  clines++;
	} 
	else if(indexnum!=1){
	  countmatrix[clines][0]=countmatrix[clines-1][0]+20;
	  countmatrix[clines][1]=countmatrix[clines][0]+19;
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
	         sum=persistmatrix[prow][1]; 
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
         printf("%d %d %d\n",countmatrix[x][0],countmatrix[x][1],countmatrix[x][2]);
//	   printf("%d %d\n",persistmatrix[x][0],persistmatrix[x][1]);
      }
 
}//end main-loop
