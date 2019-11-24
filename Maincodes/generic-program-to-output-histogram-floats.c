#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
   int indexnum,num,dlines,x,y,count,clines;
   int crow,drow,sum,numbins,numlines,hold;
   int snap,Oxy,Hyd,z,coln;
   float weight,dist,threshold;
   FILE *input,*output;

   numbins=23; //# of desired bins
   numlines=143000; //# of lines in file
   threshold = 0.0000001;

   float countmatrix[numbins][3][26]; //countmatrix[# of predefined bins+1][# colns]
   float distwgt[numlines][2]; //distwgt[# of lines in inputfile][# colns]

//Initialize matrix
   for(x=0; x<numbins; x++){
      for(y=0; y<3; y++){
         for(z=0; z<26; z++){
	   countmatrix[x][y][z]=0.000;
	 }
      }
   }
   for(x=0; x<numlines; x++){
      distwgt[x][0]=0;
      distwgt[x][1]=0;
   }
   clines=1;
   dlines=1;
   sum=0;

//   output=fopen("hotmap-cent-rep-Zundel-individ-weight-dist.txt","w");
   input=fopen("reorganized-data.txt","r");
//Create matrix with proper x and y markers
    for(indexnum=1; indexnum<numbins; indexnum++){
       if(indexnum==1){


          countmatrix[clines][0][0]=1.081; //pick initial value;
          countmatrix[clines][1][0]=countmatrix[clines][0][0]+0.009; //should be (bin value subtracted by 1, ex. 20-1=19)

          countmatrix[clines][0][1]=0.250;
          countmatrix[clines][0][2]=0.281;
          countmatrix[clines][0][3]=0.312;
          countmatrix[clines][0][4]=0.344;
          countmatrix[clines][0][5]=0.375;
          countmatrix[clines][0][6]=0.406;
          countmatrix[clines][0][7]=0.438;
          countmatrix[clines][0][8]=0.469;
          countmatrix[clines][0][9]=0.500;
          countmatrix[clines][0][10]=0.531;
          countmatrix[clines][0][11]=0.562;
          countmatrix[clines][0][12]=0.594;
          countmatrix[clines][0][13]=0.625;
          countmatrix[clines][0][14]=0.656;
          countmatrix[clines][0][15]=0.688;
          countmatrix[clines][0][16]=0.719;
          countmatrix[clines][0][17]=0.750;
          countmatrix[clines][0][18]=0.781;
          countmatrix[clines][0][19]=0.812;
          countmatrix[clines][0][20]=0.844;
          countmatrix[clines][0][21]=0.875;
          countmatrix[clines][0][22]=0.906;
          countmatrix[clines][0][23]=0.938;
          countmatrix[clines][0][24]=0.969;
          countmatrix[clines][0][25]=1.000;
          clines++;
       }
       else if(indexnum==2){
          countmatrix[clines][0][0]=1.091; //pick initial value;
          countmatrix[clines][1][0]=1.100; //should be (bin value subtracted by 1, ex. 20-1=19)

          countmatrix[clines][0][1]=0.250;
          countmatrix[clines][0][2]=0.281;
          countmatrix[clines][0][3]=0.312;
          countmatrix[clines][0][4]=0.344;
          countmatrix[clines][0][5]=0.375;
          countmatrix[clines][0][6]=0.406;
          countmatrix[clines][0][7]=0.438;
          countmatrix[clines][0][8]=0.469;
          countmatrix[clines][0][9]=0.500;
          countmatrix[clines][0][10]=0.531;
          countmatrix[clines][0][11]=0.562;
          countmatrix[clines][0][12]=0.594;
          countmatrix[clines][0][13]=0.625;
          countmatrix[clines][0][14]=0.656;
          countmatrix[clines][0][15]=0.688;
          countmatrix[clines][0][16]=0.719;
          countmatrix[clines][0][17]=0.750;
          countmatrix[clines][0][18]=0.781;
          countmatrix[clines][0][19]=0.812;
          countmatrix[clines][0][20]=0.844;
          countmatrix[clines][0][21]=0.875;
          countmatrix[clines][0][22]=0.906;
          countmatrix[clines][0][23]=0.938;
          countmatrix[clines][0][24]=0.969;
          countmatrix[clines][0][25]=1.000;
          clines++;
	} 
	else if(indexnum>2){
	  countmatrix[clines][0][0]=countmatrix[clines-1][0][0]+0.010; //change accordingly
	  countmatrix[clines][1][0]=countmatrix[clines-1][1][0]+0.010; //change

          countmatrix[clines][0][1]=0.250;
          countmatrix[clines][0][2]=0.281;
          countmatrix[clines][0][3]=0.312;
          countmatrix[clines][0][4]=0.344;
          countmatrix[clines][0][5]=0.375;
          countmatrix[clines][0][6]=0.406;
          countmatrix[clines][0][7]=0.438;
          countmatrix[clines][0][8]=0.469;
          countmatrix[clines][0][9]=0.500;
          countmatrix[clines][0][10]=0.531;
          countmatrix[clines][0][11]=0.562;
          countmatrix[clines][0][12]=0.594;
          countmatrix[clines][0][13]=0.625;
          countmatrix[clines][0][14]=0.656;
          countmatrix[clines][0][15]=0.688;
          countmatrix[clines][0][16]=0.719;
          countmatrix[clines][0][17]=0.750;
          countmatrix[clines][0][18]=0.781;
          countmatrix[clines][0][19]=0.812;
          countmatrix[clines][0][20]=0.844;
          countmatrix[clines][0][21]=0.875;
          countmatrix[clines][0][22]=0.906;
          countmatrix[clines][0][23]=0.938;
          countmatrix[clines][0][24]=0.969;
          countmatrix[clines][0][25]=1.000;
	  clines++;

	}
      }//end index-loop
   
//      for(x=1; x<clines; x++){ 
//printf("%0.3f %0.3f\n",countmatrix[x][0][0],countmatrix[x][1][0]);
//        for(y=1; y<26; y++){
//            printf("%0.3f\n",countmatrix[0][0][y]);
//        }
//      }

//Read in lifetime vs. count input file
      while(fscanf(input,"%d %d %d %f %f\n",&snap,&Oxy,&Hyd,&dist,&weight)==5){
	 distwgt[dlines][0]=dist;
	 distwgt[dlines][1]=weight;
         dlines++;
      }//end while-loop
      fclose(input);

      for(drow=1; drow<dlines; drow++){
		 if(distwgt[drow][0]<100.0){
            for(crow=1; crow<clines; crow++){
               if(distwgt[drow][0] > countmatrix[crow][0][0] && distwgt[drow][0] < countmatrix[crow][1][0]){
//if distance is within range
	              for(coln=1; coln<26; coln++){
	                 if(distwgt[drow][1]-threshold < countmatrix[crow][0][coln] && distwgt[drow][1]+threshold > countmatrix[crow][0][coln]){ //if weight match
			           hold=coln;
			           countmatrix[crow][1][hold]=countmatrix[crow][1][hold]+1.000;
			           distwgt[drow][0]=100.0;
			           distwgt[drow][1]=100.0;
      		        }
                 }//end coln-loop
              }
			  else if(distwgt[drow][0]-threshold <= countmatrix[crow][0][0] && distwgt[drow][0]+threshold >= countmatrix[crow][0][0] ){
	              for(coln=1; coln<26; coln++){
	                 if(distwgt[drow][1]-threshold < countmatrix[crow][0][coln] && distwgt[drow][1]+threshold > countmatrix[crow][0][coln]){
				  //if distance is equal to min. range
	                   hold=coln;
	                   countmatrix[crow][1][hold]=countmatrix[crow][1][hold]+1.000;
	                   distwgt[drow][0]=100.0;
	                   distwgt[drow][1]=100.0;
				   }
			    }
			  }
			  else if(distwgt[drow][0]-threshold <= countmatrix[crow][1][0] && distwgt[drow][0]+threshold >= countmatrix[crow][1][0] ){
	              for(coln=1; coln<26; coln++){
	                 if(distwgt[drow][1]-threshold < countmatrix[crow][0][coln] && distwgt[drow][1]+threshold > countmatrix[crow][0][coln]){
				  //if distance is equal to max. range
	                    hold=coln;
	                    countmatrix[crow][1][hold]=countmatrix[crow][1][hold]+1.000;
	                    distwgt[drow][0]=100.0;
	                    distwgt[drow][1]=100.0;
					}
				 }	
			  }
           }		    
	    }//end coln-loop
     }//end drow-loop
 	 
//        for(x=1; x<dlines; x++){
//          if(distwgt[x][0]<100.0){ 
//            printf("%0.3f %0.3f\n",distwgt[x][0],distwgt[x][1]);
//	       }
//         }


 //     for(x=1; x<clines; x++){ 
 //        printf("%d %0.3f %0.3f\n",x,countmatrix[x][0][0],countmatrix[x][1][0]);
 //     }

    for(x=1; x<clines; x++){
	    for(y=1; y<26; y++){
		   printf("%0.3f %0.3f %0.3f %0.1f\n",countmatrix[x][0][0],countmatrix[x][1][0],countmatrix[x][0][y],countmatrix[x][1][y]);
       }
    }

//      fclose(output);

}//end main-loop
