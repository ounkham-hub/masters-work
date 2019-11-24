#include<stdio.h>
#include<stdlib.h>

//*********************************************************************************************
//Program created to identify true PT based on Zundel switch events and from a predetermiend
//list of original and new O index list.
//
//Input: uniq-exchange-OH-indices.txt
//	 #-Zundel-switches.txt
//
//Created by Lelee Ounkham
//Last modified February 25, 2019
//
//*********************************************************************************************

int main(){

   int Hyd,O1,O2,startsnap,finalsnap,Hmax,plines;
   int sharedH,origOindex,nextOindex;
   int tempH,tempO1,tempO2,diff;
   int proton,slines,y,z,cutoff,srow,prow;
   FILE *input,*pinput,*output;
   char ifile[100];

   int switchlist[10000][5];
   int pairlist[10000][4];


//Initialize matrix and counters
   for(y=0; y<10000; y++){
      for(z=0; z<6; z++){
	 switchlist[y][z]=0;
      }
      for(z=0; z<4; z++){
	pairlist[y][z]=0;
      }
   }
   plines=0;
   slines=0;
   Hmax=314;
   cutoff=0;
   
//*********************************************************************************************
//Section I: Store Zundel switch information and unique O pair indices into appropriate 
//matrix.
//*********************************************************************************************
   output=fopen("complete-list-of-true-PTs-NOpersist.txt","a");
   pinput=fopen("uniq-exchange-OH-indices.txt","r");
   while(fscanf(pinput,"%d %d %d\n",&sharedH,&origOindex,&nextOindex)==3){
      pairlist[plines][0]=sharedH;
      pairlist[plines][1]=origOindex;
      pairlist[plines][2]=nextOindex;
      plines++;
   }
   fclose(pinput);
/*
  for(y=0; y<plines; y++){
      printf("%d %d %d\n",
      pairlist[y][0],
      pairlist[y][1],
      pairlist[y][2]);
   }
*/
   for(proton=103; proton<=Hmax; proton++){
      sprintf(ifile,"%d-Zundel-switches.txt",proton);
      input=fopen(ifile,"r");   

      while(fscanf(input,"%d %d %d %d %d\n",&Hyd,&O1,&O2,&startsnap,&finalsnap)==5){
         switchlist[slines][0]=Hyd;
         switchlist[slines][1]=O1;
         switchlist[slines][2]=O2;
         switchlist[slines][3]=startsnap;
         switchlist[slines][4]=finalsnap;
         slines++;
      }
      fclose(input);
/*
      for(y=0; y<slines; y++){
         printf("%d %d %d %d %d\n",
         switchlist[y][0],
         switchlist[y][1],
         switchlist[y][2],
	 switchlist[y][3],
	 switchlist[y][4]);
      }
*/

/*
      for(prow=0; prow<plines; prow++){
         if(prow==0){ //store first pair 
	    tempH=pairlist[prow][0];
	    tempO1=pairlist[prow][1];
            tempO2=pairlist[prow][2];
	  
//	    printf("%d %d %d\n",tempH,tempO1,tempO2);
	  }
       }

*/

      for(prow=0; prow<plines; prow++){
         if(prow==0){ //store first pair 
	    tempH=pairlist[prow][0];
	    tempO1=pairlist[prow][1];
            tempO2=pairlist[prow][2];
	  
//	    printf("%d %d %d\n",tempH,tempO1,tempO2);
	    for(srow=0; srow<slines; srow++){
	       if(switchlist[srow][0]!=0){
	          if(switchlist[srow][0]==tempH){ //If shared proton matches 
 	             if(switchlist[srow][1]==tempO1 && switchlist[srow][2]==tempO2){ //if O1 and O2 matches	
		        if(switchlist[srow-1][1]==tempO2 && switchlist[srow-1][2]==tempO1){ //if previously bound to diff. O atom
		
			    diff=(switchlist[srow][3]-switchlist[srow-1][4])-1; //diff between time PT and then formed Zundel

			     if(diff >= cutoff){ //If proton resides on new partner for longer than cutoff time, successful PT

			      fprintf(output,"%d %d %d %d %d %d\n",
			      switchlist[srow][0], //H index
			      switchlist[srow][1], //current and new O index after transfer
			      switchlist[srow][2], //original O index before transfer
			      switchlist[srow-1][4]+1, //time when proton tranferred
			      switchlist[srow][3]-1, //snap befoe it formed Zundel
			      diff); //persistence

			      switchlist[srow][0]=0;

/*		           printf("Success: diff %d - %d = %d pair:%d %d %d previous: %d %d\n",
			   switchlist[srow][3],
			   switchlist[srow-1][4],
			   diff,
		           switchlist[srow][0],
		           switchlist[srow][1],
		           switchlist[srow][2],
		           switchlist[srow-1][1],
		           switchlist[srow-1][2]);	
*/
		           
			      }
		              else if(diff < cutoff){ //Unsuccessful PT
			         switchlist[srow][0]=0;

/*		           printf("Fail: %d - %d = %d pair: %d %d %d previous: %d %d\n",
			   switchlist[srow][3],
			   switchlist[srow-1][4],
			   diff,
		           switchlist[srow][0],
		           switchlist[srow][1],
		           switchlist[srow][2],
		           switchlist[srow-1][1],
		           switchlist[srow-1][2]);	
*/
			      
		            }		
	                 }//end previous O1 and O2 matches
		       }//end O1 and O2 matches
	           }//end if H index matches
	       }//end nonzero element
            }//end srow-loop
//	    tempH=pairlist[prow][0];
//	    tempO1=pairlist[prow][1];
//            tempO2=pairlist[prow][2];

	 }
	 else if(prow>0){
	    tempH=pairlist[prow][0];
	    tempO1=pairlist[prow][1];
            tempO2=pairlist[prow][2];
//	    printf("%d %d %d\n",tempH,tempO1,tempO2);

	    for(srow=0; srow<slines; srow++){
	       if(switchlist[srow][0]!=0){
	          if(switchlist[srow][0]==tempH){ //If shared proton matches 
 	             if(switchlist[srow][1]==tempO1 && switchlist[srow][2]==tempO2){ //if O1 and O2 matches	
		        if(switchlist[srow-1][1]==tempO2 && switchlist[srow-1][2]==tempO1){ //if previously bound to diff. O atom
		
			    diff=(switchlist[srow][3]-switchlist[srow-1][4])-1; //diff between time PT and then formed Zundel

			     if(diff >= cutoff){ //If proton resides on new partner for longer than cutoff time, successful PT

			      fprintf(output,"%d %d %d %d %d %d\n",
			      switchlist[srow][0], //H index
			      switchlist[srow][1], //current and new O index after transfer
			      switchlist[srow][2], //original O index before transfer
			      switchlist[srow-1][4]+1, //time when proton tranferred
			      switchlist[srow][3]-1, //snap befoe it formed Zundel
			      diff); //persistence

			      switchlist[srow][0]=0;

/*		           printf("Success: diff %d - %d = %d pair:%d %d %d previous: %d %d\n",
			   switchlist[srow][3],
			   switchlist[srow-1][4],
			   diff,
		           switchlist[srow][0],
		           switchlist[srow][1],
		           switchlist[srow][2],
		           switchlist[srow-1][1],
		           switchlist[srow-1][2]);	
*/
		           
			}
		        else if(diff < cutoff){ //Unsuccessful PT
			   switchlist[srow][0]=0;

/*		           printf("Fail: %d - %d = %d pair: %d %d %d previous: %d %d\n",
			   switchlist[srow][3],
			   switchlist[srow-1][4],
			   diff,
		           switchlist[srow][0],
		           switchlist[srow][1],
		           switchlist[srow][2],
		           switchlist[srow-1][1],
		           switchlist[srow-1][2]);	
*/
			}
		     }		
	          }//end if O1 and O2 matches
	       }//end if H index matches
	    }//end nonzero element
         }//end srow-loop
      }//end if prow>1
   }//end prow-loop 


//Initialize matrix and counters for next snapshot
   for(y=0; y<10000; y++){
      for(z=0; z<6; z++){
	 switchlist[y][z]=0;
      }
   }//end y-loop
   slines=0;

   }//end proton-loop
   fclose(output);
        
}//end main -loop
