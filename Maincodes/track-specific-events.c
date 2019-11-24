#include<stdio.h>
#include<stdlib.h>

//*********************************************************************
//Program created to track a specific proton in time output all
//associated structure changes especially 21 or Zundel to Eigen 
//flag types.
//
//Note the flags indicate the following:
//    - 0, water
//    - 1, eigen
//    - 2, zundel
//
//Input: #-event-list.txt
//Output: H#-PT-events.txt
//
//Created by Lelee Ounkham
//Last modified: January 29, 2018
//
//*********************************************************************


int main(){

//Parameters to read in input file
   int type,startsnap,nextsnap,Oindex;
   int startH1,startH2,startH3,startZH;
   int nextH1,nextH2,nextH3,nextZH;
   int y,z,row,coln,lines,bead,proton;
   int origlist[3000][12];
   FILE *input;
   char ifile[100];

     
   for(bead=32; bead<=32; bead++){
      sprintf(ifile,"%d-event-list.txt",bead);
      input=fopen(ifile,"r");

//Initialize matrices and counters 
      for(y=0; y<3000; y++){
         for(z=0; z<13; z++){
            origlist[y][z]=0;
         }//end z-loop
      }//end y-loop
      lines=0;
      proton=106; //CHANGE ACCORDINGLY

//Store input file into an independent matrix,origlist
      while(fscanf(input,"%d %d %d %d %d %d %d %d %d %d %d %d\n",
      &type,&startsnap,&nextsnap,&Oindex,&startH1,&startH2,&startH3,
      &nextH1,&nextH2,&nextH3,&startZH,&nextZH)==12){
         origlist[lines][0]=type;
         origlist[lines][1]=startsnap;
         origlist[lines][2]=nextsnap;
         origlist[lines][3]=Oindex;
         origlist[lines][4]=startH1;
         origlist[lines][5]=startH2;
         origlist[lines][6]=startH3;
         origlist[lines][7]=nextH1;
         origlist[lines][8]=nextH2;
         origlist[lines][9]=nextH3;;
         origlist[lines][10]=startZH;
         origlist[lines][11]=nextZH;
         lines++;
      }//end while-loop
   }//end bead-loop

//*********************************************************************
//Section to identify all events associated with specified H atom
//*********************************************************************

   for(row=0; row<lines; row++){
      for(coln=5; coln<12; coln++){
         if(origlist[row][coln]==proton){//If SPECIFIED H index is identified
	    if(origlist[row][0]==2){//if water to zundel transfer
	       if(origlist[row][11]==proton){ //if shared H inzundel 
                  printf("%d %d %d %d %d %d %d %d %d %d %d %d\n",
                  origlist[row][0],
	          origlist[row][1],
	          origlist[row][2],
	          origlist[row][3],
	          origlist[row][4],
                  origlist[row][4],
	          origlist[row][5],
	          origlist[row][6],
	          origlist[row][7],
	          origlist[row][8],
	          origlist[row][9],
	          origlist[row][10],
	          origlist[row][11]);

		  for(y=0; y<12; y++){
		     origlist[row][y]=0;
		  }
	       }
	    }
	    else if(origlist[row][0]==20){//Zundel to water transfer
	       if(origlist[row][10]==106){
	          printf("%d %d %d %d %d %d %d %d %d %d %d %d\n",
                  origlist[row][0],
	          origlist[row][1],
	          origlist[row][2],
	          origlist[row][3],
	          origlist[row][4],
                  origlist[row][4],
	          origlist[row][5],
	          origlist[row][6],
	          origlist[row][7],
	          origlist[row][8],
	          origlist[row][9],
	          origlist[row][10],
	          origlist[row][11]);

		  for(y=0; y<12; y++){
		     origlist[row][y]=0;
		  }
	       }
	    }
	    else if(origlist[row][0]==12){//Eigen to Zundel transfer
	       if(origlist[row][11]==106){
	          printf("%d %d %d %d %d %d %d %d %d %d %d %d\n",
                  origlist[row][0],
	          origlist[row][1],
	          origlist[row][2],
	          origlist[row][3],
	          origlist[row][4],
                  origlist[row][4],
	          origlist[row][5],
	          origlist[row][6],
	          origlist[row][7],
	          origlist[row][8],
	          origlist[row][9],
	          origlist[row][10],
	          origlist[row][11]);

		  for(y=0; y<12; y++){
		     origlist[row][y]=0;
		  }
	       }
	    }
	    else if(origlist[row][0]==21){//Zundel to Eigen transfer
	       if(origlist[row][10]==106){
	          printf("%d %d %d %d %d %d %d %d %d %d %d %d\n",
                  origlist[row][0],
	          origlist[row][1],
	          origlist[row][2],
	          origlist[row][3],
	          origlist[row][4],
                  origlist[row][4],
	          origlist[row][5],
	          origlist[row][6],
	          origlist[row][7],
	          origlist[row][8],
	          origlist[row][9],
	          origlist[row][10],
	          origlist[row][11]);

		  for(y=0; y<12; y++){
		     origlist[row][y]=0;
		  }
	       }
	    }
	    else if(origlist[row][0]<=10){//All non-zundel type proton transfers
	       if(origlist[row][0]!=2){//exludes water to Zundel transfers, designated as 2
	             printf("%d %d %d %d %d %d %d %d %d %d %d %d\n",
                     origlist[row][0],
	             origlist[row][1],
	             origlist[row][2],
	             origlist[row][3],
	             origlist[row][4],
                     origlist[row][4],
	             origlist[row][5],
	             origlist[row][6],
	             origlist[row][7],
	             origlist[row][8],
	             origlist[row][9],
	             origlist[row][10],
	             origlist[row][11]);
 
		  for(y=0; y<12; y++){
		     origlist[row][y]=0;
		  }
	       }
	    }
	 }
      }//end coln-loop
   }//end row-loop

}//end main-loop
