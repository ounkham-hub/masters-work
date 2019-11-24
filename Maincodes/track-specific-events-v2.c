#include<stdio.h>
#include<stdlib.h>

//******************************************************************************
//Program created to track a specific proton in time output all
//associated structure changes especially 21 or Zundel to Eigen 
//flag types.
//
//An additional output will contain information on 
//successful PTs along with corresponding replica #'s. This will allow us to 
//determine which replicas have or have not undergone PT.
//
//To use this code, you need to know the H index and the two O atoms that the 
//proton is oscillating between. Time frames selected based on previous analysis.
//In the case of H106, O2 and O38 the time frame selected was 150 snapshots 
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
//Last modified: January 31, 2018
//
//****************************************************************************

int main(){

//Parameters to read in input file
   int type,startsnap,nextsnap,Oindex1;
   int startH1,startH2,startH3,startZH;
   int nextH1,nextH2,nextH3,nextZH,Oindex2;
   int y,z,row,coln,lines,bead,proton;
   int initialtimeframe,finaltimeframe;
   int newEigenOindex;
   int origlist[350000][13];
   FILE *input,*output,*output2,*output3;
   char ifile[100],ofile[100];

//CHANGE ACCORDINGLY
   proton=106; 
   newEigenOindex=38;
   initialtimeframe=230; //snapshot 230 to 530
   finaltimeframe=530;


   output2=fopen("all-rep-H106-snaps-230-530.txt","a");
   output3=fopen("all-ZtoE-H106-snaps-230-530","a");
   for(bead=1; bead<=32; bead++){
      sprintf(ifile,"%d-event-list.txt",bead);
      sprintf(ofile,"H%d-rep%d-PTevents.txt",proton,bead); 
      input=fopen(ifile,"r");
      output=fopen(ofile,"w");
  
      fprintf(output2,"\nProton transfer events for Replica %d\n\n",bead);
      fprintf(output3,"\nZundel to Eigen or Successful PT events for Replica %d\n\n",bead);

//Initialize matrices and counters 
      for(y=0; y<350000; y++){
         for(z=0; z<14; z++){
            origlist[y][z]=0;
         }//end z-loop
      }//end y-loop
      lines=0;

//Store input file into an independent matrix,origlist
      while(fscanf(input,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
      &type,&startsnap,&nextsnap,&Oindex1,&startH1,&startH2,&startH3,
      &Oindex2,&nextH1,&nextH2,&nextH3,&startZH,&nextZH)==13){
         origlist[lines][0]=type;
         origlist[lines][1]=startsnap;
         origlist[lines][2]=nextsnap;
         origlist[lines][3]=Oindex1;
         origlist[lines][4]=startH1;
         origlist[lines][5]=startH2;
         origlist[lines][6]=startH3;
	 origlist[lines][7]=Oindex2;
         origlist[lines][8]=nextH1;
         origlist[lines][9]=nextH2;
         origlist[lines][10]=nextH3;;
         origlist[lines][11]=startZH;
         origlist[lines][12]=nextZH;
         lines++;
      }//end while-loop
      fclose(input);

//*********************************************************************
//Section to identify all events associated with specified H atom
//*********************************************************************

      for(row=0; row<lines; row++){
         for(coln=5; coln<14; coln++){
            if(origlist[row][coln]==proton){//If SPECIFIED H index is identified
	       if(origlist[row][0]==2){//if water to zundel transfer
	          if(origlist[row][12]==proton){ //if shared H inzundel 
                     fprintf(output,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                     origlist[row][0],
	             origlist[row][1],
	             origlist[row][2],
	             origlist[row][3],
                     origlist[row][4],
	             origlist[row][5],
	             origlist[row][6],
	             origlist[row][7],
	             origlist[row][8],
	             origlist[row][9],
	             origlist[row][10],
	             origlist[row][11],
	 	     origlist[row][12]);

	             if(origlist[row][1]>=initialtimeframe && origlist[row][1]<=finaltimeframe){
                        fprintf(output2,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			bead,
                        origlist[row][0],
	                origlist[row][1],
	                origlist[row][2],
	                origlist[row][3],
                        origlist[row][4],
	                origlist[row][5],
	                origlist[row][6],
	                origlist[row][7],
	                origlist[row][8],
	                origlist[row][9],
	                origlist[row][10],
	                origlist[row][11],
	 	        origlist[row][12]);
		     }

		     for(y=0; y<13; y++){
		        origlist[row][y]=0;
		     }   
		  }
	       }
	       else if(origlist[row][0]==20){//Zundel to water transfer
	          if(origlist[row][11]==proton){
	             fprintf(output,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                     origlist[row][0],
	             origlist[row][1],
	             origlist[row][2],
	             origlist[row][3],
	             origlist[row][4],
	             origlist[row][5],
	             origlist[row][6],
	             origlist[row][7],
	             origlist[row][8],
	             origlist[row][9],
	             origlist[row][10],
	             origlist[row][11],
		     origlist[row][12]);

	             if(origlist[row][1]>=initialtimeframe && origlist[row][1]<=finaltimeframe){
                        fprintf(output2,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			bead,
                        origlist[row][0],
	                origlist[row][1],
	                origlist[row][2],
	                origlist[row][3],
                        origlist[row][4],
	                origlist[row][5],
	                origlist[row][6],
	                origlist[row][7],
	                origlist[row][8],
	                origlist[row][9],
	                origlist[row][10],
	                origlist[row][11],
	 	        origlist[row][12]);
		     }


		     for(y=0; y<13; y++){
		        origlist[row][y]=0;
		     }
	          }
	       } 
	       else if(origlist[row][0]==12){//Eigen to Zundel transfer
	          if(origlist[row][12]==proton){
	             fprintf(output,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                     origlist[row][0],
	             origlist[row][1],
	             origlist[row][2],
	             origlist[row][3],
                     origlist[row][4],
	             origlist[row][5],
	             origlist[row][6],
	             origlist[row][7],
	             origlist[row][8],
	             origlist[row][9],
	             origlist[row][10],
	             origlist[row][11],
	   	     origlist[row][12]);

	             if(origlist[row][1]>=initialtimeframe && origlist[row][1]<=finaltimeframe){
                        fprintf(output2,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			bead,
                        origlist[row][0],
	                origlist[row][1],
	                origlist[row][2],
	                origlist[row][3],
                        origlist[row][4],
	                origlist[row][5],
	                origlist[row][6],
	                origlist[row][7],
	                origlist[row][8],
	                origlist[row][9],
	                origlist[row][10],
	                origlist[row][11],
	 	        origlist[row][12]);
		     }

		     for(y=0; y<13; y++){
		        origlist[row][y]=0;
		     }
	          }
	       }
	       else if(origlist[row][0]==21){//Zundel to Eigen transfer
	          if(origlist[row][11]==proton){
	             fprintf(output,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                     origlist[row][0],
	             origlist[row][1],
	             origlist[row][2],
	             origlist[row][3],
	             origlist[row][4],
	             origlist[row][5],
	             origlist[row][6],
	             origlist[row][7],
	             origlist[row][8],
	             origlist[row][9],
	             origlist[row][10],
	             origlist[row][11],
	 	     origlist[row][12]);

	             if(origlist[row][1]>=initialtimeframe && origlist[row][1]<=finaltimeframe){
                        fprintf(output2,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			bead,
                        origlist[row][0],
	                origlist[row][1],
	                origlist[row][2],
	                origlist[row][3],
                        origlist[row][4],
	                origlist[row][5],
	                origlist[row][6],
	                origlist[row][7],
	                origlist[row][8],
	                origlist[row][9],
	                origlist[row][10],
	                origlist[row][11],
	 	        origlist[row][12]);

		        if(origlist[row][3]==38){ //If Zundel dissociated to O38 then a transfer occurred!
                           fprintf(output3,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			   bead,
                           origlist[row][0],
	                   origlist[row][1],
	                   origlist[row][2],
	                   origlist[row][3],
                           origlist[row][4],
	                   origlist[row][5],
	                   origlist[row][6],
	                   origlist[row][7],
	                   origlist[row][8],
	                   origlist[row][9],
	                   origlist[row][10],
	                   origlist[row][11],
	 	           origlist[row][12]);
		        }
		     }

		     for(y=0; y<13; y++){
		        origlist[row][y]=0;
		     }
	          }
	       }
	       else if(origlist[row][0]<=10){//All non-zundel type proton transfers
	          if(origlist[row][0]!=2){//exludes water to Zundel transfers, designated as 2
	             fprintf(output,"%d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                     origlist[row][0],
	             origlist[row][1],
	             origlist[row][2],
	             origlist[row][3],
	             origlist[row][4],
	             origlist[row][5],
	             origlist[row][6],
	             origlist[row][7],
	             origlist[row][8],
	             origlist[row][9],
	             origlist[row][10],
	             origlist[row][11],
		     origlist[row][12]);

 	             if(origlist[row][1]>=initialtimeframe && origlist[row][1]<=finaltimeframe){
                        fprintf(output2,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
			bead,
                        origlist[row][0],
	                origlist[row][1],
	                origlist[row][2],
	                origlist[row][3],
                        origlist[row][4],
	                origlist[row][5],
	                origlist[row][6],
	                origlist[row][7],
	                origlist[row][8],
	                origlist[row][9],
	                origlist[row][10],
	                origlist[row][11],
	 	        origlist[row][12]);
		     }

		     for(y=0; y<13; y++){
		        origlist[row][y]=0;
		     }
	          }
	       }
	    }
         }//end coln-loop
      }//end row-loop
      fclose(output);
   }//end bead-loop
   fclose(output2);
   fclose(output3);
}//end main-loop
