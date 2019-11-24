#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//**************************************************************************
//Program created to compare weighted networks to networks created by the
//weighing function. This is done by taking the difference between weights
//for specific bonds and count how many bonds are excluded. 
//
//Created by Lelee Ounkham
//Last modified October 2, 2018
//
//**************************************************************************

int main(){
   int snap,Oxy,Hyd,Ox,Hy,flines,nlines,nrow,frow,index,y,z;
   int wncount,wfcount;
   float dist,wnweight,wfweight,diff,wf,wn;
   FILE *WFinput,*WNinput,*output,*output2;
   char WNfile[100],WFfile[100];

   int WNmatrix[1000][3];
   int WFmatrix[1000][3];
   float WFweight[1000];
   float WNweight[1000];


   output=fopen("total-WF-WN-weight-difference.txt","a");
   output2=fopen("WF-WN-disagreement.txt","a");

   for(snap=1; snap<=61555; snap++){
//Initialize matrix
      for(y=0; y<1000; y++){
         for(z=0; z<3; z++){
	    WFmatrix[y][z]=0;
	    WNmatrix[y][z]=0;
            WFweight[y]=0.0;
	    WNweight[y]=0.0;
         }
      }
      flines=0;
      nlines=0;
      wncount=0;
      wfcount=0;

      sprintf(WFfile,"R-avg-O%d.H%d.wGraphGeod",snap,snap); //open current snap
      WFinput=fopen(WFfile,"r");

      while(fscanf(WFinput,"%d %d %d %f\n",&index,&Ox,&Hy,&wfweight)==4){
	 WFmatrix[flines][0]=Ox;
	 WFmatrix[flines][1]=Hy;
	 WFweight[flines]=wfweight;
	 flines++;
      }//end while-loop
      fclose(WFinput);
    
      sprintf(WNfile,"HCl.covalent.O%d.xyz.H%d.xyz.wGraphGeod",snap,snap); //open current snap
      WNinput=fopen(WNfile,"r");

      while(fscanf(WNinput,"%d %d %f %f\n",&Oxy,&Hyd,&dist,&wnweight)==4){
         WNmatrix[nlines][0]=Oxy;
	 WNmatrix[nlines][1]=Hyd;
	 WNweight[nlines]=wnweight;
         nlines++;
      }//end while-loop
      fclose(WNinput);

//Take the difference in weights
      for(nrow=0; nrow<nlines; nrow++){
         for(frow=0; frow<flines; frow++){
	    if(WFmatrix[frow][0]!=0 && WNmatrix[nrow][0]!=0){
	       if(WNmatrix[nrow][0]==WFmatrix[frow][0] && WNmatrix[nrow][1]==WFmatrix[frow][1]){
 //If O and H index matches
 	          diff=WNweight[nrow]-WFweight[frow];
		  fprintf(output,"%d %d %0.3f\n",WNmatrix[nrow][0],WFmatrix[frow][1],diff);
		  WFmatrix[frow][0]=0;
		  WNmatrix[nrow][0]=0;
	       }    	
	    }
	 }//end frow-loop
      }//end nrow-loop

//Output the number of disagreed weights 
      for(nrow=0; nrow<nlines; nrow++){
	if(WNmatrix[nrow][0]!=0){
           wncount++;
	}
      }
 
       for(frow=0; frow<flines; frow++){
	if(WFmatrix[frow][0]!=0){
           wfcount++;
	}
      }       
      printf(output2,"%d %d %d\n",snap,wncount,wfcount); 
   }//end snap-loop
}//end main-loop
