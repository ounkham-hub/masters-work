#include<stdio.h>
#include<stdlib.h>

//*****************************************************************
//Program created to output O+-Ow, Ow-Ow and Oz-Oz distances.
//
//Created by Lelee Ounkham
//Last Modified on December 4, 2017
//
//*****************************************************************  

int main(){
//Setting up parameters
   int Oindex,O1index,O2index,Hindex,PBC1,PBC2,PBC3,PBC4,PBC5;
   int Zindex,pairno;
   char label;
   float dist,ang,xcoord,ycoord,zcoord;
   int snap,elines,olines,wlines,zlines;
   int erow,orow,wrow,zrow,y,z;
   int ZGraphGeod[300][2];
   int EGraphGeod[300];
   int WGraphGeod[300];
   int OOGraphGeod[300][2];
   float OOdistance[300];
   FILE *einput,*input,*zinput;
   FILE *eoutput,*zoutput,*woutput;
   char efile[100],zfile[100],ifile[100],ofile[100];

   eoutput=fopen("Eigen-OO-distances.txt","a");
   zoutput=fopen("Zundel-OO-distances.txt","a");
   woutput=fopen("Water-OO-distances.txt","a");

   for(snap=1; snap<=55001; snap++){
      sprintf(zfile,"Zpair-OO-%d.xyz",snap);
      sprintf(ifile,"HCl.input.O%d.xyz.O%d.xyz.GraphGeod",snap,snap);
      sprintf(efile,"E-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",snap,snap);
      einput=fopen(efile,"r");
      input=fopen(ifile,"r");
      zinput=fopen(zfile,"r");

      elines=0;
      olines=0;
      wlines=0;
      zlines=0;

      for(y=0; y<=300; y++){
	 EGraphGeod[y]=0;
	 WGraphGeod[y]=0;
	 ZGraphGeod[y][0]=0;
         ZGraphGeod[y][1]=0;
	 OOGraphGeod[y][0]=0;
	 OOGraphGeod[y][1]=0;
	 OOdistance[y]=0.0;
      }

      while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&O1index,&O2index,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&ang)==9){
         OOGraphGeod[olines][0]=O1index;
	 OOGraphGeod[olines][1]=O2index;
	 OOdistance[olines]=dist;
         olines++;

//         printf("%d %d\n",O1index,O2index);
      }//end while-loop
      fclose(input);

      while(fscanf(zinput,"%d %s %f %f %f %d\n",&Zindex,&label,&xcoord,&ycoord,&zcoord,&pairno)==6){
         if(pairno!=77){
            ZGraphGeod[zlines][0]=Zindex;
	    ZGraphGeod[zlines][1]=pairno;
            zlines++;

//            printf("%d %d\n",Zindex,pairno);
	 }
	 else if(pairno==77){
	    WGraphGeod[wlines]=Zindex;
	    wlines++;
	 }
      }//end while-loop
      fclose(zinput);

//**************************************************************
//Zundels:
//Section to output O-O distances for H3O+ to H3O+
//
//**************************************************************

      for(zrow=0; zrow<zlines; zrow++){
         for(orow=0; orow<olines; orow++){
	    if(OOGraphGeod[orow][0]!=0){
               if(ZGraphGeod[zrow][0]==OOGraphGeod[orow][0]){//if ZO1 matches O1 in O1-O2 GraphGeod
		  for(z=0; z<zlines; z++){
		     if(ZGraphGeod[z][0]!=ZGraphGeod[zrow][0]){ //If O index DIFFERS from ZO1
			if(ZGraphGeod[z][1]==ZGraphGeod[zrow][1]){ //and pairno matches, Zundel identified
                           if(OOGraphGeod[orow][1]==ZGraphGeod[z][0]){
  	                      fprintf(zoutput,"%d %d %f\n",
	                      OOGraphGeod[orow][0],
	                      OOGraphGeod[orow][1],
	                      OOdistance[orow]);

	                      OOGraphGeod[orow][0]=0;
	                      OOGraphGeod[orow][1]=0;
	                      OOdistance[orow]=0.0;
			   }
		        }
		     }
		  }//end z-loop
	       }
	    }
	 }//end orow-loop
      }//end zrow-loop

      while(fscanf(einput,"%d %d %d %d %d %d %d %f %f\n",&Oindex,&Hindex,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&dist,&ang)==9){
         EGraphGeod[elines]=Oindex;
         elines++;
//printf("%d\n",Oindex);
      }//end while-loop
      fclose(einput);      

//**************************************************************
//Eigens:
//Section to output O-O distances for H3O+ to H2O
//
//**************************************************************

      for(erow=0; erow<elines; erow++){
	 for(orow=0; orow<olines; orow++){
	    if(OOGraphGeod[orow][0]!=0){
	       if(EGraphGeod[erow]==OOGraphGeod[orow][0]){ //if O index from Eigen GG matches O1 from OO GG
	          fprintf(eoutput,"%d %d %f\n",
	          OOGraphGeod[orow][0],
	          OOGraphGeod[orow][1],
	          OOdistance[orow]);

	          OOGraphGeod[orow][0]=0;
	          OOGraphGeod[orow][1]=0;
	          OOdistance[orow]=0.0;
	       }
               else if(EGraphGeod[erow]==OOGraphGeod[orow][1]){
	          fprintf(eoutput,"%d %d %f\n",
	          OOGraphGeod[orow][0],
	          OOGraphGeod[orow][1],
	          OOdistance[orow]);

	          OOGraphGeod[orow][0]=0;
	          OOGraphGeod[orow][1]=0;
	          OOdistance[orow]=0.0;
	       }
            }
         }//end orow-loop
      }//end erow-loop

//**************************************************************
//Waters:
//Section to output O-O distances for H2O to H2O
//
//**************************************************************

      for(wrow=0; wrow<wlines; wrow++){
	 for(orow=0; orow<olines; orow++){
            if(OOGraphGeod[orow][0]!=0){
	       if(WGraphGeod[wrow]==OOGraphGeod[orow][0]){
	          fprintf(woutput,"%d %d %f\n",
	          OOGraphGeod[orow][0],
	          OOGraphGeod[orow][1],
	          OOdistance[orow]);
	       }
	    }
	 }//end orow-loop
      }//end wrow-loop

   } //end snap-loop
   fclose(eoutput);
   fclose(zoutput);
   fclose(woutput);  
}//end main

