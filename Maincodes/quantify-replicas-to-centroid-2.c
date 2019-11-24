#include<stdio.h>
#include<stdlib.h>


//**************************************************************************
//Program created to quantify how many replicas per O-H bond agree with
//the centroid dataset. 
//
//Input: R-HCl.input.solB#.xyz.solC#.xyz.GraphGeod
//       count-32rep-snapshot-#.txt
//Output:
//
//Created by Lelee Ounkham
//Last modified on July 28, 2018
//
//**************************************************************************


int main(){
  int Oindex,Hindex,PBC1,PBC2,PBC3,PBC4,PBC5;
  int Oatom,Hatom,count,clines,rlines,snap;
  int countones,countzeros,x,y,crow,row;
  float distance,angle,dist; 
  FILE *gginput,*cinput,*coutput,*routput;
  char gfile[100],cfile[100];
  
  int centroidmatrix[300][3];
  int replicamatrix[300][4];
  float centdist[300][1];
  float repdist[300][1];

  clines=0;
  rlines=0; 
  countones=0;
  countzeros=0;

  coutput=fopen("centroid-num-beads.txt","a");
  routput=fopen("replica-num-beads.txt","a");

  for(snap=1; snap<=61555; snap++){
      sprintf(gfile,"R-HCl.input.solB%d.xyz.solC%d.xyz.GraphGeod",snap,snap);	
      sprintf(cfile,"dist-count-32rep-%d.txt",snap);
      gginput=fopen(gfile,"r");
      cinput=fopen(cfile,"r");

      while(fscanf(gginput,"%d %d %d %d %d %d %d %f %f\n",&Oindex,&Hindex,&PBC1,&PBC2,&PBC3,&PBC4,&PBC5,&distance,&angle)==9){
	  centroidmatrix[clines][0]=Oindex;
	  centroidmatrix[clines][1]=Hindex;
	  centdist[clines][0]=distance;

	  clines++;
      }
      fclose(gginput);

      while(fscanf(cinput,"%d %d %d %f\n",&count,&Oatom,&Hatom,&dist)==4){
	  replicamatrix[rlines][0]=Oatom;
	  replicamatrix[rlines][1]=Hatom;
	  replicamatrix[rlines][2]=count;

	  repdist[rlines][0]=dist;
	  rlines++;
      }
      fclose(cinput);

      for(crow=0; crow<clines; crow++){
         for(row=0; row<rlines; row++){
	    if(centroidmatrix[crow][0]==replicamatrix[row][0] && centroidmatrix[crow][1]==replicamatrix[row][1]){
	       if(centdist[x][0]==repdist[x][0]){ //If distances matches
	 	  centroidmatrix[crow][2]=1;
  		  replicamatrix[row][3]=1; //turn on flags

	           fprintf(coutput,"%d %d %d %f %d\n",
		   snap,
		   centroidmatrix[crow][0], //Oindex
		   centroidmatrix[crow][1], //Hindex
		   centdist[crow][0], //distance
		   replicamatrix[row][2]); //# of replicas 
		   countones++;
               }
	    }
	 }//end row-loop
      }//end crow-loop 

      for(crow=0; crow<clines; crow++){
	 if(centroidmatrix[crow][2]==0){
	    countzeros++;
	 }
      }//end crow-loop
 
      for(row=0; row<rlines; row++){
	 if(replicamatrix[row][3]==0){
	    fprintf(routput,"%d %d %d %f %d\n",
	    snap,
	    replicamatrix[row][0],
	    replicamatrix[row][1],
	    repdist[row][0],
	    replicamatrix[row][2]);
	 }
      }

//      printf("centroid: snapshot %d agree: %d disagree: %d out of %d\n",snap,countones,countzeros,clines);
//      printf("replica: snapshot %d agree: %d disagree: %d out of %d\n",snap,countones,rlines-countones,rlines);


//      fprintf(coutput,"%d %d %d %d\n",snap,countones,countzeros,clines);
//      fprintf(routput,"%d %d %d %d\n",snap,countones,rlines-countones,rlines);

//Initialize centroidmatrix and replicamatrix
      for(x=0; x<clines; x++){
         for(y=0; y<3; y++){
	    centroidmatrix[x][y]=0;
	    centdist[x][0]=0.0;
	 }
      }//end x-loop

      for(x=0; x<rlines; x++){
         for(y=0; y<4; y++){
	    replicamatrix[x][y]=0;
	    repdist[x][0]=0.0;
	 }
      }//end x-loop

      clines=0;
      rlines=0;
      countzeros=0;
      countones=0; 

  }//end snap-loop
  fclose(coutput);
  fclose(routput);
}//end main-loop
