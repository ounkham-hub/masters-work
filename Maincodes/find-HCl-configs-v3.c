#include<stdio.h>
#include<stdlib.h>

//*************************************************************
//Program written to identify the over-coordinated oxygens 
//and the associated configurations found in aqueous HCl. The 
//following configurations can be identified:
//Zundels (two water covalently bonded to a proton)
//Eigens (an H3O+ solvated by three waters)
//
//**0.0-1.3 GraphGeod files MUST be organized in ascending order
//**prior to running this program.
//
//Made by Lelee Ounkham and Lance Edens
//Last Updated: 10/25/16
//*************************************************************

void findzundel(int i, int m, int **indexarray, float **distarray);
void findeigen(int i, int m, int n, int **indexarray, int **HBarray);
float **f2d_Matrix(unsigned long n, unsigned long l);
void free_f2d_Matrix(float **m);
int **i2d_Matrix(unsigned long n, unsigned long l);
void free_i2d_Matrix(int **z);

int main ()
{
//Setting up parameters for indexarray
	int l=400;
	int Oindice,Hindice,dum1,dum2,dum3,dum4,dum5;
	float dist,ang;
	int Ocount, j;
//4 variables treated as temp. memory
	int Hhold1,Hhold2,m,p;
	float dhold1,dhold2;
//Setting up parameters for HBarray
	int OHBindice,HHBindice,x,n;
	float HBdist,HBang;
//variables to read input and output
	int i,q;
	char ifile[100],ofile[100],efile[100];
	FILE *input,*output,*einput; // declare a file pointer

	int **indexarray; //sets pointer for indexarray
	int **HBarray;
	float **distarray;

//	input=fopen("zundel-test2.GraphGeod","r");
//	output=fopen("z2.txt","w");
	
//Filling in indice and distance array
	indexarray=i2d_Matrix(l,4); //allocates memory for indexarray
	distarray=f2d_Matrix(l,4);
	HBarray=i2d_Matrix(l,2);	
	
	for(i=33986; i<=33986; i++)
   	{
	   sprintf(ifile,"CB-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",i,i);
           sprintf(ofile,"OC-O-list%d.txt",i);
           input=fopen(ifile,"r"); 
           output=fopen(ofile,"w");
           j=1;
	   m=0;
	   x=0;

        rewind(input);

	   while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oindice,&Hindice,&dum1,&dum2,&dum3,&dum4,&dum5,&dist,&ang)==9)
	   { 		
	       	   if(j==1)
		   { 
	              Ocount=Oindice;
		      Hhold1=Hindice;
		      dhold1=dist;
	           }
	   	   else if(j==2)
		   {
		      if(Ocount==Oindice) //if O is the same on next line
		      {
		         Hhold2=Hindice;
		         dhold2=dist;
		      }
		      else
		      {   
		         j=1;
		         Ocount=Oindice;
		         Hhold1=Hindice;
		         dhold1=dist;
		      }
		   }
		   else if(j==3) //if there is a 3rd occurrence then send to memory/array

		   {
		      if(Ocount==Oindice)
		      {
		         indexarray[m][0]=Oindice;
		         indexarray[m][1]=Hhold1;
		         indexarray[m][2]=Hhold2;
		         indexarray[m][3]=Hindice;
 //converts indice to float
		         distarray[m][0]=(float)Oindice;
  		         distarray[m][1]=(float)dhold1;
		         distarray[m][2]=(float)dhold2;
		         distarray[m][3]=(float)dist;
		         m++;
		         j=0; //restarts Oindice counter which reads next line
		      }
		   	
		      else
                      {   
                         j=1;
                         Ocount=Oindice;
                         Hhold1=Hindice;
                         dhold1=dist;
 		      }
		   }
		   j++; //increments counter of j from 0 to 1 (which is where the loop starts)
	   } //closes while loop

	   fclose(input); 
	
	   for(p=0; p<m; p++) //might need m-1
	   {
	      fprintf(output,"%d %d %d %d %f %f %f\n",indexarray[p][0],indexarray[p][1],indexarray[p][2],indexarray[p][3],distarray[p][1],distarray[p][2],distarray[p][3]);
	   }
	
	fclose(output);
	
	findzundel(i,m,indexarray,distarray); //executes subroutine with the associated information

//Creating HBarray containing O and H index in H-bonded (1.3-2.5) GraphGeod files
	sprintf(efile,"HB-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",i,i);
	einput=fopen(efile,"r");
	n=0;
	 
	while(fscanf(einput,"%d %d %d %d %d %d %d %f %f\n",&OHBindice,&HHBindice,&dum1,&dum2,&dum3,&dum4,&dum5,&HBdist,&HBang)==9)
	{
	   HBarray[n][0]=OHBindice;
	   HBarray[n][1]=HHBindice;
	   n++;
	}
	
	fclose(einput);

	/*for(x=0; x<n; x++)
	{
	   printf("%d %d\n",HBarray[x][0],HBarray[x][1]);
	}
*/
	findeigen(i,m,n,indexarray,HBarray);

	//Zero out arrays after identifying Zundels 
	for(q=0; q<l; q++)
	{
	   indexarray[q][0]=0;
	   indexarray[q][1]=0;
	   indexarray[q][2]=0;
	   indexarray[q][3]=0;
	   distarray[q][0]=0.0;
	   distarray[q][1]=0.0;
	   distarray[q][2]=0.0;
	   distarray[q][3]=0.0;
	}
	   dhold1=0.0;
           dhold1=0.0;
           Hhold1=0;
           Hhold2=0;

	} // closes i-loop

	free_i2d_Matrix(indexarray);
	free_f2d_Matrix(distarray); //removes all previously allocated memory	
	free_i2d_Matrix(HBarray);
	return(0);
}

//Subroutine created to find Zundel structures
void findzundel(int i, int m, int **indexarray, float **distarray) // subroutines go AFTER main()
{
   FILE *Zoutput;
   char Zofile[100];
   int n,p,q,r,z,row1,row2,coln1,coln2;
   int zundelarray[m][4];
     
    
	n=0; //counter for zundelarray

	for(row1=0; row1<m-1; row1++) // n-loop for rows in array because index starts at 0 
	{
	   if(indexarray[row1][0]==0)
	   { 
 	      continue;
	   }
           for(coln1=1; coln1<4; coln1++) // o-loop selects column in array 
	   {
	      for(row2=row1+1; row2<m; row2++)   // p selects the first row after "primary element"
	      {
		for(coln2=1; coln2<4; coln2++) // r selects the second column after "primary element"
	        {
		  if(indexarray[row1][coln1]==indexarray[row2][coln2] && indexarray[row1][coln1] !=0)
 // setting n,o equal to q,r under conditions
		  {
		     for(q=0; q<4; q++) // Column one, oxgen indice
		     {
			
			zundelarray[n][q]=indexarray[row1][q]; //copies entire row1 in new array n
			zundelarray[n+1][q]=indexarray[row2][q]; //n+1 because the partner will always be next to the HB1
			
			indexarray[row1][q]=0; 
			indexarray[row2][q]=0; //inserting zeros into main indexarray
      		     }
			n=n+2; //increments every 2 lines
		  } // end if
		} // end 4
 	      } // end 3
	   } //end for 2
	} //end for 1
	
	sprintf(Zofile,"zundel%d.txt",i);
	Zoutput=fopen(Zofile,"w");


	   for(p=0; p<n; p++) //n starts at 0, n+1 goes to the end of the array 
	   {	
	      fprintf(Zoutput,"%d %d %d %d\n",zundelarray[p][0],zundelarray[p][1],zundelarray[p][2],zundelarray[p][3]);
	   }	
	fclose(Zoutput);
 }

//Subroutine created to find eigen structures	
void findeigen(int i, int m, int n, int **indexarray, int **HBarray)
{

   FILE *Eoutput;
   char Eofile[100];
   int y,z,eigenhit,Enum,row1,HBrow1,coln1;
   int eigenarray[m][4][3];
   int holdarray[3][2];
    	 
    	Enum=0;
	eigenhit=0; //counter for eigen occurrences
	for(y=0; y<3; y++) //zero out holdarray
	    {
	        holdarray[y][0]=0;
	        holdarray[y][1]=0;
	    }
	for(row1=0; row1<m; row1++) // row1-loop for rows in array because index starts at 0 
	{
	   if(indexarray[row1][0]==0)
	   { 
 	      continue;
	   }
           for(coln1=1; coln1<4; coln1++) // coln1-loop selects column in indexarray 
	   {
	      for(HBrow1=0; HBrow1<n; HBrow1++)   // row2 selects the first row after HBarray
	      {
	         if(indexarray[row1][coln1]==HBarray[HBrow1][1] && indexarray[row1][coln1] !=0)
//if the two H's in indexarray and HBarray are equivalent then continue
	         {     		     
		     if(holdarray[coln1-1][0]!=0)
	   	     {
	 		holdarray[coln1-1][1]=HBarray[HBrow1][0];
		     }		      	
		     else
		     {
			holdarray[coln1-1][0]=HBarray[HBrow1][0];
		     }
		     eigenhit++;
			  }//end if
	       }//end HBrow1-loop
	  
	    if(coln1==3)
	    {
	       if(eigenhit==3)
	       {
		  eigenarray[Enum][0][0]=indexarray[row1][0]; //Corresponds to OC-O
		  eigenarray[Enum][1][0]=indexarray[row1][1]; //CBH1 to OC-O
		  eigenarray[Enum][1][1]=holdarray[0][0]; //O1 to HB1 
		  eigenarray[Enum][1][2]=holdarray[0][1]; //O2 to HB1
		  eigenarray[Enum][2][0]=indexarray[row1][2]; //CBH2 to OC-O
		  eigenarray[Enum][2][1]=holdarray[1][0]; //O1 to HB2
		  eigenarray[Enum][2][2]=holdarray[1][1]; //O2 to HB2
		  eigenarray[Enum][3][0]=indexarray[row1][3]; //CBH3 to OC-O 
		  eigenarray[Enum][3][1]=holdarray[2][0]; //O1 to HB3
		  eigenarray[Enum][3][2]=holdarray[2][1]; //O2 to HB3
		  Enum++;
	       }
    	    }
	  }//end coln1-loop
	  eigenhit=0; //zero out counter eigenhit
	  for(y=0; y<3; y++) //zero out holdarray
	    {
	        holdarray[y][0]=0;
	        holdarray[y][1]=0;
	    }
	}//end row1-loop
	
	sprintf(Eofile,"eigen%d.txt",i);
	Eoutput=fopen(Eofile,"w");
//Abbreviations: 
//OC-O (over-coordinated oxygen, has 3 CHBs)
//CBH1,CBH2,CBH3 (covalently-bonded hydrogen)
//HB1O1,HB2O1,HB2O1 (hydrogen-bonded oxygen, typically only 1 but can be 2) 
	fprintf(Eoutput,"OC-O CBH1 HB1O1 CBH2 HB2O1 CBH3 HB3O1\n");
	for(z=0; z<Enum; z++)  
	{	 
	  if(eigenarray[z][0][0]!=0)
	   {
	     fprintf(Eoutput,"%d    %d   %d   ",
	     eigenarray[z][0][0], //OCO 
	     eigenarray[z][1][0], //covalently bonded H1 to OCO
             eigenarray[z][1][1]); //O1 H-bonded to H1
	   }
	   if(eigenarray[z][1][2]!=0)
	   {
	     fprintf(Eoutput,"%d   ",eigenarray[z][1][2]); //O2 H-bonded to H1, if ever
	   }   
	  
	   fprintf(Eoutput,"%d   %d   ",
	   eigenarray[z][2][0], //covalently bonded H2 to OCO 
	   eigenarray[z][2][1]); //O1 H-bonded to H2

	   if(eigenarray[z][2][2]!=0) 
	   {
	      fprintf(Eoutput,"%d   ",eigenarray[z][2][2]); //O2 H-bonded to H2, if ever
	   }	     

	   fprintf(Eoutput,"%d   %d   ",
	   eigenarray[z][3][0], //covalently bonded to H3 to OCO
	   eigenarray[z][3][1]); //O1 H-bonded to H3
	 
	   if(eigenarray[z][3][2]!=0)
	   {
	      fprintf(Eoutput,"%d ",eigenarray[z][3][2]); //O2 H-bonded to H3, if ever
	   }
	   fprintf(Eoutput,"\n");
	}
	
	fclose(Eoutput);	
} 

float **f2d_Matrix(unsigned long n, unsigned long l)
{
        float **m;
	unsigned long i;

	m = (float **)calloc(n, sizeof(float *));
        if(m==NULL){printf("Memory allocation error: m. Exiting.\n"); exit(EXIT_FAILURE);}

	*m = (float *)calloc(n*l, sizeof(float));
	if(*m==NULL){printf("Memory allocation error: *m. Exiting.\n"); exit(EXIT_FAILURE);}

        for(i=1; i<n; i++)
	        m[i] = m[i-1] + l;

        return(m);
}
void free_f2d_Matrix(float **m)
{
        free(*m);
        free(m);
}

int **i2d_Matrix(unsigned long n, unsigned long l)
{
	int **m;
	unsigned long i;

	m = (int **)calloc(n, sizeof(int *));
	if(m==NULL){printf("Memory allocation error: m. Exiting.\n"); exit(EXIT_FAILURE);}
	*m = (int *)calloc(n*l, sizeof(int));
	if(*m==NULL){printf("Memory allocation error: *m. Exiting.\n"); exit(EXIT_FAILURE);}

	for(i=1; i<n; i++)
		m[i] = m[i-1] + l;

	return(m);
}
void free_i2d_Matrix(int **z)
{
        free(*z);
        free(z);
}

