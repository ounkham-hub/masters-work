#include<stdio.h>
#include<stdlib.h>

//*************************************************************
//Program written to identify the over-coordinated oxygens 
//and the associated configurations found in aqueous HCl. The 
//following configurations can be identified:
//Zundels (two water covalently bonded to a proton)
//Eigens (an H3O+ solvated by three waters)
//
//**ALL GraphGeod files MUST be organized in ascending order
//**prior to running this program.
//
//Made by Lelee Ounkham and Lance Edens
//Last Updated: 10/13/16
//*************************************************************

void findzundel(int i, int m, int **indexarray, float **distarray);
float **f2d_Matrix(unsigned long n, unsigned long l);
void free_f2d_Matrix(float **m);
int **i2d_Matrix(unsigned long n, unsigned long l);
void free_i2d_Matrix(int **z);

int main ()
{
//Setting up parameters
	int l=400;
	char chr;
	int Oindice,Hindice,dum1,dum2,dum3,dum4,dum5;
	float dist,ang,dhold1,dhold2;
	int Ocount, j;
//4 variables treated as temp. memory
	int Hhold1,Hhold2,m,p;
//variables to read input and output
	int i,q;
	char ifile[100],ofile[100];
	FILE *input, *output; // declare a file pointer

	int **indexarray; //sets pointer for indexarray
	float **distarray;

//	input=fopen("zundel-test2.GraphGeod","r");
//	output=fopen("z2.txt","w");
	
//Filling in indice and distance array
	indexarray=i2d_Matrix(l,4); //allocates memory for indexarray
	distarray=f2d_Matrix(l,4);
	
	
	for(i=1; i<=31136; i++)
   	{
           sprintf(ifile,"R-HCl.input.O%d.xyz.H%d.xyz.GraphGeod",i,i);
           sprintf(ofile,"OC-O-list%d.txt",i);
           input=fopen(ifile,"r");
           output=fopen(ofile,"w");
           m=0;
           j=1;

        rewind(input);

	   while(fscanf(input,"%d %d %d %d %d %d %d %f %f\n",&Oindice,&Hindice,	&dum1,&dum2,&dum3,&dum4,&dum5,&dist,&ang)==9)
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

	return(0);
}

void findzundel(int i, int m, int **indexarray, float **distarray) // subroutines go AFTER main()
{
   FILE *Zoutput;
   char Zofile[100];
   int n,p,q,r,row1,row2,coln1,coln2;
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

