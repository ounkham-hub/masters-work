#include<stdio.h>
#include<stdlib.h>

//Program created to read and output multiple files. 
//Code works for the following example.

int numfiles=3;
int i, j, line;
int MAXSIZE=100;
FILE *input;
FILE *output;
char ifile[50],ofile[50];

int main()
{

  for(i=1; i<=3; i++)
  {
	printf("%d\n",i);         
       sprintf(ifile,"test%d.txt",i); 
//Because i goes from 1-3, test0.txt to test3.txt will be written into i files
       sprintf(ofile, "output%d.txt",i);
      input=fopen(ifile,"r");
       output=fopen(ofile,"w");
	printf("%s\n",ifile);
	printf("%s\n",ofile);
      
       fscanf(input,"%d",&line);
       fprintf(output,"%d\n",line);
  
      fclose(input);
      fclose(output);     
 
  } //end of for loop
	return(0);
}
