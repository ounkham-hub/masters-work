/*
 * This is the code for post-selecting the O..O graph based on angle criterion
 * 2017.01.09 by Tiecheng Zhou
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define Xsize 14.926
#define Ysize 14.926
#define Zsize 14.926
#define PI 3.14159265

float pbcshift(float ref, float input, char lab)  // considering the periodic boundary condition
{
	float result;
	result=input;

	if(lab == 'x')
	{
		if(input - ref > Xsize/2)
			result = input - Xsize;
		if(ref - input > Xsize/2)
			result = input + Xsize;
	}

	if(lab == 'y')
	{
		if(input - ref > Ysize/2)
			result = input - Ysize;
		if(ref - input > Ysize/2)
			result = input + Ysize;
	}

	if(lab == 'z')
	{
		if(input - ref > Zsize/2)
			result = input - Zsize;
		if(ref - input > Zsize/2)
			result = input + Zsize;
	}

	return(result);
}


int main(int argc, char *argv[])
{
	if(argc != 6)
	{
		printf("Use this script:\n%s GraphGeod-file O-xyz-file H-xyz-file\n\n",argv[0]);
		exit(-1);
	}

	FILE *fp1,*fp2,*fp3;
	if( (fp1=fopen(argv[1],"r"))==NULL )
	{
		printf("Error: can not open %s\n",argv[1]);
		exit(-1);
	}
	if( (fp2=fopen(argv[2],"r"))==NULL )
	{
		printf("Error: can not open %s\n",argv[2]);
		exit(-1);
	}
	if( (fp3=fopen(argv[3],"r"))==NULL )
	{
		printf("Error: can not open %s\n",argv[3]);
		exit(-1);
	}

	float anglow,angup; // angle criterion with lower boundary and upper boundary
	anglow = atof(argv[4]);
	angup = atof(argv[5]);

	printf(" Use the code with angle selection of (%f %f)\n\n",anglow,angup);

	int Graphpairs[10000][2]; // the interaction pairs from GraphGeod file
	float Oxycoord[1000][3]; // the coordinates from Oxygen-xyz-file
	float Hydcoord[1000][3]; // the coordinates from Hydrogen-xyz-file
	int npairs,nOxy,nHyd; // number of paiers, number of Oxygen atoms, number of Hydrogen atoms

	char buffer[1000];
	int node1,node2;

	rewind(fp1);
	npairs=0;
	while( fscanf(fp1,"%d %d",&node1,&node2) == 2) // read the molecular index, (the index of Oxygen atoms)
	{
		npairs++;
		Graphpairs[npairs][0] = node1;
		Graphpairs[npairs][1] = node2;
//		printf("%d %d %d\n",npairs,node1,node2);

		fgets(buffer,sizeof(buffer),fp1);
	}
	fclose(fp1);

	char label,*tok;
	float tempx,tempy,tempz;

	rewind(fp2);
	nOxy=0;
	fgets(buffer,sizeof(buffer),fp2);
	fgets(buffer,sizeof(buffer),fp2);   // skip the two head lines in xyz-file
	while( fgets(buffer,sizeof(buffer),fp2) != NULL ) // read the xyz coordinates from Oxygen-xyz-file
	{
		tok = strtok(buffer," "); // break the whole line using " " as separation, the first one is atomic label
		label = tok[0];
		tok = strtok(NULL, " ");
		tempx = atof(tok);
		tok = strtok(NULL, " ");
		tempy = atof(tok);
		tok = strtok(NULL, " ");
		tempz = atof(tok);

		nOxy++;
		Oxycoord[nOxy][0] = tempx;
		Oxycoord[nOxy][1] = tempy;
		Oxycoord[nOxy][2] = tempz;
//		printf(" %d %f %f %f\n",nOxy,tempx,tempy,tempz);

	}
	fclose(fp2);

	rewind(fp3);
	nHyd=0;
	fgets(buffer,sizeof(buffer),fp3);
	fgets(buffer,sizeof(buffer),fp3);   // skip the two head lines in xyz-file
	while( fgets(buffer,sizeof(buffer),fp3) != NULL ) // read the xyz coordinates from Hydrogen-xyz-file
	{
		tok = strtok(buffer," "); // break the whole line using " " as separation, the first one is atomic label
		label = tok[0];
		tok = strtok(NULL, " ");
		tempx = atof(tok);
		tok = strtok(NULL, " ");
		tempy = atof(tok);
		tok = strtok(NULL, " ");
		tempz = atof(tok);

		nHyd++;
		Hydcoord[nHyd][0] = tempx;
		Hydcoord[nHyd][1] = tempy;
		Hydcoord[nHyd][2] = tempz;
//		printf(" %d %f %f %f\n",nHyd,tempx,tempy,tempz);

	}
	fclose(fp3);


//	printf("%d %d %d\n",npairs,nOxy,nHyd);

	// end of reading input files, the following is calculating the angle and generating the new GraphGeod file
	FILE *fop;
	char filename[1000];
	strcpy(filename,argv[1]);
	strcat(filename,"-updated");  // the new name of the GraphGeod file is with "-updated"
	fop=fopen(filename,"w");

	int i,j,nodei,nodej;
	float O1x,O1y,O1z,O2x,O2y,O2z; // the coordinates of the two Oxygen atoms
	float Hx,Hy,Hz; // the coordinates of the Hydrogen atom
	float distOO, distO1H, distO2H, cosang, angle; // the distance between O1..O2, O1..H, O2..H, and the O1.H.O2 angle

	for(i=1;i<=npairs;i++)
	{
		nodei = Graphpairs[i][0];
		nodej = Graphpairs[i][1];
		O1x = Oxycoord[nodei][0];
		O1y = Oxycoord[nodei][1];
		O1z = Oxycoord[nodei][2];
		O2x = pbcshift(O1x,Oxycoord[nodej][0],'x');
		O2y = pbcshift(O1y,Oxycoord[nodej][1],'y');
		O2z = pbcshift(O1z,Oxycoord[nodej][2],'z');
		distOO = sqrt( pow(O1x-O2x,2) + pow(O1y-O2y,2) + pow(O1z-O2z,2) );

//		printf("%d %d: %f\n",nodei,nodej,distOO);

		for(j=1;j<=nHyd;j++)
		{
			Hx = pbcshift(O1x,Hydcoord[j][0],'x');
			Hy = pbcshift(O1y,Hydcoord[j][1],'y');
			Hz = pbcshift(O1z,Hydcoord[j][2],'z');

			distO1H = sqrt( pow(O1x-Hx,2) + pow(O1y-Hy,2) + pow(O1z-Hz,2) );
			distO2H = sqrt( pow(O2x-Hx,2) + pow(O2y-Hy,2) + pow(O2z-Hz,2) );

			if(distO1H < 3.0 && distO2H < 3.0) // rough selection, make sure the H atom is near the two O atoms
			{
				cosang = (pow(distO1H,2) + pow(distO2H,2) - pow(distOO,2)) / (2*distO1H*distO2H);   // cosine theorem to calculate the angle
				angle = acos(cosang) * 180.0 / PI;

				if(angle >= anglow && angle <= angup)
					fprintf(fop, "%d %d %d %f %f %f %f\n",nodei,nodej,j,distOO,distO1H,distO2H,angle);

			}
		}
 
	}

	fclose(fop);

	return(0);
}




