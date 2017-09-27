/* This is the code to calculate pagerank from the trajectory files
 * it can read the xyz trajectory files, construct the graphs and calculate the pagerank
 *    the graph is constructed by Al center with its neighbors within a cutoff shell
 * Before running it, make sure the Box size is correct, and the number of atoms does not exceed the array size
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define FLN 1000
#define Size 10000   // maximum number of atoms in the system
#define Boxx 14.3   // system box size of this 'test.xyz' is 14.3 Ang
#define Boxy 14.3
#define Boxz 14.3
#define Factor 1.0  // the factor in fermi function or power function, make pagerank smooth
#define CNfactor 20  // this is used in determining coordination number
#define Coordcut 3.3 // the cutoff range in definiting coordination number of Na, I rename Na as Al (Jun 22,2017)
#define CovOHcut 1.3 // the cutoff of covalent OH bond, used to define charge state
#define Convtol 0.00001  // the tolerance of convergence in pagerank calculation


// the functions

// considering pbc, shift the coordinates
double pbcshift(double inp, double center, char dim);
// calculate the coordination number, charge state, and pageranks
int myanalysis(double res[10], int snap, int idAl, int ntot, char Labels[Size][2], double Coords[Size][3], double clustersize, int graphAlO, double cutoffAlO, int graphOO, double cutoffOO, int graphOH, double cutoffOH, double dampingf, int Ntot, int sindex[Size], int getforce, double target, double kbias, double fkxyz[Size][3]);
// the fermi function
double myfermi(double x, double cutoff, double factor);
// here is the pagerank analysis based on all the atoms in the system
int wholeboxanalysis(int num, int ntot,char label[Size][2],double xyz[Size][3], int graphAlO, double cutoffAlO, int graphOO, double cutoffOO, int graphOH, double cutoffOH, double damping);
// this is a function to understand the type of interaction between node-m and node-n
int mygetype(int idm, int idn, int graphid[Size], char Labels[Size][2]);  // this is a function to understand the type of interaction between node-m and node-n

int main(int argc, char *argv[])
{

	// input file and the parameters

	if(argc != 2)
	{
		printf("Usage %s xyz-file\n\n",argv[0]);
		exit(-1);
	}

	char filename[FLN];
	sprintf(filename,"%s",argv[1]);
	FILE *fip,*fop,*fpp,*ffp;
	if( (fip=fopen(filename,"r"))==NULL )
	{
		printf("Error: cannot find %s\n",filename);
		exit(-1);
	}
	fop=fopen("Cluster.xyz","w");
	fpp=fopen("Output","w");


	int graphAlO, graphOO, graphOH;  // indicator whether to construct graph for Al..O, O..O and O..H
	graphAlO=graphOO=graphOH=1;    // default to construct all the graphs
	double cutoffAlO, cutoffOO, cutoffOH;  // the cutoffs, at the 0.5 weight
	cutoffAlO=cutoffOO=cutoffOH=0.0;
	double clustersize;  // the distance to Al center, serve as the cluster size for constructing graphs
	double damping=0.85;  // the damping factor for pagerank calculation
	int wholebox=0;
	int getforce;   // 2017.09.19, whether to calculate force or not
	double target, kbias; // the hamonic bias, U = kbias * (x - target)^2

	printf("Enter local Al cluster size\n");
	scanf("%lf",&clustersize);

	printf("To create Al..O graph?\n Enter 1(yes) or 0 (no)\n");
	scanf("%d",&graphAlO);
	if(graphAlO == 1)
	{
		printf(" graph Al..O, enter the cutoff\n");
		scanf("%lf",&cutoffAlO);
	}

	printf("To create O..O graph?\n Enter 1(yes) or 0 (no)\n");
	scanf("%d",&graphOO);
	if(graphOO == 1)
	{
		printf(" graph O..O, enter the cutoff\n");
		scanf("%lf",&cutoffOO);
	}

	printf("To create O..H graph?\n Enter 1(yes) or 0 (no)\n");
	scanf("%d",&graphOH);
	if(graphOH == 1)
	{
		printf(" graph O..H, enter the cutoff\n");
		scanf("%lf",&cutoffOH);
	}

	printf("Enter damping factor for pagerank calculation, (0 ~ 1)\n");
	scanf("%lf",&damping);
	if(damping<0 || damping>1)
	{
		printf("Wrong range of the damping factor, it should be a number from 0 to 1\n\n");
		exit(-1);
	}

	printf("To calculate PageRank of whole system?\n Enter 1(yes) or 0(no)\n");  // 2017.03.13, added pagerank analysis for the whole coordinates
	scanf("%d",&wholebox);

	target = 0.0; kbias = 0.0;
	printf("To calculate force due to a bias?\n Enter 1(yes) or 0(no)\n");   // 2017.09.19, add force calculation
	scanf("%d",&getforce);
	if(getforce == 1)
	{
		printf(" enter the target value:\n");
		scanf("%lf",&target);
		printf(" enter the bias constant:\n");
		scanf("%lf",&kbias);
	}

	printf("Summary: Boxsize %lf %lf %lf, Array size %d \n",Boxx,Boxy,Boxz,Size);
	printf("Summary: cluster size to Al center is %f\n",clustersize);
	printf("Summary: Graph-Al-O %d, cutoff %lf\n",graphAlO,cutoffAlO);
	printf("Summary: Graph-O-O %d, cutoff %lf\n",graphOO,cutoffOO);
	printf("Summary: Graph-O-H %d, cutoff %lf\n",graphOH,cutoffOH);
	printf("Summary: PageRank damping factor %lf\n",damping);
	printf("Summary: PageRank calculation of whole box %d\n",wholebox);
	printf("\n\n");

	// reading the xyz files, it is a loop over all snapshots
	char buffer[FLN];
	char *token;
	int ntot,i,j,k,num; // the number of atoms in a snapshot

	char label[Size][2];
	double xyz[Size][3];
	char slabel[Size][2];
	double sxyz[Size][3];
	int nn,cn,nAl,nO,nH,test;
	double sx,sy,sz,cx,cy,cz,jx,jy,jz,temp,Hx,Hy,Hz,tempOH;
	double fcn; // coordination number
	int id[Size];    // ids of the oxygens and hydrogens that are within clustersize of Al species
	double results[10];  // results of coordination number, charge state, pagerank of center, pagerank of oxygens,
	double fxyz[Size][3],fkxyz[Size][3];   // this is the force on each atom, 2017.09.18
	int sindex[Size];  // this is an array indicating the atom index of the selected cluster for calculating force

	num=1;
	while( fgets(buffer,sizeof(buffer),fip) != NULL)
	{
//		printf(" snap %d\n", num);
		token=strtok(buffer," ");
		ntot = atoi(token);    // the total number of atoms in the system
//		printf("%d\n",ntot);
		if(ntot >= Size)
		{
			printf("Error in main: too many atoms, increase 'Size'\n");
			break;   // this will break the while loop
		}

		fgets(buffer,sizeof(buffer),fip); // skip the second line in xyz format

		for(i=1;i<=ntot;i++)   // atom id goes from 1 to ntot
		{
			fgets(buffer,sizeof(buffer),fip);
			token = strtok(buffer," \n");
//printf(" %d %s %d\n",i,token,strlen(token));
			label[i][0] = token[0];    // read the atom label
			if( strlen(token)>1 )
				label[i][1] = token[1];
			else
				label[i][1] = ' ';
			j=0;
			while(token != NULL && j<3)   // read the xyz coordinates
			{
				token = strtok(NULL," \n");
				xyz[i][j] = atof(token);
				j++;
			}
		}

		for(i=0;i<=ntot;i++)
			sindex[i] = 0; // initialize this matrix for every frame

		nn=0; nAl=0;
		for(i=1;i<=ntot;i++)  // pick Al
		{
			if(label[i][0]=='A' && label[i][1]=='l')
			{
				nn++; nAl++;
				slabel[nn][0]='A'; slabel[nn][1]='l';
				sxyz[nn][0]=xyz[i][0]; sxyz[nn][1]=xyz[i][1]; sxyz[nn][2]=xyz[i][2];
				sindex[nn] = i;  // 2017.09.18, keep record of the original atom index
			}
		}

		for(i=0;i<Size;i++) 
			id[i]=0;
		for(k=1;k<=nAl;k++)    // for each Al center, pick the oxygens and hydrogens within a given range
		{
			cx=sxyz[k][0]; cy=sxyz[k][1]; cz=sxyz[k][2];
			nO=0; nH=0; cn=0; fcn=0.0;
			for(i=1;i<=ntot;i++) 
			{
				if(label[i][0]=='O' || label[i][0]=='o')   // loop over all oxygens and pick those with clustersize of Al center
				{
					sx=xyz[i][0]; sy=xyz[i][1]; sz=xyz[i][2];
					sx=pbcshift(sx,cx,'x'); sy=pbcshift(sy,cy,'y'); sz=pbcshift(sz,cz,'z');
					temp=sqrt( pow(sx-cx,2) + pow(sy-cy,2) + pow(sz-cz,2) );
//printf(" %d %d %f\n",k,i,temp);
					if(temp < Coordcut)  // coordination number using Al..O cutoff of 2.3
					{
						cn++;
						nO++;			
					}
					fcn = fcn + 1.0/(1.0+exp(Factor*(temp-Coordcut)/Coordcut));  // coordination number using a continuous fermi function at 2.3, from RDF

					if(temp < clustersize && id[i]==0) // pick oxygen if it is within clustersize of Al center
					{	
						nn++; 
						id[i]=1;   // this is to avoid duplication of atoms if it is with clustersize of two Al centers
						slabel[nn][0]='O'; slabel[nn][1]=' '; // oxygen label
						sxyz[nn][0]=sx; sxyz[nn][1]=sy; sxyz[nn][2]=sz;
						sindex[nn] = i;  // 2017.09.18, keep record of the original atom index, for Oxygen atoms

						for(j=1;j<=ntot;j++)
						{
							if(label[j][0]=='H' || label[j][0]=='h') // loop over all hydrogens and pick those near the oxygen atoms
							{
								Hx=xyz[j][0]; Hy=xyz[j][1]; Hz=xyz[j][2];
								Hx=pbcshift(Hx,sx,'x'); Hy=pbcshift(Hy,sy,'y'); Hz=pbcshift(Hz,sz,'z');
								tempOH=sqrt( pow(Hx-sx,2) + pow(Hy-sy,2) + pow(Hz-sz,2) );  // distance of hydrogen to oxygen
								Hx=pbcshift(Hx,cx,'x'); Hy=pbcshift(Hy,cy,'y'); Hz=pbcshift(Hz,cz,'z');

								if(temp < Coordcut && tempOH < CovOHcut)
									nH++;   // nH is number of H atoms that is bonded to the oxygens that are coordinated to Al center

								if(tempOH < CovOHcut && id[j]==0)   // this is to pick H atoms that is bonded to oxygens that are within clustersize of Al center
								{
									nn++;
									id[j]=1;
									slabel[nn][0]='H'; slabel[nn][1]=' '; // hydrogen label
									sxyz[nn][0]=Hx; sxyz[nn][1]=Hy; sxyz[nn][2]=Hz;
									sindex[nn] = j;   // 2017.09.18, keep record of the original atom index for Hydrogen atoms
								}	
							}
						}
					}  // end of pick oxygens within clustersize of Al center
				}
			} // end of loop over i

	//		printf(" snap %d CN %d %f CS %d nAl %d nO %d nH %d\n",num,cn,fcn,nAl*3-nO*2+nH,nAl,nO,nH);
		} // end of loop over k for Al centers

		// finally, print the cooridnates of the cluster;
		fprintf(fop,"    %d \n\n",nn);
		for(i=1;i<=nn;i++)
			fprintf(fop," %c%c  %f  %f  %f  \n",slabel[i][0],slabel[i][1],sxyz[i][0],sxyz[i][1],sxyz[i][2]);

		for(i=0;i<=ntot;i++)     // 2017.09.18, initialize the force on each atom
			fxyz[i][0] = fxyz[i][1] = fxyz[i][2] = 0.0;

		// here is the analysis based on the selected Al clusters, calculating the coordination number; the charge state; and the pagerank
		for(k=1;k<=nAl;k++)
		{
			for(i=0;i<=ntot;i++)
				fkxyz[i][0] = fkxyz[i][0] = fkxyz[i][0] = 0.0; // this array is the force due to k-th Al
			results[0]=results[1]=results[2]=results[3]=results[4]=results[5]=0.0;
			test=myanalysis(results,num,k,nn,slabel,sxyz,clustersize,graphAlO,cutoffAlO,graphOO,cutoffOO,graphOH,cutoffOH,damping,ntot,sindex,getforce,target,kbias,fkxyz);
			if(test != 0)   // there is error
				break;

			for(i=1; i<=ntot; i++)
			{
				fxyz[i][0] += fkxyz[i][0];   // the array of fxyz keeps the total force due to all the Al clusters
				fxyz[i][1] += fkxyz[i][1];
				fxyz[i][2] += fkxyz[i][2];
			}

			printf(" snap %d Al-id %d : CN %f CS %f PR(Al) %e PR(O)-avg %e PR(O)-diff %e fCN %f \n",num,k,results[0],results[1],results[2],results[3],results[4],results[5]);
			fprintf(fpp," snap %d Al-id %d : CN %f CS %f PR(Al) %e PR(O)-avg %e PR(O)-diff %e fCN %f \n",num,k,results[0],results[1],results[2],results[3],results[4],results[5]);
		}
		if(test != 0)
			break;

		sprintf(filename,"snap%d.total-force",num);     // here is to print the force at each frame
		ffp = fopen(filename,"w");
		for(i=1;i<=ntot;i++)
			fprintf(ffp, "%e  %e  %e  \n",fxyz[i][0], fxyz[i][1], fxyz[i][2]);

		// here is the analysis for the whole box
		if(wholebox == 1)
		{
			test=wholeboxanalysis(num,ntot,label,xyz,graphAlO,cutoffAlO,graphOO,cutoffOO,graphOH,cutoffOH,damping);
			if(test != 0)
				break;	
		}

		num = num + 1;

	}   // this is the end of while loop

	fclose(fip);
	fclose(fop);
	fclose(fpp);
	fclose(ffp);

	return(0);
}

// calculate the coordination number, charge state, and pageranks
int myanalysis(double res[10], int snap, int idAl, int ntot, char Labels[Size][2], double Coords[Size][3], double clustersize, int graphAlO, double cutoffAlO, int graphOO, double cutoffOO, int graphOH, double cutoffOH, double dampingf, int Ntot, int sindex[Size], int getforce, double target, double kbias, double fkxyz[Size][3])
{

	int iO,iH,numO,numH,iH3;
	double Alx,Aly,Alz,Ox,Oy,Oz,O2x,O2y,O2z,Hx,Hy,Hz,distAlO,distOO,distOH,tempdist;
	double fCN;    // 2017.06.27, smooth coordination number using fermi function
	int idCoord[Size]; // the id of oxygens that is in the coordination shell of Al center
	for(iO=0;iO<Size;iO++)
		idCoord[iO]=0;

	if(!(Labels[idAl][0]=='A' && Labels[idAl][1]=='l'))
	{
		printf("Error in 'getcncs' the id %d is not an Al atom\n",idAl);
		return(-1);
	}

//printf("   %d %d : CN & CS\n",snap,idAl);
	// first to calculate the coordination number and the charges state of the Al center
	Alx=Coords[idAl][0]; Aly=Coords[idAl][1]; Alz=Coords[idAl][2];
	numO=0;numH=0; fCN=0.0;
	for(iO=1;iO<=ntot;iO++)
	{
		if(Labels[iO][0]=='O' || Labels[iO][0]=='o')
		{
			Ox=Coords[iO][0]; Oy=Coords[iO][1]; Oz=Coords[iO][2];
			Ox=pbcshift(Ox,Alx,'x'); Oy=pbcshift(Oy,Aly,'y'); Oz=pbcshift(Oz,Alz,'z');
			distAlO=sqrt( pow(Ox-Alx,2)+pow(Oy-Aly,2)+pow(Oz-Alz,2) );
			fCN += myfermi(distAlO,Coordcut,CNfactor);    // continuous coordination number
			if(distAlO <= Coordcut)
			{
				numO++;
				idCoord[iO]=1;
				iH3=0;  // this is to find whether there is a H3O+, that is an oxygen with 3 H atoms within CovOHcut
				for(iH=1;iH<=ntot;iH++)
				{
					if(Labels[iH][0]=='H' || Labels[iH][0]=='h')
					{
						Hx=Coords[iH][0]; Hy=Coords[iH][1]; Hz=Coords[iH][2];
						Hx=pbcshift(Hx,Ox,'x'); Hy=pbcshift(Hy,Oy,'y'); Hz=pbcshift(Hz,Oz,'z');
						distOH=sqrt( pow(Hx-Ox,2)+pow(Hy-Oy,2)+pow(Hz-Oz,2) );
						if(distOH <= CovOHcut)
						{
							iH3++;
							numH++;
						}
					}
				}
				if(iH3>=3)
					printf(" Warning in 'getcncs': there is H3O+ from oxygen %f %f %f, that is %f away from Al %d\n",Ox,Oy,Oz,distAlO,idAl);
			}
		}
	}
	res[0]= numO+0.0;  // this is coordination number
	res[1]= 3-2*numO+numH+0.0; // this is the charge state
	res[5]= fCN;

//printf("   %d %d : graphs\n",snap,idAl);
	// second, to calculate the pageranks for the chosen graphs, there are two steps, one is to create the graph, the other is to calculate pagerank
	FILE *fgraph,*fPageRank; // the output file for the graphs and the PageRank vector
	char filename[FLN];
	int graphid[Size]; // the list of the atom id in the graph
	int numatom,mi,mj,id1,id2,tempnum;
	double **Amatrix,**Smatrix,**Gmatrix; // the adjacency matrix for the graph, and google matrix for pagerank calculation
//	double Alx,Aly,Alz,Ox,Oy,Oz,O2x,O2y,O2z,Hx,Hy,Hz,distAlO,distOO,distOH,tempdist;   // coordinates has been declared at the beginning of this function
	double *Lweight,totweight;  // the sum of the weight of connections with each node
	double **distmatrix, **xmatrix, **ymatrix, **zmatrix;   // this is the distance between nodes

	sprintf(filename,"snap%d.Al%d.Graph-%.2f-%.2f-%.2f",snap,idAl,cutoffAlO,cutoffOO,cutoffOH);   // graph file
	fgraph = fopen(filename,"w");

	// this part is to determine the number of atoms in the graph
	numatom=1;      // numatom goes from 1 to the actual number
	graphid[numatom]=idAl; // the first id is the Al center
	Alx=Coords[idAl][0]; Aly=Coords[idAl][1]; Alz=Coords[idAl][2];

	if(graphAlO==1 || graphOO==1)  // add oxygen ids if the oxygen graphs are needed
	{
		for(mi=1;mi<=ntot;mi++)
		{
			if(Labels[mi][0]=='O' || Labels[mi][0]=='o')   // choose oxygens
			{
				Ox=Coords[mi][0]; Oy=Coords[mi][1]; Oz=Coords[mi][2];
				Ox=pbcshift(Ox,Alx,'x'); Oy=pbcshift(Oy,Aly,'y'); Oz=pbcshift(Oz,Alz,'z');
				distAlO=sqrt( pow(Ox-Alx,2)+pow(Oy-Aly,2)+pow(Oz-Alz,2) );
				if(distAlO < clustersize) // only include semi-graph, oxygens and hydrogens within a range
				{
					numatom++;
					graphid[numatom]=mi;
				}
			}
		}   // put oxygen ids first and then hydrogen ids
	}
	tempnum = numatom;  // therefore, the oxygen index goes from 2 to tempnum

	if(graphOH==1)   // add hydrogen ids if the hydrogen graphs are needed
	{
		for(mi=1;mi<=ntot;mi++)
		{
			if(Labels[mi][0]=='H' || Labels[mi][0]=='h')   // examine all the hydrogen atoms
			{
				Hx=Coords[mi][0]; Hy=Coords[mi][1]; Hz=Coords[mi][2];
				for(mj=2;mj<=tempnum;mj++)   // find the hydrogens from all the oxygens
				{
					Ox=Coords[graphid[mj]][0]; Oy=Coords[graphid[mj]][1]; Oz=Coords[graphid[mj]][2];
					Ox=pbcshift(Ox,Hx,'x'); Oy=pbcshift(Oy,Hy,'y'); Oz=pbcshift(Oz,Hz,'z');
					distOH=sqrt( pow(Hx-Ox,2)+pow(Hy-Oy,2)+pow(Hz-Oz,2) );
					if(distOH < CovOHcut)   // includes the hydrogen atoms, that belongs to the chosen oxygens
					{
						numatom++;
						graphid[numatom]=mi;
						break;     // once the hydrogen is include in the graph, skip from looping oxygen, this is to prevent duplication of this hydrogen when it is close to more than one oxygen (although it usually will not happen)
					}
				}
			}
		}
	}

	
	Amatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		Amatrix[mi] = (double *)calloc(numatom,sizeof(double));
	Smatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		Smatrix[mi] = (double *)calloc(numatom,sizeof(double));
	Gmatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		Gmatrix[mi] = (double *)calloc(numatom,sizeof(double));
	Lweight = (double *)calloc(numatom,sizeof(double));

	distmatrix = (double **)malloc(numatom*sizeof(double *));    // distmatrix, xmatrix, ymatrix, zmatrix, record the actual distances between two nodes
	for(mi=0;mi<numatom;mi++)
		distmatrix[mi] = (double *)calloc(numatom,sizeof(double));
	xmatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		xmatrix[mi] = (double *)calloc(numatom,sizeof(double));
	ymatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		ymatrix[mi] = (double *)calloc(numatom,sizeof(double));
	zmatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		zmatrix[mi] = (double *)calloc(numatom,sizeof(double));

	for(mi=0;mi<numatom;mi++)    // this is to build the Amatrix
	{
		for(mj=mi+1;mj<numatom;mj++)
		{
			id1=graphid[mi+1]; id2=graphid[mj+1];   // the atom index goes from 1 to numatom
			tempdist=0.0;
			if(graphAlO==1)   // construct Al..O graph
			{
				if((Labels[id1][0]=='A' && Labels[id1][1]=='l') && (Labels[id2][0]=='O' || Labels[id2][0]=='o'))
				{
					Alx=Coords[id1][0]; Aly=Coords[id1][1]; Alz=Coords[id1][2];
					Ox=Coords[id2][0]; Oy=Coords[id2][1]; Oz=Coords[id2][2];
					Ox=pbcshift(Ox,Alx,'x'); Oy=pbcshift(Oy,Aly,'y'); Oz=pbcshift(Oz,Alz,'z');
					tempdist = distAlO = sqrt( pow(Ox-Alx,2)+pow(Oy-Aly,2)+pow(Oz-Alz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distAlO,cutoffAlO,Factor);   // undirect graph is used, treat as double-sized graph
					distmatrix[mi][mj] = distmatrix[mj][mi] = distAlO;
					xmatrix[mi][mj] = Alx - Ox;  xmatrix[mj][mi] = Ox - Alx;
					ymatrix[mi][mj] = Aly - Oy;  ymatrix[mj][mi] = Oy - Aly;
					zmatrix[mi][mj] = Alz - Oz;  zmatrix[mj][mi] = Oz - Alz;
				}
				if((Labels[id2][0]=='A' && Labels[id2][1]=='l') && (Labels[id1][0]=='O' || Labels[id1][0]=='o'))
				{
					Alx=Coords[id2][0]; Aly=Coords[id2][1]; Alz=Coords[id2][2];
					Ox=Coords[id1][0]; Oy=Coords[id1][1]; Oz=Coords[id1][2];
					Ox=pbcshift(Ox,Alx,'x'); Oy=pbcshift(Oy,Aly,'y'); Oz=pbcshift(Oz,Alz,'z');
					tempdist = distAlO = sqrt( pow(Ox-Alx,2)+pow(Oy-Aly,2)+pow(Oz-Alz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distAlO,cutoffAlO,Factor);   // undirect graph is used, treat as double-sized graph
					distmatrix[mi][mj] = distmatrix[mj][mi] = distAlO;
					xmatrix[mi][mj] = Ox - Alx; xmatrix[mi][mj] = Alx - Ox;
					ymatrix[mi][mj] = Oy - Aly; ymatrix[mi][mj] = Aly - Oy;
					zmatrix[mi][mj] = Oz - Alz; zmatrix[mi][mj] = Alz - Oz;
				}
			}
			if(graphOO==1)  // construct O..O graph
			{
				if((Labels[id1][0]=='O' || Labels[id1][0]=='o') && (Labels[id2][0]=='O' || Labels[id2][0]=='o'))
				{
					Ox=Coords[id1][0]; Oy=Coords[id1][1]; Oz=Coords[id1][2];   
					O2x=Coords[id2][0]; O2y=Coords[id2][1]; O2z=Coords[id2][2];
					O2x=pbcshift(O2x,Ox,'x'); O2y=pbcshift(O2y,Oy,'y'); O2z=pbcshift(O2z,Oz,'z');
					tempdist = distOO = sqrt( pow(O2x-Ox,2)+pow(O2y-Oy,2)+pow(O2z-Oz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distOO,cutoffOO,Factor);   // undirect graph is used, treat as double-sized graph
					distmatrix[mi][mj] = distmatrix[mj][mi] = distOO;
					xmatrix[mi][mj] = Ox - O2x; xmatrix[mj][mi] = O2x - Ox;
					ymatrix[mi][mj] = Oy - O2y; ymatrix[mj][mi] = O2y - Oy;
					zmatrix[mi][mj] = Oz - O2z; zmatrix[mj][mi] = O2z - Oz;
				}
			}
			if(graphOH==1)   // construct O..H graph
			{
				if((Labels[id1][0]=='O' || Labels[id1][0]=='o') && (Labels[id2][0]=='H' || Labels[id2][0]=='h'))
				{
					Ox=Coords[id1][0]; Oy=Coords[id1][1]; Oz=Coords[id1][2];   
					Hx=Coords[id2][0]; Hy=Coords[id2][1]; Hz=Coords[id2][2];
					Hx=pbcshift(Hx,Ox,'x'); Hy=pbcshift(Hy,Oy,'y'); Hz=pbcshift(Hz,Oz,'z');
					tempdist = distOH = sqrt( pow(Hx-Ox,2)+pow(Hy-Oy,2)+pow(Hz-Oz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distOH,cutoffOH,Factor);   // undirect graph is used, treat as double-sized graph
					distmatrix[mi][mj] = distmatrix[mj][mi] = distOH;
					xmatrix[mi][mj] = Ox - Hx; xmatrix[mj][mi] = Hx - Ox;
					ymatrix[mi][mj] = Oy - Hy; ymatrix[mj][mi] = Hy - Oy;
					zmatrix[mi][mj] = Oz - Hz; zmatrix[mj][mi] = Hz - Oz;
				}
				if((Labels[id2][0]=='O' || Labels[id2][0]=='o') && (Labels[id1][0]=='H' || Labels[id1][0]=='h'))
				{
					Ox=Coords[id2][0]; Oy=Coords[id2][1]; Oz=Coords[id2][2];   
					Hx=Coords[id1][0]; Hy=Coords[id1][1]; Hz=Coords[id1][2];
					Hx=pbcshift(Hx,Ox,'x'); Hy=pbcshift(Hy,Oy,'y'); Hz=pbcshift(Hz,Oz,'z');
					tempdist = distOH = sqrt( pow(Hx-Ox,2)+pow(Hy-Oy,2)+pow(Hz-Oz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distOH,cutoffOH,Factor);   // undirect graph is used, treat as double-sized graph
					distmatrix[mi][mj] = distmatrix[mj][mi] = distOH;
					xmatrix[mi][mj] = Hx - Ox; xmatrix[mj][mi] = Ox - Hx;
					ymatrix[mi][mj] = Hy - Oy; ymatrix[mj][mi] = Oy - Hy;
					zmatrix[mi][mj] = Hz - Oz; zmatrix[mj][mi] = Oz - Hz;
				}

			}
			if( tempdist>0.0 )
				fprintf(fgraph,"%d %d %d %d %c%c %c%c %f %f \n",mi+1,mj+1,id1,id2,Labels[id1][0],Labels[id1][1],Labels[id2][0],Labels[id2][1],tempdist,Amatrix[mi][mj]);  // print this connection
		}
	}

	fclose(fgraph);

//	for(mi=0; mi<numatom; mi++)
//	{
//		for(mj=0; mj<numatom; mj++)
//			printf("%f ",Amatrix[mi][mj]);
//		printf("\n");
//	}
//	printf("\n");


	for(mi=0;mi<numatom;mi++)   // this is to construct the total weight of each node in graph
	{
		totweight=0.0;
		for(mj=0;mj<numatom;mj++)
			totweight += Amatrix[mi][mj];
		Lweight[mi] = totweight;
	}

	for(mi=0;mi<numatom;mi++) // this is to construct the matrix S, which is a column-stochastic matrix
	{
		for(mj=0;mj<numatom;mj++)
		{
			if(Lweight[mj] > 0.0)      // this is to make S is column-stochastic
				Smatrix[mi][mj] = Amatrix[mi][mj] / Lweight[mj];
			else        // if there is a dangling node, then use this normalized term
				Smatrix[mi][mj] = 1.0/numatom;
		}
	}

	for(mi=0;mi<numatom;mi++)
		for(mj=0;mj<numatom;mj++)
			Gmatrix[mi][mj] = dampingf * Smatrix[mi][mj] + (1.0-dampingf) * 1.0 / numatom;   // this is the definition of Google matrix

//	for(mi=0; mi<numatom; mi++)
//	{
//		for(mj=0; mj<numatom; mj++)
//			printf("%f ",Gmatrix[mi][mj]);
//		printf("\n");
//	}
//	printf("\n");

/*	// 2017.09.12, this part has been changed, I'm using matrix-S to construct matrix-G
	for(mi=0;mi<numatom;mi++)
	{
		for(mj=0;mj<numatom;mj++)
		{
			if(Lweight[mj] > 0.0)
				Gmatrix[mi][mj] = dampingf*Amatrix[mi][mj]/Lweight[mj]+(1.0-dampingf)/numatom;
			else
				Gmatrix[mi][mj] = 1.0/numatom;   // if there is a sink, i.e. a node without any connection, it will probably not happen
		}
	}
*/


//printf("   %d %d : pageranks\n",snap,idAl);
	// this part is to calculate the pagerank by solving the eigenvalue problem of Gmatrix, I use power method, i.e. multiple many times of Gmatrix until it converges
	double epsilon,sum,avg,diff;
	int iter,ii,jj;
	double *PRi,*PRj;
	double **GJmatrix;
	double mtemp;

	PRi = (double *)calloc(numatom,sizeof(double));
//	PRj = (double *)calloc(numatom,sizeof(double));

	for(ii=0;ii<numatom;ii++)
		PRi[ii]=1.0/numatom;

	// 2017.09.14, here I want to use matrix inverse to solve the eigenvalue problem of Google matrix, i.e. to find the PageRank vector
	GJmatrix = (double **)malloc(numatom * sizeof(double *));    // I'm using Gauss-Jordan method, so this is N * (2N) matrix, of the original matrix and an indentify matrix
	for(ii=0;ii<numatom;ii++)
		GJmatrix[ii] = (double *)calloc(numatom*2, sizeof(double));

	for(ii=0;ii<numatom;ii++)
	{
		for(jj=0;jj<numatom;jj++)
		{
			if(jj == ii)
				GJmatrix[ii][jj] = 1.0 - dampingf * Smatrix[ii][jj];    // in this problem, I want to get the inverse of matrix: (I-alpha*S)
			else
				GJmatrix[ii][jj] = 0.0 - dampingf * Smatrix[ii][jj];
		}

		for(jj=numatom;jj<numatom*2;jj++)    // this is the Identity matrix
		{
			if(jj == ii+numatom)
				GJmatrix[ii][jj] = 1.0;
			else
				GJmatrix[ii][jj] = 0.0;
		}
	}

//	for(ii=0;ii<numatom;ii++)
//	{
//		for(jj=0;jj<numatom*2;jj++)
//			printf("%f ",GJmatrix[ii][jj]);
//		printf("\n");
//	}
//	printf("\n");

	// Below, I'm using Gauss-Jordan method to calculate the matrix inverse
	// In this method, GJmatrix (A | I) is reducted using row operations to be (I | A-1)
	for(ii=0; ii<numatom; ii++)
	{
		if(GJmatrix[ii][ii] == 0)  // if this row, ii, starts with 0, find another row_j, and swap them
		{
			for(jj=ii+1; jj<numatom; jj++) // searching the following rows
				if(GJmatrix[jj][ii] != 0)
					break;    // find a row that starts with non-zero
			if(jj >= numatom)
			{
				printf("Error, matrix is not invertible\n");
				break;
			}
			else     // swap two rows
			{
				for(iter=0; iter<numatom*2; iter++)
				{
					mtemp = GJmatrix[ii][iter];
					GJmatrix[ii][iter] = GJmatrix[jj][iter];
					GJmatrix[jj][iter] = mtemp;
				}
			}
		}   // after this loop, it is ensured that the first element of row ii, is not zero

		mtemp = GJmatrix[ii][ii];  // it is not zero now
		for(iter=0; iter<numatom*2; iter++)
			GJmatrix[ii][iter] = GJmatrix[ii][iter] / mtemp;  // make the first element to 1

		for(jj=0;jj<numatom;jj++) // to make all other rows zero in column ii, expect for row ii
		{
			if(jj == ii)
				continue;
			else
			{
				mtemp = GJmatrix[jj][ii];
				for(iter=0; iter<numatom*2; iter++)
					GJmatrix[jj][iter] = GJmatrix[jj][iter] - GJmatrix[ii][iter] * mtemp;
			}
		}

	}  // this is the end of Gauss-Jordan method, now the GJmatrix is [ I | (I-alpha*S)-1 ], the inverse is obtained

//	for(ii=0;ii<numatom;ii++)
//	{
//		for(jj=0;jj<numatom*2;jj++)
//			printf("%f ",GJmatrix[ii][jj]);
//		printf("\n");
//	}
//	printf("\n");

	for(ii=0; ii<numatom; ii++)   // PageRank vector can be calculated from: PR = (I-alpha*S)-1 * (1.0-alpha)1.0/N e
	{
		PRi[ii] = 0.0;
		for(jj=numatom; jj<numatom*2; jj++)
			PRi[ii] += GJmatrix[ii][jj] * (1.0-dampingf) * 1.0/numatom;
	}


/*	// 2017.09.12, this part uses Power-method to solve the eigenvalue problem of Google matrix, changed it to using matrix inverse
	epsilon=1.0; iter=1;
	while(epsilon > Convtol)
	{
		for(ii=0;ii<numatom;ii++)   // update new PRj vector by multiply Gmatrix to the left of old PRi vector
		{
			PRj[ii]=0.0;
			for(jj=0;jj<numatom;jj++)
				PRj[ii] += Gmatrix[ii][jj]*PRi[jj];      // see Gmatrix definition before
		}
		sum=0.0;  // to normalize the PageRank vector so that its sum equals to 1
		for(ii=0;ii<numatom;ii++)
			sum += PRj[ii];
		for(ii=0;ii<numatom;ii++)
			PRj[ii] = PRj[ii]/sum;

		epsilon=0.0; // check the convergence criterion
		for(ii=0;ii<numatom;ii++)
		{
			epsilon += fabs(PRi[ii]-PRj[ii]);
//printf("   %f %f\n",PRi[ii],PRj[ii]);
		}
		epsilon = epsilon*numatom;

//printf("Iteration %d epsilon %e\n",iter,epsilon);
		iter +=1;
		if(iter > 1000)    // if it does not converge
			break;

		for(ii=0;ii<numatom;ii++)
			PRi[ii]=PRj[ii];
	}

	if(iter > 1000)
	{
		printf("Error: pagerank calculation not converge for snap %d Al-id %d\n",snap,idAl);
		printf("  Iteration %d epsilon %e\n",iter,epsilon);
		for(ii=0;ii<numatom;ii++)
			printf("   %d %d %c%c %f %f \n",ii+1,graphid[ii+1],Labels[graphid[ii+1]][0],Labels[graphid[ii+1]][1],PRi[ii],PRj[ii]);
		return(-1);
	}
//printf("   %d %d : print\n",snap,idAl);
	*/	// 2017.09.12, change this Power-method part to matrix inverse


	sprintf(filename,"snap%d.Al%d.PR-%.3f",snap,idAl,dampingf);    // pagerank file
	fPageRank = fopen(filename,"w");

	for(ii=0;ii<numatom;ii++)
		fprintf(fPageRank,"%d %d %c%c %f \n",ii+1,graphid[ii+1],Labels[graphid[ii+1]][0],Labels[graphid[ii+1]][1],PRi[ii]);	

	fclose(fPageRank);


	// 2017.09.19, this part, I'm calculating the derivative of PageRank vector with respect to atom coordinate
	double **dPRdx, **dPRdy, **dPRdz;
	double ***dSmatrixdx, ***dSmatrixdy, ***dSmatrixdz;
	int dk,di,dm,dn;
	int typemn, typekn;
	double cutmn, cutkn, sumkx, sumky, sumkz;
	double *Nforcex, *Nforcey, *Nforcez;  // this is the array of force for each atom
	FILE *ffp; 

	Nforcex = (double *) calloc(Ntot+1, sizeof(double));
	Nforcey = (double *) calloc(Ntot+1, sizeof(double));
	Nforcez = (double *) calloc(Ntot+1, sizeof(double));

	dPRdx = (double **) malloc( numatom * sizeof(double *) );
	for(dk = 0; dk<numatom; dk++)
		dPRdx[dk] = (double *) calloc(numatom, sizeof(double));
	dPRdy = (double **) malloc( numatom * sizeof(double *) );
	for(dk = 0; dk<numatom; dk++)
		dPRdy[dk] = (double *) calloc(numatom, sizeof(double));
	dPRdz = (double **) malloc( numatom * sizeof(double *) );
	for(dk = 0; dk<numatom; dk++)
		dPRdz[dk] = (double *) calloc(numatom, sizeof(double));

	dSmatrixdx = (double ***) malloc( numatom * sizeof(double **) );
	for(dk=0; dk<numatom; dk++)
	{
		dSmatrixdx[dk] = (double **) malloc( numatom * sizeof(double *) );
		for(di=0; di<numatom; di++)
			dSmatrixdx[dk][di] = (double *) calloc(numatom, sizeof(double) );
	}
	dSmatrixdy = (double ***) malloc( numatom * sizeof(double **) );
	for(dk=0; dk<numatom; dk++)
	{
		dSmatrixdy[dk] = (double **) malloc( numatom * sizeof(double *) );
		for(di=0; di<numatom; di++)
			dSmatrixdy[dk][di] = (double *) calloc(numatom, sizeof(double) );
	}
	dSmatrixdz = (double ***) malloc( numatom * sizeof(double **) );
	for(dk=0; dk<numatom; dk++)
	{
		dSmatrixdz[dk] = (double **) malloc( numatom * sizeof(double *) );
		for(di=0; di<numatom; di++)
			dSmatrixdz[dk][di] = (double *) calloc(numatom, sizeof(double) );
	}

if(getforce == 1)  // 
{
	sprintf(filename,"snap%d.Al%d.force",snap,idAl);
	ffp = fopen(filename,"w");


	// first calculate dSmatrix / dx ...
	for(dm=0; dm<numatom; dm++)
	{
		for(dn=0; dn<numatom; dn++)
		{
			if(Lweight[dn] == 0.0)  // if the sum of column-n is zero in Amatrix
			{
				for(di=0; di<numatom; di++)
				{
					dSmatrixdx[dm][dn][di] = 0.0;
					dSmatrixdy[dm][dn][di] = 0.0;
					dSmatrixdz[dm][dn][di] = 0.0;
				}
			}
			else
			{
				typemn = mygetype(dm,dn,graphid,Labels);  // this is a function to understand the type of interaction between node-m and node-n
				if(typemn == 1 && graphAlO == 1)
					cutmn = cutoffAlO;
				if(typemn == 2 && graphOO == 1)
					cutmn = cutoffOO;
				if(typemn == 3 && graphOH == 1)
					cutmn = cutoffOH;

				for(di=0; di<numatom; di++)
				{
					if(dm == di)   // this if statement and the below if(dn == di) are the first term in dSmatrix / dx
					{
						if(Amatrix[dm][dn] != 0.0)
							dSmatrixdx[dm][dn][di] += 1.0 / Lweight[dn] * (-1.0 * Amatrix[dm][dn] * Amatrix[dm][dn] * (1.0/Amatrix[dm][dn] - 1.0) * Factor/cutmn * xmatrix[dm][dn]/distmatrix[dm][dn]);
						else
							dSmatrixdx[dm][dn][di] += 0.0;

						if(Amatrix[dm][dn] != 0.0)
							dSmatrixdy[dm][dn][di] += 1.0 / Lweight[dn] * (-1.0 * Amatrix[dm][dn] * Amatrix[dm][dn] * (1.0/Amatrix[dm][dn] - 1.0) * Factor/cutmn * ymatrix[dm][dn]/distmatrix[dm][dn]);
						else
							dSmatrixdy[dm][dn][di] += 0.0;

						if(Amatrix[dm][dn] != 0.0)
							dSmatrixdz[dm][dn][di] += 1.0 / Lweight[dn] * (-1.0 * Amatrix[dm][dn] * Amatrix[dm][dn] * (1.0/Amatrix[dm][dn] - 1.0) * Factor/cutmn * zmatrix[dm][dn]/distmatrix[dm][dn]);
						else
							dSmatrixdz[dm][dn][di] += 0.0;
						
					}
					if(dn == di)   // this along with the previous if, are the first term in dSmatrix / dx
					{
						if(Amatrix[dm][dn] != 0.0)
							dSmatrixdx[dm][dn][di] += -1.0 / Lweight[dn] * (-1.0 * Amatrix[dm][dn] * Amatrix[dm][dn] * (1.0/Amatrix[dm][dn] - 1.0) * Factor/cutmn * xmatrix[dm][dn]/distmatrix[dm][dn]);
						else
							dSmatrixdx[dm][dn][di] += 0.0;

						if(Amatrix[dm][dn] != 0.0)
							dSmatrixdy[dm][dn][di] += -1.0 / Lweight[dn] * (-1.0 * Amatrix[dm][dn] * Amatrix[dm][dn] * (1.0/Amatrix[dm][dn] - 1.0) * Factor/cutmn * ymatrix[dm][dn]/distmatrix[dm][dn]);
						else
							dSmatrixdy[dm][dn][di] += 0.0;

						if(Amatrix[dm][dn] != 0.0)
							dSmatrixdz[dm][dn][di] += -1.0 / Lweight[dn] * (-1.0 * Amatrix[dm][dn] * Amatrix[dm][dn] * (1.0/Amatrix[dm][dn] - 1.0) * Factor/cutmn * zmatrix[dm][dn]/distmatrix[dm][dn]);
						else
							dSmatrixdz[dm][dn][di] += 0.0;
					}

					// below I'm getting the second term in dSmatrix / dx
					sumkx = sumky = sumkz = 0.0;
					for(dk=0; dk<numatom; dk++)
					{
						typekn = mygetype(dk,dn,graphid,Labels);  // this is a function to understand the type of interaction between node-m and node-n
						if(typekn == 1 && graphAlO == 1)
							cutkn = cutoffAlO;
						if(typekn == 2 && graphOO == 1)
							cutkn = cutoffOO;
						if(typekn == 3 && graphOH == 1)
							cutkn = cutoffOH;

						if(dk == di)
						{
							if(Amatrix[dk][dn] != 0.0)
								sumkx += -1.0 * Amatrix[dk][dn] * Amatrix[dk][dn] * (1.0/Amatrix[dk][dn] - 1.0) * Factor/cutkn * xmatrix[dk][dn]/distmatrix[dk][dn];
							else
								sumkx += 0.0;

							if(Amatrix[dk][dn] != 0.0)
								sumky += -1.0 * Amatrix[dk][dn] * Amatrix[dk][dn] * (1.0/Amatrix[dk][dn] - 1.0) * Factor/cutkn * ymatrix[dk][dn]/distmatrix[dk][dn];
							else
								sumkx += 0.0;
						
							if(Amatrix[dk][dn] != 0.0)
								sumkz += -1.0 * Amatrix[dk][dn] * Amatrix[dk][dn] * (1.0/Amatrix[dk][dn] - 1.0) * Factor/cutkn * zmatrix[dk][dn]/distmatrix[dk][dn];
							else
								sumkz += 0.0;

						}
						if(dn == di)
						{
							if(Amatrix[dk][dn] != 0.0)
								sumkx += 1.0 * Amatrix[dk][dn] * Amatrix[dk][dn] * (1.0/Amatrix[dk][dn] - 1.0) * Factor/cutkn * xmatrix[dk][dn]/distmatrix[dk][dn];
							else
								sumkx += 0.0;

							if(Amatrix[dk][dn] != 0.0)
								sumky += 1.0 * Amatrix[dk][dn] * Amatrix[dk][dn] * (1.0/Amatrix[dk][dn] - 1.0) * Factor/cutkn * ymatrix[dk][dn]/distmatrix[dk][dn];
							else
								sumkx += 0.0;
						
							if(Amatrix[dk][dn] != 0.0)
								sumkz += 1.0 * Amatrix[dk][dn] * Amatrix[dk][dn] * (1.0/Amatrix[dk][dn] - 1.0) * Factor/cutkn * zmatrix[dk][dn]/distmatrix[dk][dn];
							else
								sumkz += 0.0;
						}

					}  // end of the for loop of dk

					dSmatrixdx[dm][dn][di] += -1.0 * Amatrix[dm][dn]/(Lweight[dn] * Lweight[dn]) * sumkx;
					dSmatrixdy[dm][dn][di] += -1.0 * Amatrix[dm][dn]/(Lweight[dn] * Lweight[dn]) * sumky;
					dSmatrixdz[dm][dn][di] += -1.0 * Amatrix[dm][dn]/(Lweight[dn] * Lweight[dn]) * sumkz;
 
				}  // end of the loop for di

			}  // end of if if(Lweight[dn] == 0.0) else

		}  // end of loop for dn
	}  // end of loop for dm

	// now start to get the derivative of pagerank with respect to coordinate, i.e. dPRdx, dPRdy, dPRdz
	for(dk=0; dk<numatom; dk++)   // dk and di are variable, dm and dn are dummy number in summation
	{
		for(di=0; di<numatom; di++)
		{
			for(dm=0; dm<numatom; dm++)
			{
				for(dn=0; dn<numatom; dn++)
				{
					dPRdx[dk][di] += GJmatrix[dk][dm+numatom] * dampingf * dSmatrixdx[dm][dn][di] * PRi[di];
					dPRdy[dk][di] += GJmatrix[dk][dm+numatom] * dampingf * dSmatrixdy[dm][dn][di] * PRi[di];
					dPRdz[dk][di] += GJmatrix[dk][dm+numatom] * dampingf * dSmatrixdz[dm][dn][di] * PRi[di];
				}
			}
		}
	}

	// then to calculate the force, if the bias is on PageRank of Al (index 0 in PR[] vector)
	for(di=0; di<numatom; di++)   // di is the node id in graph, PR[] vector, ...
	{
		dk = graphid[di+1];  // this is the node id in the selected cluster
		dm = sindex[dk];   // this is the original atom id

		Nforcex[dm] = kbias * (PRi[0] - target) * dPRdx[0][di];
		Nforcey[dm] = kbias * (PRi[0] - target) * dPRdy[0][di];
		Nforcez[dm] = kbias * (PRi[0] - target) * dPRdz[0][di];

		fprintf(ffp,"%d %d %d %e %e %e\n",di+1,dk,dm,Nforcex[dm],Nforcey[dm],Nforcez[dm]);
	}
	fclose(ffp);

	for(di=1; di<=Ntot; di++) // here I get the force on each atom, if the atom is not in Al cluster, its force is zero
	{
		fkxyz[di][0] = Nforcex[di];
		fkxyz[di][1] = Nforcey[di];
		fkxyz[di][2] = Nforcez[di];
	}


} // end of if (getforce == 1)


	avg=0.0; jj=0;
	for(ii=2;ii<=tempnum;ii++)
	{
		avg += PRi[ii];
		jj += 1;
	}
	if(jj>0)
		avg = avg/jj;
	else
		avg = -1;
/*
	for(ii=1;ii<=ntot;ii++)   // this loop is to calculate the PageRank of oxygens that directly coordinate Al center
	{
		if(idCoord[ii]==1)    // this array indicates whether it is the oxygen that coordinates the Al center
		{
			avg += PRi[ii];
			jj += 1;
		}
	}
	if(jj>0)
		avg = avg/jj;
	else
		avg = -1;
*/


	diff=0.0; jj=0;
	for(ii=2;ii<=tempnum;ii++)
	{
		diff += fabs(PRi[ii]-avg);
		jj += 1;
	}
	if(jj>0)
		diff = diff/jj;
	else
		diff = -1;
/*
	for(ii=1;ii<=ntot;ii++)  // this loop is to calculate the mean absolute deviation of PageRank of oxygens that directly coordinate Al center
	{
		if(idCoord[ii]==1)    // this array indicates whether it is the oxygen that coordinates the Al center
		{
			diff += fabs(PRi[ii]-avg);
			jj += 1;
		}
	}
	if(jj>0)
		diff = diff/jj;
	else
		diff = -1;
*/

	res[2]=PRi[0];   // the PageRank value of the Al center
	res[3]=avg;   // the average PageRank value of the oxygens that directly coordinate Al center
	res[4]=diff;  // mean absolute deviation of PageRank of oxygens that directly coordinate AL center


	// here I free all the dynamically allocated arrays, matrixes
	for(dm=0; dm<numatom; dm++)
	{
		free(Amatrix[dm]);
		free(Smatrix[dm]);
		free(Gmatrix[dm]);
		free(distmatrix[dm]);
		free(xmatrix[dm]);
		free(ymatrix[dm]);
		free(zmatrix[dm]);
		free(GJmatrix[dm]);
		free(dPRdx[dm]);
		free(dPRdy[dm]);
		free(dPRdz[dm]);
	}
	free(Amatrix);
	free(Smatrix);
	free(Gmatrix);
	free(distmatrix);
	free(xmatrix);
	free(ymatrix);
	free(zmatrix);
	free(GJmatrix);
	free(dPRdx);
	free(dPRdy);
	free(dPRdz);

	free(Lweight);
	free(PRi);
	free(Nforcex);
	free(Nforcey);
	free(Nforcez);

	for(dm=0; dm<numatom; dm++)
	{
		for(dn=0; dn<numatom; dn++)
		{
			free(dSmatrixdx[dm][dn]);
			free(dSmatrixdy[dm][dn]);
			free(dSmatrixdz[dm][dn]);
		}
		free(dSmatrixdx[dm]);
		free(dSmatrixdy[dm]);
		free(dSmatrixdz[dm]);
	}
	free(dSmatrixdx);
	free(dSmatrixdy);
	free(dSmatrixdz);


	return(0);

} 

// here is the pagerank analysis based on all the atoms in the system, copied and modified 2017.03.13
int wholeboxanalysis(int num, int ntot,char label[Size][2],double xyz[Size][3], int graphAlO, double cutoffAlO, int graphOO, double cutoffOO, int graphOH, double cutoffOH, double damping)
{
	// to calculate the pageranks for the whole graphs, there are two steps, one is to create the graph, the other is to calculate pagerank
	FILE *fgraph,*fPageRank; // the output file for the graphs and the PageRank vector
	char filename[FLN];
	int graphid[Size]; // the list of the atom id in the graph
	int numatom,mi,mj,id1,id2;
	double **Amatrix,**Smatrix,**Gmatrix; // the adjacency matrix for the graph, and google matrix for pagerank calculation
	double Alx,Aly,Alz,Ox,Oy,Oz,O2x,O2y,O2z,Hx,Hy,Hz,distAlO,distOO,distOH,tempdist;  
	double *Lweight,totweight;  // the sum of the weight of connections with each node

	sprintf(filename,"snap%d.whole.Graph-%.2f-%.2f-%.2f",num,cutoffAlO,cutoffOO,cutoffOH);   // graph file
	fgraph = fopen(filename,"w");

	// there is no need to re-order the atoms, and keep the order as in the original xyz file
	numatom=0;
	for(mi=1;mi<=ntot;mi++)
	{
		numatom++;
		graphid[numatom]=mi;
	}

//printf("    %d\n",numatom);
	
	Amatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		Amatrix[mi] = (double *)calloc(numatom,sizeof(double));
	Smatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		Smatrix[mi] = (double *)calloc(numatom,sizeof(double));
	Gmatrix = (double **)malloc(numatom*sizeof(double *));
	for(mi=0;mi<numatom;mi++)
		Gmatrix[mi] = (double *)calloc(numatom,sizeof(double));
	Lweight = (double *)calloc(numatom,sizeof(double));

	for(mi=0;mi<numatom;mi++)    // this is to build the Amatrix
	{
		for(mj=mi+1;mj<numatom;mj++)
		{
			id1=graphid[mi+1]; id2=graphid[mj+1];   // the atom index goes from 1 to numatom
			tempdist=0.0;
			if(graphAlO==1)   // construct Al..O graph
			{
				if((label[id1][0]=='A' && label[id1][1]=='l') && (label[id2][0]=='O' || label[id2][0]=='o'))
				{
					Alx=xyz[id1][0]; Aly=xyz[id1][1]; Alz=xyz[id1][2];
					Ox=xyz[id2][0]; Oy=xyz[id2][1]; Oz=xyz[id2][2];
					Ox=pbcshift(Ox,Alx,'x'); Oy=pbcshift(Oy,Aly,'y'); Oz=pbcshift(Oz,Alz,'z');
					tempdist = distAlO = sqrt( pow(Ox-Alx,2)+pow(Oy-Aly,2)+pow(Oz-Alz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distAlO,cutoffAlO,Factor);   // undirect graph is used, treat as double-sized graph
				}
				if((label[id2][0]=='A' && label[id2][1]=='l') && (label[id1][0]=='O' || label[id1][0]=='o'))
				{
					Alx=xyz[id2][0]; Aly=xyz[id2][1]; Alz=xyz[id2][2];
					Ox=xyz[id1][0]; Oy=xyz[id1][1]; Oz=xyz[id1][2];
					Ox=pbcshift(Ox,Alx,'x'); Oy=pbcshift(Oy,Aly,'y'); Oz=pbcshift(Oz,Alz,'z');
					tempdist = distAlO = sqrt( pow(Ox-Alx,2)+pow(Oy-Aly,2)+pow(Oz-Alz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distAlO,cutoffAlO,Factor);   // undirect graph is used, treat as double-sized graph
				}
			}
			if(graphOO==1)  // construct O..O graph
			{
				if((label[id1][0]=='O' || label[id1][0]=='o') && (label[id2][0]=='O' || label[id2][0]=='o'))
				{
					Ox=xyz[id1][0]; Oy=xyz[id1][1]; Oz=xyz[id1][2];   
					O2x=xyz[id2][0]; O2y=xyz[id2][1]; O2z=xyz[id2][2];
					O2x=pbcshift(O2x,Ox,'x'); O2y=pbcshift(O2y,Oy,'y'); O2z=pbcshift(O2z,Oz,'z');
					tempdist = distOO = sqrt( pow(O2x-Ox,2)+pow(O2y-Oy,2)+pow(O2z-Oz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distOO,cutoffOO,Factor);   // undirect graph is used, treat as double-sized graph
				}
			}
			if(graphOH==1)   // construct O..H graph
			{
				if((label[id1][0]=='O' || label[id1][0]=='o') && (label[id2][0]=='H' || label[id2][0]=='h'))
				{
					Ox=xyz[id1][0]; Oy=xyz[id1][1]; Oz=xyz[id1][2];   
					Hx=xyz[id2][0]; Hy=xyz[id2][1]; Hz=xyz[id2][2];
					Hx=pbcshift(Hx,Ox,'x'); Hy=pbcshift(Hy,Oy,'y'); Hz=pbcshift(Hz,Oz,'z');
					tempdist = distOH = sqrt( pow(Hx-Ox,2)+pow(Hy-Oy,2)+pow(Hz-Oz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distOH,cutoffOH,Factor);   // undirect graph is used, treat as double-sized graph
				}
				if((label[id2][0]=='O' || label[id2][0]=='o') && (label[id1][0]=='H' || label[id1][0]=='h'))
				{
					Ox=xyz[id2][0]; Oy=xyz[id2][1]; Oz=xyz[id2][2];   
					Hx=xyz[id1][0]; Hy=xyz[id1][1]; Hz=xyz[id1][2];
					Hx=pbcshift(Hx,Ox,'x'); Hy=pbcshift(Hy,Oy,'y'); Hz=pbcshift(Hz,Oz,'z');
					tempdist = distOH = sqrt( pow(Hx-Ox,2)+pow(Hy-Oy,2)+pow(Hz-Oz,2) );
					Amatrix[mi][mj] = Amatrix[mj][mi] = myfermi(distOH,cutoffOH,Factor);   // undirect graph is used, treat as double-sized graph
				}
			}
			if( tempdist>0.0 )
				fprintf(fgraph,"%d %d %d %d %c%c %c%c %f %f \n",mi+1,mj+1,id1,id2,label[id1][0],label[id1][1],label[id2][0],label[id2][1],tempdist,Amatrix[mi][mj]);  // print this connection
		}
	}

	fclose(fgraph);


	for(mi=0;mi<numatom;mi++)   // this is to construct the total weight of each node in graph
	{
		totweight=0.0;
		for(mj=0;mj<numatom;mj++)
			totweight += Amatrix[mi][mj];
		Lweight[mi] = totweight;
	}

	for(mi=0;mi<numatom;mi++) // this is to construct the matrix S, which is a column-stochastic matrix
	{
		for(mj=0;mj<numatom;mj++)
		{
			if(Lweight[mj] > 0.0)      // this is to make S is column-stochastic
				Smatrix[mi][mj] = Amatrix[mi][mj] / Lweight[mj];
			else        // if there is a dangling node, then use this normalized term
				Smatrix[mi][mj] = 1.0/numatom;
		}
	}

	for(mi=0;mi<numatom;mi++)
		for(mj=0;mj<numatom;mj++)
			Gmatrix[mi][mj] = damping * Smatrix[mi][mj] + (1.0-damping) * 1.0 / numatom;   // this is the definition of Google matrix

/*
	for(mi=0;mi<numatom;mi++)
	{
		for(mj=0;mj<numatom;mj++)
		{
			if(Lweight[mj] > 0.0)
				Gmatrix[mi][mj] = damping*Amatrix[mi][mj]/Lweight[mj]+(1.0-damping)/numatom;
			else
				Gmatrix[mi][mj] = 1.0/numatom;   // if there is a sink, i.e. a node without any connection, it will probably not happen
		}
	}
*/

	// this part is to calculate the pagerank by solving the eigenvalue problem of Gmatrix, I use power method, i.e. multiple many times of Gmatrix until it converges
	double epsilon,sum,avg,diff;
	int iter,ii,jj;
	double *PRi,*PRj;
	double **GJmatrix;
	double mtemp;

	PRi = (double *)calloc(numatom,sizeof(double));
//	PRj = (double *)calloc(numatom,sizeof(double));

	for(ii=0;ii<numatom;ii++)
		PRi[ii]=1.0/numatom;

	// 2017.09.14, here I want to use matrix inverse to solve the eigenvalue problem of Google matrix, i.e. to find the PageRank vector
	GJmatrix = (double **)malloc(numatom * sizeof(double *));    // I'm using Gauss-Jordan method, so this is N * (2N) matrix, of the original matrix and an indentify matrix
	for(ii=0;ii<numatom;ii++)
		GJmatrix[ii] = (double *)calloc(numatom*2, sizeof(double));

	for(ii=0;ii<numatom;ii++)
	{
		for(jj=0;jj<numatom;jj++)
		{
			if(jj == ii)
				GJmatrix[ii][jj] = 1.0 - damping * Smatrix[ii][jj];    // in this problem, I want to get the inverse of matrix: (I-alpha*S)
			else
				GJmatrix[ii][jj] = 0.0 - damping * Smatrix[ii][jj];
		}

		for(jj=numatom;jj<numatom*2;jj++)    // this is the Identity matrix
		{
			if(jj == ii+numatom)
				GJmatrix[ii][jj] = 1.0;
			else
				GJmatrix[ii][jj] = 0.0;
		}
	}

//	for(ii=0;ii<numatom;ii++)
//	{
//		for(jj=0;jj<numatom*2;jj++)
//			printf("%f ",GJmatrix[ii][jj]);
//		printf("\n");
//	}
//	printf("\n");

	// Below, I'm using Gauss-Jordan method to calculate the matrix inverse
	// In this method, GJmatrix (A | I) is reducted using row operations to be (I | A-1)
	for(ii=0; ii<numatom; ii++)
	{
	//	printf("calculating inverse: row %d\n",ii);
		if(GJmatrix[ii][ii] == 0)  // if this row, ii, starts with 0, find another row_j, and swap them
		{
			for(jj=ii+1; jj<numatom; jj++) // searching the following rows
				if(GJmatrix[jj][ii] != 0)
					break;    // find a row that starts with non-zero
			if(jj >= numatom)
			{
				printf("Error, matrix is not invertible\n");
				break;
			}
			else     // swap two rows
			{
				for(iter=0; iter<numatom*2; iter++)
				{
					mtemp = GJmatrix[ii][iter];
					GJmatrix[ii][iter] = GJmatrix[jj][iter];
					GJmatrix[jj][iter] = mtemp;
				}
			}
		}   // after this loop, it is ensured that the first element of row ii, is not zero

		mtemp = GJmatrix[ii][ii];  // it is not zero now
		for(iter=0; iter<numatom*2; iter++)
			GJmatrix[ii][iter] = GJmatrix[ii][iter] / mtemp;  // make the first element to 1

		for(jj=0;jj<numatom;jj++) // to make all other rows zero in column ii, expect for row ii
		{
			if(jj == ii)
				continue;
			else
			{
				mtemp = GJmatrix[jj][ii];
				for(iter=0; iter<numatom*2; iter++)
					GJmatrix[jj][iter] = GJmatrix[jj][iter] - GJmatrix[ii][iter] * mtemp;
			}
		}

	}  // this is the end of Gauss-Jordan method, now the GJmatrix is [ I | (I-alpha*S)-1 ], the inverse is obtained

//	for(ii=0;ii<numatom;ii++)
//	{
//		for(jj=0;jj<numatom*2;jj++)
//			printf("%f ",GJmatrix[ii][jj]);
//		printf("\n");
//	}
//	printf("\n");

	for(ii=0; ii<numatom; ii++)   // PageRank vector can be calculated from: PR = (I-alpha*S)-1 * (1.0-alpha)1.0/N e
	{
		PRi[ii] = 0.0;
		for(jj=numatom; jj<numatom*2; jj++)
			PRi[ii] += GJmatrix[ii][jj] * (1.0-damping) * 1.0/numatom;
	}




	sprintf(filename,"snap%d.whole.PR-%.3f",num,damping);    // pagerank file
	fPageRank = fopen(filename,"w");

	for(ii=0;ii<numatom;ii++)
		fprintf(fPageRank,"%d %d %c%c %e \n",ii+1,graphid[ii+1],label[graphid[ii+1]][0],label[graphid[ii+1]][1],PRi[ii]);	

	fclose(fPageRank);

	
	// here I free all the dynamically allocated arrays, matrixes
	for(ii=0; ii<numatom; ii++)
	{
		free(Amatrix[ii]);
		free(Smatrix[ii]);
		free(Gmatrix[ii]);
		free(GJmatrix[ii]);
	}
	free(Amatrix);
	free(Smatrix);
	free(Gmatrix);
	free(GJmatrix);

	free(Lweight);
	free(PRi);


	return(0);
}


// considering pbc, shift the coordinates
double pbcshift(double inp, double center, char dim)
{
	double out;

	if(dim == 'x' || dim == 'X')
	{
		if(inp - center > Boxx/2)
			out=inp-Boxx;
		else if(inp - center < -1*Boxx/2)
			out=inp+Boxx;
		else
			out=inp;
	}
	else if (dim == 'y' || dim == 'Y')
	{
		if(inp - center > Boxy/2)
			out=inp-Boxy;
		else if(inp - center < -1*Boxy/2)
			out=inp+Boxy;
		else
			out=inp;
	}
	else if(dim == 'z' || dim == 'Z')
	{
		if(inp - center > Boxz/2)
			out=inp-Boxz;
		else if(inp - center < -1*Boxz/2)
			out=inp+Boxz;
		else
			out=inp;
		
	}

	return(out);
}


// the fermi function
double myfermi(double x, double cutoff, double factor)  // the fermi function
{
	double result;

	result = (1.0)/(1.0 + exp(factor*(x-cutoff)/cutoff));

	return result;
}

int mygetype(int idm, int idn, int graphid[Size], char Labels[Size][2])  // this is a function to understand the type of interaction between node-m and node-n
{
	int output=0;
	int sm,sn;

	sm = graphid[idm+1];  // the graphid starts from 1 to numatoms
	sn = graphid[idn+1];

	if(Labels[sm][0] == 'A' && Labels[sm][1] == 'l')
	{
		if(Labels[sn][0] == 'O' || Labels[sn][0] == 'o')
			output = 1;     // return 1 if it is Al..O interaction
		else if(Labels[sn][0] == 'H' || Labels[sn][0] == 'h')
			output = 4;     // return 4 if it is Al..H interaction, which we do not use
		else if(Labels[sn][0] == 'A' && Labels[sn][1] == 'l')
			output = 7;    // return 7 if it is Al..Al interaction, which we do not use
	}
	else if(Labels[sm][0] == 'O' || Labels[sm][0] == 'o')
	{
		if(Labels[sn][0] == 'A' && Labels[sn][1] == 'l')
			output = 1;  // return 1 if it is O..Al interaction
		else if(Labels[sn][0] == 'H' || Labels[sn][0] == 'h')
			output = 3;  // return 3 if it is O..H interaction
		else if(Labels[sn][0] == 'O' || Labels[sn][0] == 'o')
			output = 2;  // return 2 if it is O..O interaction
	}
	else if(Labels[sm][0] == 'H' || Labels[sm][0] == 'h')
	{
		if(Labels[sn][0] == 'A' && Labels[sn][1] == 'l')
			output = 4;  // return 4 if it is Al..H interaction, which we do not use
		else if(Labels[sn][0] == 'O' || Labels[sn][0] == 'o')
			output = 3;  // return 2 if it is O..H interaction
		else if(Labels[sn][0] == 'H' || Labels[sn][0] == 'h')
			output = 9;  // return 9 if it is H..H interaction, which we do not use
	} 

	return(output);
}



