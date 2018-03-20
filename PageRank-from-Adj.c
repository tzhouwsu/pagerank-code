// this is an PageRank calculation code by reading an adjacency matrix from a file
// this is "for Kerry" from Aurora, 2018.01.23

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define FLN 1000

// this function is to calculate matrix inverse using Gauss-Jordan method
int get_matrix_inverse(int num, float **M, float **InvM);


int main(int argc, char *argv[])
{
	int Nnodes;
	char filename[FLN];
	FILE *fip;
	float damping_factor;
	float temp;
	float **Adj; // this is the adjacency matrix
	float **Smatrix, **I_minus_alpha_S, **inverse_of_I_minus_alpha_S;
	float *Scolumn, *PR_final, *PR_initial;
	int i,j,syn;
	char buffer[FLN];
	char *token;

	if(argc != 4)
	{
		printf("Error, use this as: \n  %s size-of-matrix file-of-matrix damping-factor\n\n",argv[0]);
		exit(-1);
	}
	
	Nnodes = atoi(argv[1]);
	sprintf(filename,"%s",argv[2]);
	damping_factor = atof(argv[3]);

	Adj = (float **)malloc(Nnodes*sizeof(float *));  // allocate the adjacency matrix
	for(i=0;i<Nnodes;i++)
		Adj[i] = (float *)calloc(Nnodes, sizeof(float));

	Smatrix = (float **)malloc(Nnodes*sizeof(float *));  // allocate the matrix S
	for(i=0;i<Nnodes;i++)
		Smatrix[i] = (float *)calloc(Nnodes, sizeof(float));

	I_minus_alpha_S = (float **)malloc(Nnodes*sizeof(float *));  // allocate matrix I - alpha*S
	for(i=0;i<Nnodes;i++)
		I_minus_alpha_S[i] = (float *)calloc(Nnodes, sizeof(float));

	inverse_of_I_minus_alpha_S = (float **)malloc(Nnodes*sizeof(float *));  // allocate matrix (I - alpha*S)^-1
	for(i=0;i<Nnodes;i++)
		inverse_of_I_minus_alpha_S[i] = (float *)calloc(Nnodes, sizeof(float));

	Scolumn = (float *)calloc(Nnodes, sizeof(float));
	PR_initial = (float *)calloc(Nnodes, sizeof(float));
	PR_final = (float *)calloc(Nnodes, sizeof(float));

	if((fip=fopen(filename,"r"))==NULL)
	{
		printf("Error, cannot find file %s\n",filename);
	}
	else  // reading the Adjacency matrix
	{
		rewind(fip);

		i=0;
		while(fgets(buffer,sizeof(buffer),fip) != NULL)
		{
			j=0;
			token = strtok(buffer," \n");
			Adj[i][j] = atof(token);
			for(j=1;j<Nnodes;j++)
			{
				token = strtok(NULL ," \n");
				Adj[i][j] = atof(token);
			}

			i++;
		}
		fclose(fip);

//		for(i=0;i<Nnodes;i++)
//		{
//			for(j=0;j<Nnodes;j++)
//				printf("%f ",Adj[i][j]);
//			printf("\n");
//		}

		// construct matrix S, I-alpha*S, and (I-alpha*S)^-1
		for(i=0;i<Nnodes;i++)
		{
			Scolumn[i] = 0.0;
			for(j=0;j<Nnodes;j++)
			{
				Scolumn[i] += Adj[j][i];
			}

			for(j=0;j<Nnodes;j++)   // Smatrix is a column-stochastic matrix
			{
				if(Scolumn[i] == 0.0)
					Smatrix[j][i] = 1.0/Nnodes;
				else
					Smatrix[j][i] = Adj[j][i]/Scolumn[i];
			}
		}

		for(i=0;i<Nnodes;i++)
		{
			for(j=0;j<Nnodes;j++)
			{
				if(j==i)
					I_minus_alpha_S[i][j] = 1.0 - damping_factor*Smatrix[i][j];
				else
					I_minus_alpha_S[i][j] = -damping_factor*Smatrix[i][j];
			}
		}

		syn = get_matrix_inverse(Nnodes,I_minus_alpha_S,inverse_of_I_minus_alpha_S);
		if(syn == 0) // when there is no error in getting matrix inverse
		{
			for(i=0;i<Nnodes;i++)   // initial guess of PageRank
				PR_initial[i] = 1.0/Nnodes;

			for(i=0;i<Nnodes;i++)
			{
				PR_final[i] = 0.0;
				for(j=0;j<Nnodes;j++)
				{
					PR_final[i] += (1.0-damping_factor) * inverse_of_I_minus_alpha_S[i][j] * PR_initial[j];
				}
			}
		}
		else
			printf("Error, in calculating matrix inverse \n");
		
	
	}

	// print the PageRank result
	for(i=0;i<Nnodes;i++)
		printf("%f \n",PR_final[i]);

	// free the Arrays;
	for(i=0;i<Nnodes;i++)
	{
		free(Adj[i]);
		free(Smatrix[i]);
		free(I_minus_alpha_S[i]);
		free(inverse_of_I_minus_alpha_S[i]);
	}
	free(Adj);
	free(Smatrix);
	free(I_minus_alpha_S);
	free(inverse_of_I_minus_alpha_S);

	return(0);
}


// this function is to calculate matrix inverse using Gauss-Jordan method
int get_matrix_inverse(int num, float **M, float **InvM)
{
	float **GJ_matrix;
	int m,n,irow,jrow, icolumn, result;
	float temp;

	result = 0;
	GJ_matrix = (float **)malloc(num*sizeof(float *));
	for(m=0;m<num;m++)
		GJ_matrix[m] = (float *)calloc(2*num, sizeof(float)); 

	for(m=0;m<num;m++)
	{
		for(n=0;n<num;n++)
		{
			GJ_matrix[m][n] = M[m][n];
//			printf("%f ",M[m][n]);
		}
//		printf("\n");
		GJ_matrix[m][num+m] = 1.0;
	}

	for(irow=0; irow<num; irow++)
	{
		if(GJ_matrix[irow][irow] == 0)  // make sure the row_i start with non-zero
		{
			for(jrow=irow+1; jrow<num; jrow++)
				if(GJ_matrix[jrow][irow] != 0)
					break;
			if(jrow >= num)
			{
				printf("error, matrix invertible\n");
				result = -1;
				break;
			}
			else
			{
				for(icolumn=0; icolumn<num*2; icolumn++)   // switch two rows
				{
					temp = GJ_matrix[irow][icolumn];
					GJ_matrix[irow][icolumn] = GJ_matrix[jrow][icolumn];
					GJ_matrix[jrow][icolumn] = temp;
				}
			}
		}

		temp = GJ_matrix[irow][irow];
		for(icolumn=0; icolumn<num*2; icolumn++)
			GJ_matrix[irow][icolumn] = GJ_matrix[irow][icolumn] / temp;

		for(jrow=0; jrow<num; jrow++)
		{
			if(jrow ==  irow)
				continue;
			else
			{
				temp = GJ_matrix[jrow][irow];
				for(icolumn=0; icolumn<num*2; icolumn++)
					GJ_matrix[jrow][icolumn] = GJ_matrix[jrow][icolumn] - GJ_matrix[irow][icolumn] * temp;
			}
		}
	}

	for(m=0;m<num;m++)
	{
		for(n=0;n<num;n++)
		{
			InvM[m][n] = GJ_matrix[m][n+num];
//			printf("%f ",GJ_matrix[m][n+num]);
		}
//		printf("\n");
	}

	return(result);
}





