#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N 1024
#define PI acos(-1.0)
#define b 0.3

void init(double *array, double *array_n, double *arrayp, double *arrayp_n);

void main(int argc, char **argv){
	int rank, source, destination, in_number, out_number;
	
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	double dt = 0.001;
	int T = 100*N;
	int t = 0;
	int i;
	double *x;
	double *x_n;
	double *xp;
	double *xp_n;
	double *ai;
	double *ai1;
	x = (double *)malloc(N*sizeof(double));
	x_n = (double *)malloc(N*sizeof(double));
	xp = (double *)malloc(N*sizeof(double));
	xp_n = (double *)malloc(N*sizeof(double));
	ai = (double *)malloc(N*sizeof(double));
	ai1 = (double *)malloc(N*sizeof(double));
	ai[0], ai[N-1], ai1[0], ai1[N-1] = 0.0;
	init(x, x_n, xp, xp_n);
	
	while(t<T)
	{
		for(i=1;i<N-1;i++)
		{
			ai[i] = (x[i-1]+x[i+1]-2*x[i])+b*(pow(x[i+1]-x[i],3.0)-pow(x[i]-x[i-1],3));
			x_n[i] = x[i]+xp[i]*dt+0.5*ai[i]*dt*dt;
			ai1[i] = (x_n[i-1]+x_n[i+1]-2*x_n[i])+b*(pow(x_n[i+1]-x_n[i],3.0)-pow(x_n[i]-x_n[i-1],3));
			xp_n[i] = xp[i]+0.5*(ai[i]+ai1[i])*dt;
		}
		
		for(i=0;i<N;i++)
		{
			x[i] = x_n[i];
			xp[i] = xp_n[i];
		}
		
		if(t % 100 == 0)
		{
			FILE *atpos = fopen("atpos.dat", "a");
			for(i=0;i<N;i++)
			{
				fprintf(atpos,"%f ",x[i]);
				//printf("%f\n",b);
			}
			fprintf(atpos,"\n");
			fclose(atpos);
		}
		t++;
	}
	MPI_Finalize();
}

void init(double *array, double *array_n, double *arrayp, double *arrayp_n){
	int i;
	for(i=0;i<N;i++)
	{
		array[i] = sin(2*PI*i/((double)(N-1)));
		array_n[i] = sin(2*PI*i/((double)(N-1)));
		arrayp[i] = 0.0;
		arrayp_n[i] = 0.0;
	}
	array[0], array[N-1] = 0.0;
}
