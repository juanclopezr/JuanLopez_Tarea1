#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

#define N 64
#define PI acos(-1.0)
#define b 1.0

void derivative(double *pos, double *def, double h);
double freq(int mod);
double *init(double *pos, double *v);
double acceleration(int i, double *pos);
double modenergy(int mod, double *pos);

int step = 1000;

int main(int argc, char **argv)
{
    int rank, source, destination, in_number, out_number;
    double dt = 5e-3, T = 5.0*pow(N, 2.2), t = 0, a, a_, wk, D1, D2, D3;
    int i, contador = 0, j = 0, n = (T/dt)/step;
    double *x, *vn, *vhalf, E1[2], E2[2], E3[2];
    x = (double *)malloc(N*sizeof(double));
    vn = (double *)malloc(N*sizeof(double));
    vhalf = (double *)malloc(N*sizeof(double));
    double values[3];
    double dt_half = 0.5*dt, dt_squared = 0.5*pow(dt, 2);
    x, vn = init(x, vn);
    if (argc > 1)
    {
        omp_set_num_threads(atoi(argv[1]));
    }
    FILE *Ener;
    Ener = fopen("Energies.dat", "w");

    while(t<T)
    {
        #pragma omp parallel for private(a_)
        for(i=1;i<N-1;i++)
        {
            a_ = acceleration(i, x);
		}
		for(i=1;i<N-1;i++)
		{
            x[i] += vn[i]*dt + a_*dt_squared;
		}
		for(i=1;i<N-1;i++)
		{
            vn[i] += dt_half*(acceleration(i, x)+a_);
        }
	/*
        if(contador%n == 0)
        {
            printf("%d\n", j);
            E1[0] = modenergy(1, x);
            E2[0] = modenergy(2, x);
            E3[0] = modenergy(3, x);
            j += 1;
        }*/
        if (contador%n == 0)
        {
	    E1[0] = modenergy(1, vn);
	    E2[0] = modenergy(2, vn);
	    E3[0] = modenergy(3, vn);
            E1[1] = modenergy(1, x);
            E2[1] = modenergy(2, x);
            E3[1] = modenergy(3, x);/*
            derivative(E1, &D1, dt);
            derivative(E2, &D2, dt);
            derivative(E3, &D3, dt);*/
            values[0] = 0.5*(pow(E1[0], 2.0) + pow(freq(1)*E1[1], 2.0));
            values[1] = 0.5*(pow(E2[0], 2.0) + pow(freq(2)*E2[1], 2.0));
            values[2] = 0.5*(pow(E3[0], 2.0) + pow(freq(3)*E3[1], 2.0));
            fprintf(Ener,"%f %f %f %f\n", j*n*dt, values[0], values[1], values[2]);
	    j += 1;
        }
        t += dt;
        contador ++;
	}
    fclose(Ener);
    return 0;
}

double *init(double *pos, double *v){
    int i;
    for(i=0;i<N;i++)
    {
        pos[i] = sin(PI*i/((double)(N-1)));
        v[i] = 0;
    }
    pos[0], pos[N-1] = 0.0;
    return pos, v;
}

double acceleration(int i, double *pos)
{
    double ac = pos[i-1] + pos[i+1] - 2*pos[i] + b*(pow(pos[i+1] - pos[i], 2.0) - pow(pos[i] - pos[i-1], 2.0));
    return ac;
}

void derivative(double *pos, double *def, double h)
{
    *def = (pos[1] - pos[0])/h;
}

double freq(int mod)
{
    return 4.0*pow(sin((double)mod*PI)/(2*(double)N+2), 2.0);
}

double modenergy(int mod, double *pos)
{	
	double Ak = 0;
	int i;
	for(i=0;i<N;i++)
	{
		Ak += pos[i]*sin(((double)(i+1)*mod)*PI/((double)(N+1)));
    }
	return Ak*sqrt(2/((double)(N+1)));
}
