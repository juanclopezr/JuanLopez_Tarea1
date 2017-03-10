#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

#define N 64
#define PI acos(-1.0)
#define b 1.0

double freq(int mod);
double energy(int mod, double *pos, double *speed);
void *init(double *pos, double *v);
double acceleration(int i, double *pos);
double modenergy(int mod, double *pos);

int step = 1000;

int main(int argc, char **argv)
{
    double dt = 0.005, T = 5*pow(N, 2.2), t = 0, E1, E2, E3;
    int i, contador = 0, j = 0, n = (T/(dt*step)), I=0;
    double *x, *vn, *vhalf, w1, w2, w3, E1_, E2_, E3_;
    
    x = malloc(N*sizeof(double));
    vn = malloc(N*sizeof(double));
    vhalf = malloc(N*sizeof(double));
        
    init(x, vn);
    if (argc > 1)
    {
        omp_set_num_threads(atoi(argv[1]));
    }
    
    FILE *Ener;
    Ener = fopen("Energies.dat", "w");
    while(t<T)
    {
        #pragma omp parallel for
        for(i=1;i<N-1;i++)
        {
            vhalf[i] = vn[i] + 0.5*dt*acceleration(i, x);
		}
        #pragma omp parallel for
		for(i=1;i<N-1;i++)
		{
            x[i] += vhalf[i]*dt;
		}
        #pragma omp parallel for
		for(i=1;i<N-1;i++)
		{
            vn[i] = vhalf[i] + 0.5*dt*acceleration(i, x);
        }
        if(contador%n == 0)
        {
            E1 = energy(1, x, vn);
            E2 = energy(2, x, vn);
            E3 = energy(3, x, vn);
            fprintf(Ener,"%f %f %f %f\n", j*n*dt, E1, E2, E3);
            j += 1;
        }
        t += dt;
        contador += 1;
	}
    fclose(Ener);
    return 0;
}

void *init(double *pos, double *v)
{
    int i;
    for(i=0;i<N;i++)
    {
        pos[i] = sin(PI*i/(N-1.0));
        v[i] = 0.0;
    }
    pos[0] = pos[N-1] = 0.0;
}

double acceleration(int i, double *pos)
{
    return (pos[i-1] - 2*pos[i] + pos[i+1]) + b*(pow(pos[i+1] - pos[i], 2.0) - pow(pos[i] - pos[i-1], 2.0));
}

double freq(int mod)
{
    return 4.0*pow(sin(mod*PI/(2*N+2.0)), 2);
}

double energy(int mod, double *pos, double *speed)
{
    double E_p, E_s, w;
    E_p = modenergy(mod, pos);
    E_s = modenergy(mod, speed);
    w = freq(mod);
    return 0.5*(w*pow(E_p, 2) + pow(E_s, 2));
}

double modenergy(int mod, double *pos)
{	
	int i;
    double Ak = 0;
	for(i=1; i<N; i++)
	{
		Ak = Ak + pos[i]*sin(i*mod*PI/N);
    }
	return Ak*pow(2.0/N, 0.5);
}
