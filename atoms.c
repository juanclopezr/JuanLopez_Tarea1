#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"

#define N 64
#define PI acos(-1.0)
#define b 1.0

double *init(double *pos, double *v);
double acceleration(int i, double *pos);
double modenergy(int mod, double *pos);

int main(int argc, char **argv)
{
    int rank, source, destination, in_number, out_number;
    double dt = 5e-3, T = 5.0*pow(N, 2.2), t = 0, a, a_;
    int i, contador = 0, j = 0, n = (T/dt)/1000;
    double *x, *vn, *vhalf;
    x = (double *)malloc(N*sizeof(double));
    vn = (double *)malloc(N*sizeof(double));
    vhalf = (double *)malloc(N*sizeof(double));


    double dt_half = 0.5*dt, dt_squared = 0.5*pow(dt, 2);
    x, vn = init(x, vn);
    if (argc > 1)
    {
        omp_set_num_threads(atoi(argv[1]));
    }
    FILE *atpos, *Ener;
    atpos = fopen("atpos.dat", "w");
    Ener = fopen("Energies.dat", "w");

    while(t<T)
    {
        #pragma omp parallel for private(a_)
        for(i=1;i<N-1;i++)
        {
            a_ = acceleration(i, x);
            x[i] += vn[i]*dt + a_*dt_squared;
            vn[i] += dt_half*(acceleration(i, x)+a_);
        }
        if(contador%n == 0)
        {
            printf("%d\n", j);
            for(i=0; i<N; i++)
            {
                fprintf(atpos,"%f\t", x[i]);
            }
            fprintf(atpos,"\n");
            fprintf(Ener,"%f %f %f %f\n",t,modenergy(1,x),modenergy(2,x),modenergy(3,x));
            j += 1;
        }
        t += dt;
        contador ++;
	}
    fclose(atpos);
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
    double ac = pos[i-1] + pos[i+1] - 2*pos[i] + b*(pow(pos[i+1] - pos[i], 3.0) - pow(pos[i] - pos[i-1], 3.0));
    return ac;
}

double modenergy(int mod, double *pos)
{
	double wk = 4*pow(sin((double)mod*PI)/(2*(double)N+2),2.0);
	double Ak = 0.0;
	double Energy;
	int i;
	for(i=0;i<N;i++)
	{
		Ak += sqrt(2/((double)N+1))*pos[i]*sin((double)(i*mod)*PI/((double)N+1));
	}
	Energy = 0.5*(pow(Ak,2.0)+pow(wk*Ak,2.0));
	return Energy;
}
