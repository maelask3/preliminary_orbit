#include "MatlabUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <string.h>
#include "rpoly.h"


double norm(double *v)
{
        return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}

double dot(double *v, double *t)
{
	return(v[0]*t[0]+v[1]*t[1]+v[2]*t[2]);
}

double sign(double n)
{
    if(fabs(n) < 10e-12)
		return 0.0;
	else
		return ((n>0)? 1.0 : -1.0);
}

double **zeros(size_t rows, size_t cols)
{
	double **matrix;
  	matrix = (double **)malloc (rows*sizeof(double *));
    for(unsigned int i=0; i<rows; i++)
	{
    		matrix[i] = (double *) calloc (cols, sizeof(double));
  	}
	return(matrix);
}

double det(double m[][3])
{
	return m[0][0]*m[1][1]*m[2][2]+m[0][1]*m[1][2]*m[2][0]+m[0][2]*m[1][0]*m[2][1]-m[2][0]*m[1][1]*m[0][2]-m[2][1]*m[1][2]*m[0][0]-m[2][2]*m[1][0]*m[0][1];
}

double det2x2(double m[][2])
{
	return m[0][0]*m[1][1]-m[1][0]*m[0][1];
}

int roots(double *coef, double **sols_reales)
{
	double *real, *im;
    real = (double*) malloc(16 * sizeof(double));
    im = (double*) malloc(16 * sizeof(double));

	int info[15];
    int solut = rpoly(coef, 16, real, im, info);

    printf("roots() dice que %d\n", solut);
    if(solut < 0)
        return solut;
    int j = 0;
    unsigned int nr_sols = (unsigned int) solut;
    double *v = (double*) calloc(nr_sols, sizeof(double));
	for(int i=0; i<solut; i++)
	{
        if((fabs(im[i]) < 10e-12) && (fabs(real[i]) > 10e-12))
		{
            v[j] = real[i];
			j++;
            printf("%lf+%lf*i\n",real[i],im[i]);
		}
	}

    *sols_reales = v;
    return j;
}

double *cross(double *v1, double *v2)
{
	double *res = (double*) malloc(3 * sizeof(double));
	memcpy(res, (double[3]) {v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]}, sizeof(double[3]));
	return res;
}
