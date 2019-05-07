#include "MatlabUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

double **zeros(unsigned int rows, unsigned int cols)
{
	double **matrix;
    matrix = malloc (rows*sizeof(double *));
    for(unsigned int i=0; i<rows; i++)
	{
            //matrix[i] = calloc (cols, sizeof(double));
            size_t colsz = cols * sizeof(double);
            matrix[i] = malloc(colsz);
            for(unsigned int j=0; j<cols; j++)
            {
                matrix[i][j] = 0.;
            }
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
	double *real, *imagin, *aux;
    real = malloc(15*sizeof(double));
    imagin = malloc(15*sizeof(double));
    int solutions = real_poly_roots(coef, 15, real, imagin);
	
	if(solutions < 0)
		return solutions;

    aux = calloc(solutions, sizeof(double));
    int j = 0;
    for(int i = 0; i< solutions; i++)
    {

            if((fabs(imagin[i]) < 10e-12) && (fabs(real[i]) > 10e-12))
            {
                    aux[j] = real[i];
                    j++;
            }
    }
    *sols_reales = aux;

    free(real);
    free(imagin);
    return j;
}


double *cross(double *v1, double *v2)
{
    double *res = malloc(3 * sizeof(double));
	memcpy(res, (double[3]) {v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]}, sizeof(double[3]));
	return res;
}


double **productMatrix(double **m1, double **m2){
	double **product = zeros(3, 3);
	int i, j;
	double sum = 0;
    for (i = 0; i < 3; i++){
            for (j = 0; j < 3; j++){
        		sum += m1[i][j] * m2[j][i];
    		}
            product[i][j] = sum;
    		sum = 0;
	}
	return product;
}
double **sumMatrix(double **m1, double **m2){
	double **sum = zeros(3, 3);
	for(int i = 0; i<3; i++){
		for(int j = 0; j<3; j++){
			sum[i][j] = m1[i][j] + m2[i][j];
		}
	}
	return sum;
}
double **transposeMatrix(double **m){
	double **transpose = zeros(3, 3);
	for(int i = 0; i<3; i++){
		for(int j = 0; j<3; j++){
			transpose[i][j] = m[j][i];
		}
	}
	return transpose;
}
