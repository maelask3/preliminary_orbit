#ifndef __MATLABUTILS_H_
#define __MATLABUTILS_H_

#include <stdbool.h>
#include <complex.h>

double norm(double *v);
double dot(double *v1, double *v2);
double sign(double n);
double **zeros(int rows, int cols);
double det(double M[][3]);
double det2x2(double m[][2]);
double *roots(double *v);
bool isreal(double complex z);
double *cross(double *v1, double *v2);

// PRODUCTO DOS MATRICES
// SUMA MATRICES
// TRANSPUESTA MATRIZ

#endif // __MATLABUTILS_H_
