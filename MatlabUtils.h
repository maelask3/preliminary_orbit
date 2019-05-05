#ifndef __MATLABUTILS_H
#define __MATLABUTILS_H

#include <stdbool.h>


double norm(double *v);
double dot(double *v1, double *v2);
double sign(double n);
double **zeros(int rows, int cols);
double det(double M[][3]);
double det2x2(double m[][2]);
int roots(double *coef, double *sols_reales);
double *cross(double *v1, double *v2);
double **productMatrix(double **m1, double **m2);
double **sumMartix(double **m1, double **m2);
double **transposeMatrix(double **m);

// PRODUCTO DOS MATRICES
// SUMA MATRICES
// TRANSPUESTA MATRIZ

#endif // __MATLABUTILS_H_
