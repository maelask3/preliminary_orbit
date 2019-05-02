#ifndef __MATLABUTILS_H
#define __MATLABUTILS_H
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>

typedef struct {
    double *data;
    int length;
} double_arr;

double norm(double *v);
double dot(double *v1, double *v2);
double sign(double n);
double **zeros(size_t rows, size_t cols);
double det(double M[][3]);
double det2x2(double m[][2]);
double_arr roots(double *coef);
double *cross(double *v1, double *v2);

// PRODUCTO DOS MATRICES
// SUMA MATRICES
// TRANSPUESTA MATRIZ

#endif // __MATLABUTILS_H_
