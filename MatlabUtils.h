/**
 * @file MatlabUrils.h
 * @Autor Davide Pérez y Millán Santamaría
 * @brief Es el fichero de cabeceras de MAtlabUtils.c
 */
#ifndef MATLABUTILS_H
#define MATLABUTILS_H
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

double norm(double *v);
double dot(double *v1, double *v2);
double sign(double n);
double **zeros(unsigned int rows, unsigned int cols);
double det(double M[][3]);
double det2x2(double m[][2]);
int roots(double *coef, double **sols_reales);
double fix(double in);
double mod(double a, double m);
double *cross(double *v1, double *v2);
double **productMatrix(double **m1, double **m2);
double **sumMatrix(double **m1, double **m2);
double **transposeMatrix(double **m);
double *vectorProductDouble(double *v, double d);
double *sumVector(double *v1, double *v2);

#endif // MATLABUTILS_H
