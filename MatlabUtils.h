/**
 * @file MatlabUtils.h
 * @authors Davide Pérez y Millán Santamaría
 * @brief Definiciones de las funciones de MatlabUtils.c
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
double det(double **m);
double det2x2(double **m);
int roots(double *coef, double **sols_reales);
double fix(double in);
double mod(double a, double m);
double *cross(double *v1, double *v2);
double **productMatrix(double **m1, double **m2);
double **sumMatrix(double **m1, double **m2);
double **transposeMatrix(double **m);
double *vectorProductDouble(double *v, double d);
double *sumVector(double *v1, double *v2);
double *matrixProductVector(double **M, double *v);
#endif // MATLABUTILS_H
