/**
 * @file lambert_gooding.h
 * @Author Davide Pérez y Millán Santamaría
 * @brief Es el fichero de cabeceras de lambert_gooding.c
 */

#ifndef __LAMBERT_GOODING_H
#define __LAMBERT_GOODING_H
#include <stdlib.h>
#include "MatlabUtils.h"

void lambert_gooding(double *r1, double *r2, double tof, double mu, double long_way, double multi_revs, double *v1, double *v2);

double vlamb(double gm, double r1, double r2, double th, double tdelt, double *vri, double *vti, double *vrf, double *vtf);

double *tlamb(double m, double q, double qsqfm1, double x, double n);

double *xlamb(double m, double q, double qsqfm1, double tin);

double d8rt(double x);

#endif // __LAMBERT_GOODING_H_
