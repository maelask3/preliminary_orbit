/**
 * @file IERS.h
 * @authors Davide Pérez y Millán Santamaría
 * @brief Contiene las definiciones de la función IERS y la variable global eopdata
 */
#ifndef IERS_H
#define IERS_H
#include <stdlib.h>
extern double **eopdata;
extern size_t eopsize;
void IERS(double **eop, size_t eop_length, double Mjd_UTC, char interp, double *UT1_UTC, double *TAI_UTC, double *x_pole, double *y_pole, double *ddpsi, double *ddeps);
#endif //IERS_H
