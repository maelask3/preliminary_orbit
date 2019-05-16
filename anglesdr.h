/**
 * @file anglesdr.h
 * @author Davide Pérez y Millán Santamaría
 * @brief Definición de anglesdr()
 */

#ifndef __ANGLESDR_H
#define __ANGLESDR_H
#include <stdlib.h>
#include "MatlabUtils.h"
#include  "SAT_Const.h"

void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1, double decl2, double decl3, double Mjd1, double Mjd2, double Mjd3, double *rsite1, double *rsite2, double *rsite3, double *r2, double *v2);

#endif // __ANGLESDR_H_
