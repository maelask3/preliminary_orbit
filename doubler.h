/**
 * @file doubler.h
 * @author Davide Pérez y Millán Santamaría
 */
#ifndef __DOUBLER_H
#define __DOUBLER_H
#include <stdlib.h>
#include <stdbool.h>
#include "MatlabUtils.h"
#include "SAT_Const.h"

double *doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magr1in, double magr2in, double *los1, double *los2, double *los3, double *rsite1, double *rsite2, double *rsite3, double t1, double t3, char direct, double *r2, double *r3);

#endif // __DOUBLER_H_
