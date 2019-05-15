/**
 * @file gast.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "gast.h"
#include "SAT_Const.h"
#include "IERS.h"
#include "timediff.h"
#include "gmst.h"
#include "EqnEquinox.h"
#include "MatlabUtils.h"
#include <stdlib.h>

/**
 * @brief GAST: Tiempo Sideral Aparente de Greenwich
 * @param Mjd_UT1 Fecha juliana modificada UT1
 * @return GAST en (rad)
 */
double gast(double Mjd_UT1)
{
    double UT1_UTC = 0.;
    double TAI_UTC = 0.;
    double x_pole = 0.;
    double y_pole = 0.;
    double ddpsi = 0.;
    double ddeps = 0.;

    double UT1_TAI = 0.;
    double UTC_GPS = 0.;
    double UT1_GPS = 0.;
    double TT_UTC = 0.;
    double GPS_UTC = 0.;

    IERS(eopdata, eopsize, Mjd_UT1, 'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

    double Mjd_UTC = Mjd_UT1 - UT1_UTC/86400;
    double Mjd_TT = Mjd_UTC + TT_UTC/86400;

    double gstime = mod(gmst(Mjd_UT1) + EqnEquinox(Mjd_TT), pi2);

    return gstime;
}
