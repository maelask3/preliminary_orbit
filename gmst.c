/**
 * @file gmst.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "gmst.h"
#include "Frac.h"
#include "SAT_Const.h"
#include <math.h>

/**
 * @brief GMST: Tiempo sideral medio de Greenwich
 * @param Mjd_UT1 Fecha juliana modificada UT1
 * @return GMST en (rad)
 */
double gmst(double Mjd_UT1)
{
    double Secs = 86400;

    double Mjd_0 = floor(Mjd_UT1);
    double UT1 = Secs*(Mjd_UT1-Mjd_0);
    double T_0 = (Mjd_0 - MJD_J2000)/36525;
    double T = (Mjd_UT1 - MJD_J2000)/36525;

    double gmst = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1 + (0.093104-6.2e-6*T)*T*T;
    double gmstime = 2*pi*Frac(gmst/Secs);

    return gmstime;
}
