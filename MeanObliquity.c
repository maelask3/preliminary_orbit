/**
 * @file MeanObliquity.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "SAT_Const.h"
#include "MeanObliquity.h"

/**
 * @brief Calcula la oblicuidad media de la eclíptica
 * @param Mjd_TT Fecha juliana modificada (Tiempo terrestre)
 * @return Oblicuidad media de la eclíptica
 */
double MeanObliquity(double Mjd_TT)
{
    double T = (Mjd_TT - MJD_J2000)/36525;
    double MOblq = Rad*(23.43929111 - (46.8150+(0.00059-0.001813*T)*T)*T/3600);

    return MOblq;
}
