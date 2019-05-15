/**
 * @file NutMatrix.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "NutMatrix.h"
#include "R_x.h"
#include "R_z.h"
#include "NutAngles.h"
#include "MeanObliquity.h"
#include "MatlabUtils.h"

/**
 * @brief Transformación de ecuador y equinoccio medios a reales
 * @param Mjd_TT Fecha juliana modificada (Tiempo terrestre)
 * @return Matriz de nutación
 *
 * @note Esta función devuelve un puntero a memoria asignada.
 */
double **NutMatrix(double Mjd_TT)
{
    double ep = MeanObliquity(Mjd_TT);
    double dpsi = 0.;
    double deps = 0.;
    NutAngles(Mjd_TT, &dpsi, &deps);

    double **rxep = R_x(ep);
    double **rzdpsi = R_z(-dpsi);
    double **rhs = productMatrix(rzdpsi, rxep);
    free(rxep);
    free(rzdpsi);
    double **rxepdeps = R_x(-ep-deps);

    double **NutMat = productMatrix(rxepdeps, rhs);
    free(rhs);
    free(rxepdeps);

    return NutMat;
}
