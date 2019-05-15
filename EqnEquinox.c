/**
 * @file EqnEquinox.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "EqnEquinox.h"
#include "NutAngles.h"
#include "MeanObliquity.h"
#include <math.h>

/**
 * @brief Cómputo de la ecuación de los equinoccios
 * @param Mjd_TT Fecha juliana modificada (Tiempo Terrestre)
 * @return Ecuación de los equinoccios
 *
 * @note La ecuación de los equinoccios dpsi*cos(eps) es la ascensión
 * derecha del equinoccio medio referida al ecuador y equinoccio verdaderos
 * y es igual a la diferencia entre tiempo sideral aparente y medio.
 */
double EqnEquinox(double Mjd_TT)
{
    double dpsi = 0.;
    double deps = 0.;
    NutAngles(Mjd_TT, &dpsi, &deps);

    return dpsi * cos(MeanObliquity(Mjd_TT));
}
