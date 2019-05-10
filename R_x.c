/**
 * @file R_x.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "MatlabUtils.h"
#include "R_x.h"
#include <math.h>

/**
 * @brief Matriz de rotación (X)
 * @param angle Ángulo de rotación (rad)
 * @return Matriz de rotación
 *
 * @note Esta función devuelve un puntero a memoria asignada.
 */
double **R_x(double angle)
{
    double C = cos(angle);
    double S = sin(angle);
    double **rotmat = zeros(3,3);

    rotmat[0][0] = 1.; rotmat[0][1] = 0.; rotmat[0][2] = 0.;
    rotmat[1][0] = 0.; rotmat[1][1] = C; rotmat[1][2] = S;
    rotmat[2][0] = 0.; rotmat[2][1] = -1.*S; rotmat[2][2] = C;

    return rotmat;
}
