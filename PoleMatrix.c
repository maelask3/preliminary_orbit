/**
 * @file PoleMatrix.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "R_x.h"
#include "R_y.h"
#include "PoleMatrix.h"
#include "MatlabUtils.h"

/**
 * @brief Transformación de coordenadas pseudo-fijadas a la Tierra a fijadas a la Tierra para una fecha dada
 * @param xp Coordenada polar X
 * @param yp Coordenada polar Y
 * @return Matriz polar
 *
 * @note Esta función devuelve un puntero a memoria asignada.
 */
double **PoleMatrix(double xp, double yp)
{
    double **ryxp = R_y(-xp);
    double **rxyp = R_x(-yp);
    double **PoleMat = productMatrix(ryxp, rxyp);
    free(ryxp);
    free(rxyp);

    return PoleMat;
}
