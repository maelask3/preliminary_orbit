/**
 * @file GHAMatrix.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "GHAMatrix.h"
#include "R_z.h"
#include "gast.h"

/**
 * @brief Transformación de ecuador y equinoccio verdaderos a ecuador tererstre y sistema del meridiano de Greenwich
 * @param Mjd_UT1 Fecha juliana modificada UT1
 * @return Matriz del Ángulo Horario de Greenwich
 *
 * @note Esta función devuelve un puntero a memoria asignada.
 */
double **GHAMatrix(double Mjd_UT1)
{
    double **GHAmat = R_z(gast(Mjd_UT1));
    return GHAmat;
}
