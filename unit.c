/**
 * @file unit.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "unit.h"
#include "MatlabUtils.h"

/**
 * @brief Calcula un vector unitario dado el vector original. Si el vector introducido es nulo, el vector será inicializado a 0.
 * @param vec Vector
 * @return Vector unitario
 *
 * @note Esta función devuelve un puntero a memoria asignada.
 */
double *unit(double *vec)
{
    double small = 0.000001;
    double magv = norm(vec);

    double *outvec = calloc(3, sizeof(double));
    if(magv > small)
    {
        for(int i=0; i<3; i++)
            outvec[i] = vec[i]/magv;
    }

    return outvec;
}
