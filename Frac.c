/**
 * @file Frac.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "Frac.h"
#include <math.h>
/**
 * @brief Parte fraccional de un número
 * @param x Número
 * @return Devuelve la parte fraccional del número (y = x - [x])
 */
double Frac(double x)
{
    return x-floor(x);
}
