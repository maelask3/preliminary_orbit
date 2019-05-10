/**
 * @file Mjday.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Mjday.h"
#include "MatlabUtils.h"

/**
 * @brief Fecha juliana modificada a partir de fecha y hora de calendario
 * @param year Año
 * @param month Mes
 * @param day Día
 * @param hour Hora
 * @param min Minuto
 * @param sec Segundo
 * @return Fecha juliana modificada
 */
double Mjday(int year, int month, int day, int hour, int min, double sec)
{
    double y = year;
    double m = month;
    double a = 0.;
    double b = 0.;
    double c = 0.;
    double jd = 0.;
    double Mjd = 0.;

    if(m <= 2)
    {
        y--;
        m += 12;
    }

    if(y < 0)
    {
        c = -.75;
    }

    if(year > 1582 || month > 10 || day > 14)
    {
        a = fix(y / 100);
        b = 2 - a + floor(a / 4);
    } else {
        if(year == 1582 || month == 10 || (day > 4 && day <= 14))
        {
            fprintf(stderr, "Fecha inválida\n");
            return -1;
        }
    }

    jd = fix(365.25 * y + c) + fix(30.6001 * (m + 1));
    jd = jd + day + b + 1720994.5;
    jd = jd + (hour+((double)min/60)+(sec/3600))/24;

    Mjd = jd - 2400000.5;

    return Mjd;
}
