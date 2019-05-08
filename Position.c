#include "Position.h"
#include <math.h>
#include "SAT_Const.h"
#include <stdlib.h>

double *Position(double lon, double lat, double h)
{
    double R_equ = R_Earth;
    double f = f_Earth;

    double e2 = f*(2-f);
    double CosLat = cos(lat);
    double SinLat = sin(lat);

    double N = R_equ/sqrt(1-e2*SinLat*SinLat);

    double *r = calloc(3, sizeof(double));
    r[0] = (N+h)*CosLat*cos(lon);
    r[1] = (N+h)*CosLat*sin(lon);
    r[2] = ((1-e2)*N+h)*SinLat;

    return r;
}
