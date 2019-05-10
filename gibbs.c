/**
 * @file gibbs.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "gibbs.h"
#include "MatlabUtils.h"
#include "SAT_Const.h"
#include "unit.h"
#include "angl.h"
#include <math.h>
/**
 * @brief Realiza el método de gibbs de determinación de órbitas. Este método determina la velocidad en el punto medio de 3 vectores de posición dados.
 * @param r1 Vector de posición ijk No.1 (m)
 * @param r2 Vector de posición ijk No.1 (m)
 * @param r3 Vector de posición ijk No.1 (m)
 * @param v2 (Salida) Vector de velocidad ijk para r2 (m/s)
 * @param theta (Salida) Ángulo entre vectores (rad)
 * @param theta1 (Salida)
 * @param copa (Salida)
 * @param error (Salida) Cadena indicando el éxito de la operación ("ok",...)
 */
void gibbs(double *r1, double *r2, double *r3, double **v2, double *theta, double *theta1, double *copa, char **error)
{
    double small = 0.00000001;
    *theta = 0.;
    *theta1 = 0.;
    *error = "          ok";

    double magr1 = norm(r1);
    double magr2 = norm(r2);
    double magr3 = norm(r3);

    double *p = cross(r2, r3);
    double *q = cross(r3, r1);
    double *w = cross(r1, r2);

    double *pn = unit(p);
    double *r1n = unit(r1);
    *copa = asin(dot(pn, r1n));

    if(fabs(dot(r1n, pn)) > 0.017452406)
    {
        *error = "not coplanar";
    }

    double *d = sumVector(p, sumVector(q, w));
    double magd = norm(d);
    double *n = sumVector(vectorProductDouble(p, magr1),sumVector(vectorProductDouble(q, magr2), vectorProductDouble(w, magr3)));
    double magn = norm(n);
    double *nn = unit(n);
    double *dn = unit(d);

    if((fabs(magd) < small) || (fabs(magn) < small) || (dot(nn, dn) < small))
    {
        *error = "impossible";
    } else {
        *theta = angl(r1, r2);
        *theta1 = angl(r2, r3);

        double r1mr2 = magr1-magr2;
        double r3mr1 = magr3-magr1;
        double r2mr3 = magr2-magr3;
        double *s = sumVector(vectorProductDouble(r3, r1mr2), sumVector(vectorProductDouble(r2, r3mr1), vectorProductDouble(r1, r2mr3)));
        double *b = cross(d, r2);
        double l = sqrt(GM_Earth/(magd*magn));
        double tover2 = l/magr2;
        *v2 = sumVector(vectorProductDouble(b, tover2), vectorProductDouble(s, l));
        free(s);
        free(b);
    }

    free(p);
    free(q);
    free(w);
    free(pn);
    free(r1n);
    free(d);
    free(n);
}
