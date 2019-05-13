/**
 * @file rv2coe.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "rv2coe.h"
#include "newtonnu.h"
#include "MatlabUtils.h"
#include "SAT_Const.h"
#include "angl.h"
#include <stdlib.h>
#include <string.h>

/**
 * @brief Encuentra los elementos orbitales clásicos dada la posición ecuatorial geocéntrica y vectores de velocidad
 * @param r Vector de posición ijk (m)
 * @param v Vector de velocidad ijk (m/s)
 * @param p (Salida) Semilatus rectum (m)
 * @param a (Salida) Semieje mayor (m)
 * @param ecc (Salida) Excentricidad
 * @param incl (Salida) Inclinación (0 a pi rad)
 * @param omega (Salida) Longitud del nodo ascendiente (0 a 2pi rad)
 * @param argp (Salida) Argumento del perigeo (0 a 2pi rad)
 * @param nu (Salida) Anomalía real (0 a 2pi rad)
 * @param m (Salida) Anomalía mediana (0 a 2pi rad)
 * @param arglat (Salida) Argumento de la latitud (ci) (0 a 2pi rad)
 * @param truelon (Salida) Longitud real (ce) (0 a 2pi rad)
 * @param lonper (Salida) Longitud del periastro (ee) (0 a 2pi rad)
 */

void rv2coe(double *r, double *v, double *p, double *a, double *ecc, double *incl, double *omega, double *argp, double *nu, double *m, double *arglat, double *truelon, double *lonper)
{
    double mu = 398600.4418e9;
    double small = 1e-10;
    double undefined = 999999.1;

    double magr = norm(r);
    double magv = norm(v);

    double *hbar = cross(r, v);
    double magh = norm(hbar);

    if(magh > small)
    {
        double *nbar = calloc(3, sizeof(double));
        double *ebar = calloc(3, sizeof(double));

        nbar[0] = -hbar[1];
        nbar[1] = hbar[0];

        double magn = norm(nbar);
        double c1 = magv*magv - mu/magr;

        double rdotv = dot(r, v);

        for(int i=0;i<3;i++)
        {
            ebar[i] = (c1*r[i] - rdotv*v[i])/mu;
        }

        *ecc = norm(ebar);

        double sme = (magv*magv*0.5) - (mu/magr);
        if(fabs(sme) > small)
        {
            *a = (-mu)/(2.0*sme);
        } else {
            *a = HUGE_VAL;
        }

        *p = magh*magh/mu;

        double hk = hbar[2]/magh;
        *incl = acos(hk);

        char typeorbit[3] = "ei";

        if(*ecc < small)
        {
            if((*incl < small) || (fabs(*incl-pi)<small))
            {
                strncpy(typeorbit, "ce", 3);
                typeorbit[2] = '\0';
            } else {
                strncpy(typeorbit, "ci", 3);
                typeorbit[2] = '\0';
            }
        } else {
            if((*incl < small) || (fabs(*incl-pi)<small))
            {
                strncpy(typeorbit, "ee", 3);
                typeorbit[2] = '\0';
            }
        }

        if(magn > small)
        {
            double temp = nbar[0]/magn;

            if(fabs(temp) > 1.0)
                temp = sign(temp);

            *omega = acos(temp);

            if(nbar[1] < 0.0)
                *omega = pi2 - *omega;
        } else {
            *omega = undefined;
        }

        if(strcmp(typeorbit, "ei")==0)
        {
            *argp = angl(nbar, ebar);
            if(ebar[2] < 0.0)
                *argp = pi2 - *argp;
        } else {
            *argp = undefined;
        }

        if(typeorbit[0] == 'e')
        {
            *nu = angl(ebar, r);
            if(rdotv < 0.0)
                *nu = pi2 - *nu;
        } else {
            *nu = undefined;
        }

        if(strcmp(typeorbit, "ci")==0)
        {
            *arglat = angl(nbar, r);
            if(r[2] < 0.0)
                *arglat = pi2 - *arglat;

            *m = *arglat;
        } else {
            *arglat = undefined;
        }

        if((*ecc > small) && (strcmp(typeorbit, "ee")==0))
        {
            double temp = ebar[0]/(*ecc);
            if(fabs(temp) > 0.0)
                temp = sign(temp);

            *lonper = acos(temp);

            if(ebar[1] < 0.0)
                *lonper = pi2 - *lonper;

            if(*incl > halfpi)
                *lonper = pi2 - *lonper;
        } else {
            *lonper = undefined;
        }

        if((magr > small) && (strcmp(typeorbit, "ce")==0))
        {
            double temp = r[0]/magr;

            if(fabs(temp) > 1.0)
                temp = sign(temp);

            *truelon = acos(temp);

            if(r[1] < 0.0)
                *truelon = pi2 - *truelon;

            if(*incl > halfpi)
                *truelon = pi2 - *truelon;

            *m = *truelon;
        } else {
            *truelon = undefined;
        }

        if(typeorbit[0] == 'e')
        {
            double e = 0.;
            newtonnu(*ecc, *nu, &e, m);
        }

        free(nbar);
        free(ebar);
    } else {
        *p = undefined;
        *a = undefined;
        *ecc = undefined;
        *incl = undefined;
        *omega = undefined;
        *argp = undefined;
        *nu = undefined;
        *m = undefined;
        *arglat = undefined;
        *truelon = undefined;
        *lonper = undefined;
    }

    free(hbar);
}
