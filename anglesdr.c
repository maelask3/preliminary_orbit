/**
 * @file anglesdr.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "doubler.h"
#include "anglesdr.h"
#include <math.h>
#include "lambert_gooding.h"
#include <stdio.h>

/**
 * @brief Resuelve el problema de la determinación de órbitas utilizando tres observaciones ópticas
 * @param rtasc1 Ascensión derecha en t1 (rad)
 * @param rtasc2 Ascensión derecha en t2 (rad)
 * @param rtasc3 Ascensión derecha en t3 (rad)
 * @param decl1 declinación en t1 (rad)
 * @param decl2 declinación en t2 (rad)
 * @param decl3 declinación en t3 (rad)
 * @param Mjd1 Fecha juliana modificada de t1
 * @param Mjd2 Fecha juliana modificada de t2
 * @param Mjd3 Fecha juliana modificada de t3
 * @param rsite1 Vector de posición ijk de site1 (m)
 * @param rsite2 Vector de posición ijk de site2 (m)
 * @param rsite3 Vector de posición ijk de site3 (m)
 * @param r2 (Salida) Vector de posición ijk en t2 (m)
 * @param v2 (Salida) Vector de velocidad ijk en t2 (m/s)
 */
void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1, double decl2, double decl3, double Mjd1, double Mjd2, double Mjd3, double *rsite1, double *rsite2, double *rsite3, double **r2, double **v2){
	
	double magr1in = 2.01*R_Earth;
	double magr2in = 2.11*R_Earth;
	char direct  = 'y';

	double tol    = 1e-8*R_Earth;
	double pctchg = 5e-6;

	double t1 = (Mjd1 - Mjd2)*86400;
	double t3 = (Mjd3 - Mjd2)*86400;

	double los1[] = {cos(decl1)*cos(rtasc1), cos(decl1)*sin(rtasc1), sin(decl1)};
	double los2[] = {cos(decl2)*cos(rtasc2), cos(decl2)*sin(rtasc2), sin(decl2)};
	double los3[] = {cos(decl3)*cos(rtasc3), cos(decl3)*sin(rtasc3), sin(decl3)};

	double magr1old  = 99999.0e3;
	double magr2old  = 99999.0e3;

	double magrsite1 = norm(rsite1);
	double magrsite2 = norm(rsite2);
    //double magrsite3 = norm(rsite3); //error por variable no utilizada

	double cc1 = 2*dot(los1,rsite1);
	double cc2 = 2*dot(los2,rsite2);

	int ll=0;
	
    double f1, f2, q1, magr1, magr2, a, deltae32, f, g, magr1o=0., deltar1, f1delr1, f2delr1, q2, magr2o, deltar2, f1delr2, f2delr2, q3, pf1pr2, pf2pr2, delta, pf1pr1=0., pf2pr1=0., delta1, delta2;
    double *r3 = NULL;
    double *aux = NULL;
	while(fabs(magr1in-magr1old) > tol && fabs(magr2in-magr2old) > tol && ll<=3){
        ll++;
        r3 = calloc(3, sizeof(double));
        aux = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct, *r2, r3);

        f1 = aux[0];
        f2 = aux[1];
        q1 = aux[2];
        magr1 = aux[3];
        magr2 = aux[4];
        a = aux[5];
        deltae32 = aux[6];

        f = 1 - (a/magr2)*(1-cos(deltae32));
        g = t3 - sqrt((a*a*a)/GM_Earth)*(deltae32-sin(deltae32));
        *v2 = vectorProductDouble(sumVector(r3, vectorProductDouble(*r2, -f)), 1./g);

        magr1o = magr1in;
        magr1in = (1+pctchg)*magr1in;
        deltar1 = pctchg*magr1in;

        free(aux);
        aux = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct, *r2, r3);

        f1delr1 = aux[0];
        f2delr1 = aux[1];
        q2 = aux[2];
        magr1 = aux[3];
        magr2 = aux[4];
        a = aux[5];
        deltae32 = aux[6];

        pf1pr1 = (f1delr1-f1)/deltar1;
        pf2pr1 = (f2delr1-f2)/deltar1;

        magr1in = magr1o;
        deltar1 = pctchg*magr1in;
        magr2o = magr2in;
        magr2in = (1+pctchg)*magr2in;
        deltar2 = pctchg*magr2in;

        free(aux);
        aux = doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2, los3, rsite1, rsite2, rsite3, t1, t3, direct, *r2, r3);

        f1delr2 = aux[0];
        f2delr2 = aux[1];
        q2 = aux[2];
        magr1 = aux[3];
        magr2 = aux[4];
        a = aux[5];
        deltae32 = aux[6];

        pf1pr2 = (f1delr2 - f1)/deltar2;
        pf2pr2 = (f2delr2 - f2)/deltar2;

        magr2in = magr2o;
        deltar2 = pctchg*magr2in;

        delta = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
        delta1 = pf2pr2*f1 - pf1pr2*f2;
        delta2 = pf1pr1*f2 - pf2pr1*f1;

        deltar1 = -delta1/delta;
        deltar2 = -delta2/delta;

        magr1old = magr1in;
        magr2old = magr2in;

        magr1in += deltar1;
        magr2in += deltar2;

        printf("=== %d ===\nmagr1in = %lf\nmarg2in = %lf\n",ll,magr1in, magr2in);
	}

    aux = doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, *r2, r3);
	f1 = aux[0];
	f2 = aux[1];
	q1 = aux[2];
	magr1 = aux[3];
	magr2 = aux[4];
	a = aux[5]; 
	deltae32 = aux[6];

    double *v3 = calloc(3, sizeof(double));

    lambert_gooding(*r2, r3, (Mjd3-Mjd2)*86400, GM_Earth, 0, 1, *v2, v3);

	f  = 1 - a/magr2*(1-cos(deltae32));
	g  = t3 - sqrt((a*a*a)/GM_Earth)*(deltae32-sin(deltae32));
    *v2 = vectorProductDouble(sumVector(r3, vectorProductDouble(*r2, -1.0*f)),(1.0/g));
}
