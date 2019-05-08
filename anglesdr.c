/**
 * @file MatlabUtils.c
 * @Autor Davide Pérez y Millán Santamaría
 * @brief Es una libreria de funciones presentes en Matlab y no en C.
 */
#include "doubler.h"
#include <stdio.h>
#include <math.h>

void anglesdr(double rtasc1, double rtasc2, double rtasc3, double decl1, double decl2, double decl3, double Mjd1, double Mjd2, double Mjd3, double *rsite1, double *rsite2, double *rsite3, double *r2, double *v2){
	
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
	double magrsite3 = norm(rsite3);

	double cc1 = 2*dot(los1,rsite1);
	double cc2 = 2*dot(los2,rsite2);

	int ll=0;
	
	double f1, f2, q1, magr1, magr2, a, deltae32, f, g, magr1o, deltar1, magr2o, deltar2, f1delr2, f2delr2, q3, pf1pr2, pf2pr2, delta, pf1pr1, pf2pr1, delta1, delta2;
	double *r3;
	while(fabs(magr1in-magr1old) > tol && fabs(magr2in-magr2old) > tol && ll<=3){
		ll = ll+1;
		double *r3 = (double*) malloc(3*sizeof(double));
		double *aux = doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, r2, r3);
		f1 = aux[0];
		f2 = aux[1];
		q1 = aux[2];
		magr1 = aux[3];
		magr2 = aux[4];
		a = aux[5];
		deltae32 = aux[6];

		f  = 1 - a/magr2*(1-cos(deltae32));
		g  = t3 - sqrt(a*a*a/GM_Earth)*(deltae32-sin(deltae32));
		v2 = vectorProductDouble(sumVector(r3, vectorProductDouble(r2, -1.0*f)),(1.0/g));

		magr1in = magr1o;
		deltar1 = pctchg*magr1in;
		magr2o = magr2in;
		magr2in = (1+pctchg)*magr2in;
		deltar2 = pctchg*magr2in;
		aux = doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, r2, r3);
		f1delr2 = aux[0];
		f2delr2 = aux[1];
		q3 = aux[2];
		magr1 = aux[3];
		magr2 = aux[4];
		a = aux[5]; 
		deltae32 = aux [6];
		pf1pr2 = (f1delr2-f1)/deltar2;
		pf2pr2 = (f2delr2-f2)/deltar2;

		magr2in = magr2o;
		deltar2 = pctchg*magr2in;

		delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
		delta1 = pf2pr2*f1 - pf1pr2*f2;
		delta2 = pf1pr1*f2 - pf2pr1*f1;

		deltar1 = -delta1/delta;
		deltar2 = -delta2/delta;

		magr1old = magr1in;
		magr2old = magr2in;

		magr1in = magr1in + deltar1;
		magr2in = magr2in + deltar2;
	}

	double *aux = doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,direct, r2, r3);
	f1 = aux[0];
	f2 = aux[1];
	q1 = aux[2];
	magr1 = aux[3];
	magr2 = aux[4];
	a = aux[5]; 
	deltae32 = aux[6];

	f  = 1 - a/magr2*(1-cos(deltae32));
	g  = t3 - sqrt((a*a*a)/GM_Earth)*(deltae32-sin(deltae32));
	v2 = vectorProductDouble(sumVector(r3, vectorProductDouble(r2, -1.0*f)),(1.0/g));
}
