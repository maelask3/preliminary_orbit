/**
 * @file ProjectTest.c
 * @authors Davide Pérez y Millán Santamaría
 * @brief Contiene test unitarios sobre el proyecto
 */
#include "TestUtils.h"
#include "MatlabUtilsTest.h"
#include "Position.h"
#include "Mjday.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "timediff.h"
#include "Frac.h"
#include "R_x.h"
#include "R_y.h"
#include "R_z.h"
#include "unit.h"
#include "EqnEquinox.h"
#include "gmst.h"
#include "NutMatrix.h"
#include "PrecMatrix.h"
#include "PoleMatrix.h"
#include "newtonnu.h"
#include "rv2coe.h"
#include "gibbs.h"
#include "hgibbs.h"
#include "IERS.h"
#include "gast.h"
#include "GHAMatrix.h"
#include <stdlib.h>
#include <stdio.h>
#include "lambert_gooding.h"
#include "doubler.h"
#include "anglesg.h"
double **eopdata = NULL;

int main()
{
	FILE *fp = fopen("eop19620101.txt","r");
    	if(!fp)
    	{
        	fprintf(stderr, "ERROR: No se ha podido abrir eop19620101.txt\n");
        	exit(1);
    	}
    	eopdata = malloc(20026 * sizeof(double*));

    	char line[103];
    	int a = 0, b = 0, c = 0, d = 0, final = 0;
    	float e = 0.f, f =0.f, g= 0.f, h =0.f, m=0.f, j = 0.f, k = 0.f, l = 0.f;
    	for(int i=0; i<20026 && !feof(fp); i++)
    	{
        	eopdata[i] = malloc(13 * sizeof(double));
        	fgets(line, 255, fp);
        	sscanf(line, "%d %d %d %d %f %f %f %f %f %f %f %f %d ", &a,  &b,  &c,  &d,  &e, &f,
                 &g,  &h,  &m,  &j,  &k,  &l,  &final);

        	eopdata[i][0] = a;
        	eopdata[i][1] = b;
       		eopdata[i][2] = c;
        	eopdata[i][3] = d;
        	eopdata[i][4] = (double) e;
        	eopdata[i][5] = (double) f;
        	eopdata[i][6] = (double) g;
        	eopdata[i][7] = (double) h;
        	eopdata[i][8] = (double) m;
        	eopdata[i][9] = (double) j;
        	eopdata[i][10] = (double) k;
        	eopdata[i][11] = (double) l;
        	eopdata[i][12] = final;
    	}

    	fclose(fp);

	//Copiar leer el otro fichero
	
	double lat = Rad*30.5724;
	double lon = Rad*(-86.2143);
	double alt = 0.0;

	double *Rs = Position(lon, lat, alt);

	double Mjd1 = obs[0][0];
	double Mjd2 = obs[1][0];
	double Mjd3 = obs[2][0];

	double Mjd_UTC = Mjd1;
	
    	char interp = 'l';

    	double UT1_UTC = 0.;
    	double TAI_UTC = 0.;
    	double x_pole = 0.;
    	double y_pole = 0.;
    	double ddpsi = 0.;
    	double ddeps = 0.;

	IERS(eopdata, 20026, Mjd_UTC, interp, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);
	double UT1_TAI = 0.;
    	double UTC_GPS = 0.;
    	double UT1_GPS = 0.;
    	double TT_UTC = 0.;
    	double GPS_UTC = 0.;
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	double Mjd_TT = Mjd_UTC + TT_UTC/86400;
	double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

	double **P = PrecMatrix(MJD_J2000,Mjd_TT);
	double **N = NutMatrix(Mjd_TT);
	//Linea 80

	for(int i=0; i<20026; i++)
		free(eopdata[i]);
    	free(eopdata);

	return 0;
}
