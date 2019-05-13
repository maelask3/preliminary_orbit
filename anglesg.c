/**
 * @file anglesg.c
 * @authors Davide Pérez y Millán Santamaría
 */

#include "MatlabUtils.h"
#include "SAT_Const.h"
#include <math.h>
#include <string.h>
#include "gibbs.h"
#include "hgibbs.h"
#include "rv2coe.h"
#include "angl.h"
#include "lambert_gooding.h"

void anglesg(double Alpha1, double Alpha2, double Alpha3, double Delta1, double Delta2, double Delta3, double JD1, double JD2, double JD3, double *RS1, double *RS2, double *RS3, double *R2, double *V2){
	double Mu = 398600.4418e9;
	double Rad_ = 180.0/pi;

	double R1[3], R3[3], L1[3], L2[3], L3[3];

	double Tau1 = (JD1-JD2)*86400.0;
	double Tau3 = (JD3-JD2)*86400.0;

	L1[0] = cos(Delta1)*cos(Alpha1);
	L1[1] = cos(Delta1)*sin(Alpha1);
	L1[2] = sin(Delta1);

	L2[0] = cos(Delta2)*cos(Alpha2);
	L2[1] = cos(Delta2)*sin(Alpha2);
	L2[2] = sin(Delta2);

	L3[0] = cos(Delta3)*cos(Alpha3);
	L3[1] = cos(Delta3)*sin(Alpha3);
	L3[2] = sin(Delta3);

	double **LMatIi = zeros(3, 3);
	double **RSMat = zeros(3, 3);
	double **LMatI = zeros(3, 3);
	for(int i=0; i<3; i++)
	{
		LMatIi[i][0] =L1[i];
		LMatIi[i][1] =L2[i];
		LMatIi[i][2] =L3[i];
		RSMat[i][0] =RS1[i];
		RSMat[i][1] =RS2[i];
		RSMat[i][2] =RS3[i];
	}

	double D = det(LMatIi);

	LMatI[0][0] = ( L2[1]*L3[2]-L2[2]*L3[1])/D;
	LMatI[1][0] = (-L1[1]*L3[2]+L1[2]*L3[1])/D;
	LMatI[2][0] = ( L1[1]*L2[2]-L1[2]*L2[1])/D;
	LMatI[0][1] = (-L2[0]*L3[2]+L2[2]*L3[0])/D;
	LMatI[1][1] = ( L1[0]*L3[2]-L1[2]*L3[0])/D;
	LMatI[2][1] = (-L1[0]*L2[2]+L1[2]*L2[0])/D;
	LMatI[0][2] = ( L2[0]*L3[1]-L2[1]*L3[0])/D;
	LMatI[1][2] = (-L1[0]*L3[1]+L1[1]*L3[0])/D;
	LMatI[2][2] = ( L1[0]*L2[1]-L1[1]*L2[0])/D;

	double **LIR = productMatrix(LMatI, RSMat);

	double a1  = Tau3/(Tau3 - Tau1);
	double a1u = (Tau3*((Tau3-Tau1)*(Tau3-Tau1) - Tau3*Tau3 ))/(6.0*(Tau3 - Tau1));
	double a3  = -Tau1 / (Tau3 - Tau1);
	double a3u = -(Tau1*((Tau3-Tau1)*(Tau3-Tau1) - Tau1*Tau1 ))/(6.0*(Tau3 - Tau1));

	double D1 = LIR[1][0]*a1 - LIR[1][1] + LIR[1][2]*a3;
    double D2 = LIR[1][0]*a1u + LIR[1][2]*a3u; // esta wea se va por 13

	double L2DotRS= dot(L2,RS2);
	double magRS2 = norm(RS2);

	double *Poly = (double*)(double[16]) {1.0, 0.0, -(D1*D1 + 2.0*D1*L2DotRS + magRS2*magRS2), 0.0, 0.0, -2.0*Mu*(L2DotRS*D2 + D1*D2), 0.0, 0.0, -Mu*Mu*D2*D2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    double *rootarr = calloc(16, sizeof(double));
    int len_rootarr = roots(Poly, &rootarr); // se va 10^-6

	double BigR2 = 0.0;

	for(int i=0; i<len_rootarr; i++)
	{
		if(rootarr[i] > BigR2)
		{
			BigR2 = rootarr[i];
		}
	}

	double u = Mu/(BigR2*BigR2*BigR2);

	double c1 = a1+a1u*u;
	double c3 = a3+a3u*u;
	double CMat[] = {-c1, 1.0, -c3};
	double RhoMat[3];
	RhoMat[0] = LIR[0][0]*CMat[0]+LIR[0][1]*CMat[1]+LIR[0][2]*CMat[2];
	RhoMat[1] = LIR[1][0]*CMat[0]+LIR[1][1]*CMat[1]+LIR[1][2]*CMat[2];
	RhoMat[2] = LIR[2][0]*CMat[0]+LIR[2][1]*CMat[1]+LIR[2][2]*CMat[2];

	double Rhoold2 = -RhoMat[1];

	double Rho2 = 999999e3;
	int ll = 0;
	while((fabs(Rhoold2-Rho2)>1e-12) && (ll<=2))
	{
		ll++;
		Rho2 = Rhoold2;

		for(int i=0; i<3; i++)
		{
			R1[i] =  RhoMat[0]*L1[i]/c1 + RS1[i];
			R2[i] = -RhoMat[1]*L2[i]    + RS2[i];
			R3[i] =  RhoMat[2]*L3[i]/c3 + RS3[i];
		}
		double **v2 = malloc(sizeof(double*));
		double theta = 0.0;
		double theta1 = 0.0;
		double copa = 0.0;
		char **error = malloc(sizeof(char*));
		gibbs(R1, R2, R3, v2, &theta, &theta1, &copa, error);
        // esto es un bug que tiene Meysam en su propio código.
        if(strcmp(*error,"ok")!=0 && (copa < 1.0/Rad_))
		{
			hgibbs(R1,R2,R3,JD1,JD2,JD3, v2, &theta, &theta1, &copa, error);
		}
		double v1[3];
		lambert_gooding(R1,R2,(JD2-JD1)*86400,Mu,0,1, v1, V2);
		double p, a, ecc, incl, omega, argp, Nu, m, u, l, ArgPer;
		rv2coe(R2, V2, &p, &a, &ecc, &incl, &omega, &argp, &Nu, &m, &u, &l, &ArgPer);
		double magR2 = norm(R2);
		
		double U, RDot, UDot, TauSqr, f1, g1, f3, g3, Theta, Theta1, magR1, magR3;
		if(ll <= 2)
		{
			U = Mu/(magR2*magR2*magR2);
			RDot = dot(R2,V2)/magR2;
			UDot = (-3.0*Mu*RDot)/(magR2*magR2*magR2*magR2);

			TauSqr= Tau1*Tau1;
			f1 =  1.0 - 0.5*U*TauSqr -(1.0/6.0)*UDot*TauSqr*Tau1 + (1.0/24.0) * U*U*TauSqr*TauSqr + (1.0/30.0)*U*UDot*TauSqr*TauSqr*Tau1;
			g1 = Tau1 - (1.0/6.0)*U*Tau1*TauSqr - (1.0/12.0) * UDot*TauSqr*TauSqr + (1.0/120.0)*U*U*TauSqr*TauSqr*Tau1 + (1.0/120.0)*U*UDot*TauSqr*TauSqr*TauSqr;
			TauSqr = Tau3*Tau3;
			f3 =  1.0 - 0.5*U*TauSqr -(1.0/6.0)*UDot*TauSqr*Tau3 + (1.0/24.0) * U*U*TauSqr*TauSqr + (1.0/30.0)*U*UDot*TauSqr*TauSqr*Tau3;
			g3 = Tau3 - (1.0/6.0)*U*Tau3*TauSqr - (1.0/12.0) * UDot*TauSqr*TauSqr + (1.0/120.0)*U*U*TauSqr*TauSqr*Tau3 + (1.0/120.0)*U*UDot*TauSqr*TauSqr*TauSqr;
		}else{
			Theta = angl(R1, R2);
			Theta1 = angl(R2, R3);
			magR1 = norm(R1);
			magR3 = norm(R3);

			f1 = 1.0 - ( (magR1*(1.0 - cos(Theta))/p ) );
			g1 = ( magR1*magR2*sin(-theta) ) / sqrt(p);
			f3 = 1.0 - ( (magR3*(1.0 - cos(Theta1))/p ) );
			g3 = ( magR3*magR2*sin(theta1) )/sqrt(p);
		}
		double c1 = g3/(f1*g3 - f3*g1);
		double c3 = -g1/(f1*g3 - f3*g1);
		
		CMat[0] = -c1;
		CMat[1] = 1.0;
		CMat[2] = -c3;
		RhoMat[0] = LIR[0][0]*CMat[0]+LIR[0][1]*CMat[1]+LIR[0][2]*CMat[2];
		RhoMat[1] = LIR[1][0]*CMat[0]+LIR[1][1]*CMat[1]+LIR[1][2]*CMat[2];
		RhoMat[2] = LIR[2][0]*CMat[0]+LIR[2][1]*CMat[1]+LIR[2][2]*CMat[2];

		double Rhoold2 = -RhoMat[1];
	}
	for(int i=0; i<3; i++)
	{
		R1[i] =  RhoMat[0]*L1[i]/c1 + RS1[i];
		R2[i] = -RhoMat[1]*L2[i] + RS2[i];
		R3[i] =  RhoMat[2]*L3[i]/c3 + RS3[i];
	}
	free(LMatIi);
	free(RSMat);
	free(LMatI);
	free(LIR);
	free(rootarr);
}


