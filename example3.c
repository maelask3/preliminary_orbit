#include <stdio.h>
#include "IERS.h"
#include <stdlib.h>
#include "Mjday.h"
#include "SAT_Const.h"
#include "Position.h"
#include "timediff.h"
#include "PrecMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "GHAMatrix.h"
#include "MatlabUtils.h"
#include "anglesdr.h"

double **eopdata = NULL;
size_t eopsize = 0;
int main()
{
    FILE *fp;
    fp = fopen("eop19620101.txt","r");
    if(!fp)
    {
        fprintf(stderr, "ERROR: No se ha podido abrir eop19620101.txt\n");
        exit(1);
    }

    eopsize = 20026;
    eopdata = malloc(eopsize * sizeof(double*));

    char line[128];
    int c1 = 0, c2 = 0, c3 = 0, c4 = 0, c13 = 0;
    float c5 = 0.F, c6 =0.F, c7= 0.F, c8 =0.F, c9=0.F, c10 = 0.F, c11 = 0.F, c12 = 0.F;
    for(size_t i=0; i<eopsize; i++)
    {
        eopdata[i] = malloc(13 * sizeof(double));
        fgets(line, 255, fp);
        sscanf(line, "%d %d %d %d %f %f %f %f %f %f %f %f %d ", &c1,  &c2,  &c3,  &c4,  &c5, &c6,
                 &c7,  &c8,  &c9,  &c10,  &c11,  &c12,  &c13);

        eopdata[i][0] = c1;
        eopdata[i][1] = c2;
        eopdata[i][2] = c3;
        eopdata[i][3] = c4;
        eopdata[i][4] = (double) c5;
        eopdata[i][5] = (double) c6;
        eopdata[i][6] = (double) c7;
        eopdata[i][7] = (double) c8;
        eopdata[i][8] = (double) c9;
        eopdata[i][9] = (double) c10;
        eopdata[i][10] = (double) c11;
        eopdata[i][11] = (double) c12;
        eopdata[i][12] = c13;
    }

    fclose(fp);

    fp = fopen("sat3.txt","r");
    if(!fp)
    {
        fprintf(stderr, "ERROR: No se ha podido abrir sat1.txt\n");
        exit(1);
    }
    int Y = 0;
    int M = 0;
    int D = 0;
    int h = 0;
    int m = 0;
    float s = 0.F;
    double rtasc = 0.0;
    double decl = 0.0;

    size_t fsize = 3;

    double **obs = calloc(fsize, sizeof(double*));
    if(!obs)
    {
        fprintf(stderr, "ERROR: No se ha podido asignar memoria suficiente\n");
        fclose(fp);
        exit(3);
    }
    for(size_t i=0; i<fsize && !feof(fp); i++)

    {
        fgets(line, 128, fp);
        sscanf(line, "%d/%d/%d %d:%d:%f %lf %lf", &Y, &M, &D, &h, &m, &s, &rtasc, &decl);

        obs[i] = calloc(3, sizeof(double));
        if(!obs[i])
        {
            fprintf(stderr, "ERROR: No se ha podido asignar memoria suficiente\n");
            fclose(fp);
            exit(3);
        }
        obs[i][0] = Mjday(Y, M, D, h, m, (double) s);
        obs[i][1] = Rad*((double) rtasc);
        obs[i][2] = Rad*((double) decl);
    }

    fclose(fp);

    double lat = Rad*30.5724;
    double lon = Rad*(-86.2143);
    double alt = 0.0;

    double *Rs = Position(lon, lat, alt);

    double Mjd1 = obs[0][0];
    double Mjd2 = obs[1][0];
    double Mjd3 = obs[2][0];

    double Mjd_UTC = Mjd1;

    double UT1_UTC = 0.;
    double TAI_UTC = 0.;
    double x_pole = 0.;
    double y_pole = 0.;
    double ddpsi = 0.;
    double ddeps = 0.;

    IERS(eopdata, eopsize, Mjd_UTC, 'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);

    double UT1_TAI = 0.;
    double UTC_GPS = 0.;
    double UT1_GPS = 0.;
    double TT_UTC = 0.;
    double GPS_UTC = 0.;

    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

    double Mjd_TT = Mjd_UTC + (TT_UTC/86400);
    double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    double **P = PrecMatrix(MJD_J2000, Mjd_TT);
    double **N = NutMatrix(Mjd_TT);

    double **pmxy = PoleMatrix(x_pole, y_pole);
    double **ghaUT1 = GHAMatrix(Mjd_UT1);
    double **NP = productMatrix(N, P);
    double **GNP = productMatrix(ghaUT1, NP);
    free(NP);
    free(ghaUT1);
    double **E = productMatrix(pmxy, GNP);
    free(GNP);
    free(pmxy);

    //Transponer E antes de multiplicar

    double **Et = transposeMatrix(E);
    free(E);

    double *rsite1 = matrixProductVector(Et, Rs);
    free(Et);

    Mjd_UTC = Mjd2;
    IERS(eopdata, eopsize, Mjd_UTC, 'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    Mjd_TT = Mjd_UTC + (TT_UTC/86400);
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(MJD_J2000, Mjd_TT);
    N = NutMatrix(Mjd_TT);

    pmxy = PoleMatrix(x_pole, y_pole);
    ghaUT1 = GHAMatrix(Mjd_UT1);
    NP = productMatrix(N, P);
    GNP = productMatrix(ghaUT1, NP);
    free(NP);
    free(ghaUT1);
    E = productMatrix(pmxy, GNP);
    free(GNP);
    free(pmxy);

    //Transponer E antes de multiplicar

    Et = transposeMatrix(E);
    free(E);

    double *rsite2 = matrixProductVector(Et, Rs);
    free(Et);

    Mjd_UTC = Mjd3;
    IERS(eopdata, eopsize, Mjd_UTC, 'l', &UT1_UTC, &TAI_UTC, &x_pole, &y_pole, &ddpsi, &ddeps);
    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    Mjd_TT = Mjd_UTC + (TT_UTC/86400);
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(MJD_J2000, Mjd_TT);
    N = NutMatrix(Mjd_TT);

    pmxy = PoleMatrix(x_pole, y_pole);
    ghaUT1 = GHAMatrix(Mjd_UT1);
    NP = productMatrix(N, P);
    GNP = productMatrix(ghaUT1, NP);
    free(NP);
    free(ghaUT1);
    E = productMatrix(pmxy, GNP);
    free(GNP);
    free(pmxy);

    //Transponer E antes de multiplicar

    Et = transposeMatrix(E);
    free(E);

    double *rsite3 = matrixProductVector(Et, Rs);
    free(Et);

    double **r2 = malloc(sizeof(double*));
    *r2 = calloc(3, sizeof(double));
    double **v2 = malloc(sizeof(double*));
    *v2 = calloc(3, sizeof(double));

    anglesdr(obs[0][1], obs[1][1], obs[2][1], obs[0][2], obs[1][2], obs[2][2], Mjd1, Mjd2, Mjd3, rsite1, rsite2, rsite3, r2, v2);

    printf("Double-R-Iteration method\nY_apr = \n");
    for(int i=0; i<3; i++)
        printf("\t%lf\n", (*r2)[i]*1e-3);
    for(int i=0; i<3; i++)
        printf("\t%lf\n", (*v2)[i]*1e-3);

    free(*r2);
    free(r2);
    free(*v2);
    free(v2);
    for(size_t i=0; i<fsize; i++)
        free(obs[i]);
    free(obs);

    for(size_t i=0; i<eopsize; i++)
            free(eopdata[i]);
    free(eopdata);
}
