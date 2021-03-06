/**
 * @file NutAngles.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "NutAngles.h"
#include "SAT_Const.h"
#include "MatlabUtils.h"
#include <stdlib.h>
#include <math.h>

/**
 * @brief Nutación en longitud y oblicuidad
 * @param Mjd_TT Fecha juliana modificada (Tiempo terrestre)
 * @param dpsi (Salida) Ángulo de nutación
 * @param deps (Salida) Ángulo de nutación
 */
void NutAngles(double Mjd_TT, double *dpsi, double *deps)
{
    double T = (Mjd_TT-(MJD_J2000))/36525;
    double T2 = T*T;
    double T3 = T2*T;
    double rev = 360 * 3600;

    int N_coeff = 106;

    double **C = malloc(106 * sizeof(double*));
    C[0] = (double*)(double[9]) {0, 0, 0, 0, 1,-1719960,-1742,  920250,   89 };
    C[1] = (double*)(double[9]) {    0, 0, 0, 0, 2,   20620,    2,   -8950,    5 };
    C[2] = (double*)(double[9]) {   -2, 0, 2, 0, 1,     460,    0,    -240,    0 };
    C[3] = (double*)(double[9]) {    2, 0,-2, 0, 0,     110,    0,       0,    0 };
    C[4] = (double*)(double[9]) {   -2, 0, 2, 0, 2,     -30,    0,      10,    0 };
    C[5] = (double*)(double[9]) {    1,-1, 0,-1, 0,     -30,    0,       0,    0 };
    C[6] = (double*)(double[9]) {    0,-2, 2,-2, 1,     -20,    0,      10,    0 };
    C[7] = (double*)(double[9]) {    2, 0,-2, 0, 1,      10,    0,       0,    0 };
    C[8] = (double*)(double[9]) {    0, 0, 2,-2, 2, -131870,  -16,   57360,  -31 };
    C[9] = (double*)(double[9]) {    0, 1, 0, 0, 0,   14260,  -34,     540,   -1 };
    C[10] = (double*)(double[9]) {    0, 1, 2,-2, 2,   -5170,   12,    2240,   -6 };
    C[11] = (double*)(double[9]) {    0,-1, 2,-2, 2,    2170,   -5,    -950,    3 };
    C[12] = (double*)(double[9]) {    0, 0, 2,-2, 1,    1290,    1,    -700,    0 };
    C[13] = (double*)(double[9]) {    2, 0, 0,-2, 0,     480,    0,      10,    0 };
    C[14] = (double*)(double[9]) {    0, 0, 2,-2, 0,    -220,    0,       0,    0 };
    C[15] = (double*)(double[9]) {    0, 2, 0, 0, 0,     170,   -1,       0,    0 };
    C[16] = (double*)(double[9]) {    0, 1, 0, 0, 1,    -150,    0,      90,    0 };
    C[17] = (double*)(double[9]) {    0, 2, 2,-2, 2,    -160,    1,      70,    0 };
    C[18] = (double*)(double[9]) {    0,-1, 0, 0, 1,    -120,    0,      60,    0 };
    C[19] = (double*)(double[9]) {   -2, 0, 0, 2, 1,     -60,    0,      30,    0 };
    C[20] = (double*)(double[9]) {    0,-1, 2,-2, 1,     -50,    0,      30,    0 };
    C[21] = (double*)(double[9]) {    2, 0, 0,-2, 1,      40,    0,     -20,    0 };
    C[22] = (double*)(double[9]) {    0, 1, 2,-2, 1,      40,    0,     -20,    0 };
    C[23] = (double*)(double[9]) {    1, 0, 0,-1, 0,     -40,    0,       0,    0 };
    C[24] = (double*)(double[9]) {    2, 1, 0,-2, 0,      10,    0,       0,    0 };
    C[25] = (double*)(double[9]) {    0, 0,-2, 2, 1,      10,    0,       0,    0 };
    C[26] = (double*)(double[9]) {    0, 1,-2, 2, 0,     -10,    0,       0,    0 };
    C[27] = (double*)(double[9]) {    0, 1, 0, 0, 2,      10,    0,       0,    0 };
    C[28] = (double*)(double[9]) {   -1, 0, 0, 1, 1,      10,    0,       0,    0 };
    C[29] = (double*)(double[9]) {    0, 1, 2,-2, 0,     -10,    0,       0,    0 };
    C[30] = (double*)(double[9]) {    0, 0, 2, 0, 2,  -22740,   -2,    9770,   -5 };
    C[31] = (double*)(double[9]) {    1, 0, 0, 0, 0,    7120,    1,     -70,    0 };
    C[32] = (double*)(double[9]) {    0, 0, 2, 0, 1,   -3860,   -4,    2000,    0 };
    C[33] = (double*)(double[9]) {    1, 0, 2, 0, 2,   -3010,    0,    1290,   -1 };
    C[34] = (double*)(double[9]) {    1, 0, 0,-2, 0,   -1580,    0,     -10,    0 };
    C[35] = (double*)(double[9]) {   -1, 0, 2, 0, 2,    1230,    0,    -530,    0 };
    C[36] = (double*)(double[9]) {    0, 0, 0, 2, 0,     630,    0,     -20,    0 };
    C[37] = (double*)(double[9]) {    1, 0, 0, 0, 1,     630,    1,    -330,    0 };
    C[38] = (double*)(double[9]) {   -1, 0, 0, 0, 1,    -580,   -1,     320,    0 };
    C[39] = (double*)(double[9]) {   -1, 0, 2, 2, 2,    -590,    0,     260,    0 };
    C[40] = (double*)(double[9]) {    1, 0, 2, 0, 1,    -510,    0,     270,    0 };
    C[41] = (double*)(double[9]) {    0, 0, 2, 2, 2,    -380,    0,     160,    0 };
    C[42] = (double*)(double[9]) {    2, 0, 0, 0, 0,     290,    0,     -10,    0 };
    C[43] = (double*)(double[9]) {    1, 0, 2,-2, 2,     290,    0,    -120,    0 };
    C[44] = (double*)(double[9]) {    2, 0, 2, 0, 2,    -310,    0,     130,    0 };
    C[45] = (double*)(double[9]) {    0, 0, 2, 0, 0,     260,    0,     -10,    0 };
    C[46] = (double*)(double[9]) {   -1, 0, 2, 0, 1,     210,    0,    -100,    0 };
    C[47] = (double*)(double[9]) {   -1, 0, 0, 2, 1,     160,    0,     -80,    0 };
    C[48] = (double*)(double[9]) {    1, 0, 0,-2, 1,    -130,    0,      70,    0 };
    C[49] = (double*)(double[9]) {   -1, 0, 2, 2, 1,    -100,    0,      50,    0 };
    C[50] = (double*)(double[9]) {    1, 1, 0,-2, 0,     -70,    0,       0,    0 };
    C[51] = (double*)(double[9]) {    0, 1, 2, 0, 2,      70,    0,     -30,    0 };
    C[52] = (double*)(double[9]) {    0,-1, 2, 0, 2,     -70,    0,      30,    0 };
    C[53] = (double*)(double[9]) {    1, 0, 2, 2, 2,     -80,    0,      30,    0 };
    C[54] = (double*)(double[9]) {    1, 0, 0, 2, 0,      60,    0,       0,    0 };
    C[55] = (double*)(double[9]) {    2, 0, 2,-2, 2,      60,    0,     -30,    0 };
    C[56] = (double*)(double[9]) {    0, 0, 0, 2, 1,     -60,    0,      30,    0 };
    C[57] = (double*)(double[9]) {    0, 0, 2, 2, 1,     -70,    0,      30,    0 };
    C[58] = (double*)(double[9]) {    1, 0, 2,-2, 1,      60,    0,     -30,    0 };
    C[59] = (double*)(double[9]) {    0, 0, 0,-2, 1,     -50,    0,      30,    0 };
    C[60] = (double*)(double[9]) {    1,-1, 0, 0, 0,      50,    0,       0,    0 };
    C[61] = (double*)(double[9]) {    2, 0, 2, 0, 1,     -50,    0,      30,    0 };
    C[62] = (double*)(double[9]) {    0, 1, 0,-2, 0,     -40,    0,       0,    0 };
    C[63] = (double*)(double[9]) {    1, 0,-2, 0, 0,      40,    0,       0,    0 };
    C[64] = (double*)(double[9]) {    0, 0, 0, 1, 0,     -40,    0,       0,    0 };
    C[65] = (double*)(double[9]) {    1, 1, 0, 0, 0,     -30,    0,       0,    0 };
    C[66] = (double*)(double[9]) {    1, 0, 2, 0, 0,      30,    0,       0,    0 };
    C[67] = (double*)(double[9]) {    1,-1, 2, 0, 2,     -30,    0,      10,    0 };
    C[68] = (double*)(double[9]) {   -1,-1, 2, 2, 2,     -30,    0,      10,    0 };
    C[69] = (double*)(double[9]) {   -2, 0, 0, 0, 1,     -20,    0,      10,    0 };
    C[70] = (double*)(double[9]) {    3, 0, 2, 0, 2,     -30,    0,      10,    0 };
    C[71] = (double*)(double[9]) {    0,-1, 2, 2, 2,     -30,    0,      10,    0 };
    C[72] = (double*)(double[9]) {    1, 1, 2, 0, 2,      20,    0,     -10,    0 };
    C[73] = (double*)(double[9]) {   -1, 0, 2,-2, 1,     -20,    0,      10,    0 };
    C[74] = (double*)(double[9]) {    2, 0, 0, 0, 1,      20,    0,     -10,    0 };
    C[75] = (double*)(double[9]) {    1, 0, 0, 0, 2,     -20,    0,      10,    0 };
    C[76] = (double*)(double[9]) {    3, 0, 0, 0, 0,      20,    0,       0,    0 };
    C[77] = (double*)(double[9]) {    0, 0, 2, 1, 2,      20,    0,     -10,    0 };
    C[78] = (double*)(double[9]) {   -1, 0, 0, 0, 2,      10,    0,     -10,    0 };
    C[79] = (double*)(double[9]) {    1, 0, 0,-4, 0,     -10,    0,       0,    0 };
    C[80] = (double*)(double[9]) {   -2, 0, 2, 2, 2,      10,    0,     -10,    0 };
    C[81] = (double*)(double[9]) {   -1, 0, 2, 4, 2,     -20,    0,      10,    0 };
    C[82] = (double*)(double[9]) {    2, 0, 0,-4, 0,     -10,    0,       0,    0 };
    C[83] = (double*)(double[9]) {    1, 1, 2,-2, 2,      10,    0,     -10,    0 };
    C[84] = (double*)(double[9]) {    1, 0, 2, 2, 1,     -10,    0,      10,    0 };
    C[85] = (double*)(double[9]) {   -2, 0, 2, 4, 2,     -10,    0,      10,    0 };
    C[86] = (double*)(double[9]) {   -1, 0, 4, 0, 2,      10,    0,       0,    0 };
    C[87] = (double*)(double[9]) {    1,-1, 0,-2, 0,      10,    0,       0,    0 };
    C[88] = (double*)(double[9]) {    2, 0, 2,-2, 1,      10,    0,     -10,    0 };
    C[89] = (double*)(double[9]) {    2, 0, 2, 2, 2,     -10,    0,       0,    0 };
    C[90] = (double*)(double[9]) {    1, 0, 0, 2, 1,     -10,    0,       0,    0 };
    C[91] = (double*)(double[9]) {    0, 0, 4,-2, 2,      10,    0,       0,    0 };
    C[92] = (double*)(double[9]) {    3, 0, 2,-2, 2,      10,    0,       0,    0 };
    C[93] = (double*)(double[9]) {    1, 0, 2,-2, 0,     -10,    0,       0,    0 };
    C[94] = (double*)(double[9]) {    0, 1, 2, 0, 1,      10,    0,       0,    0 };
    C[95] = (double*)(double[9]) {   -1,-1, 0, 2, 1,      10,    0,       0,    0 };
    C[96] = (double*)(double[9]) {    0, 0,-2, 0, 1,     -10,    0,       0,    0 };
    C[97] = (double*)(double[9]) {    0, 0, 2,-1, 2,     -10,    0,       0,    0 };
    C[98] = (double*)(double[9]) {    0, 1, 0, 2, 0,     -10,    0,       0,    0 };
    C[99] = (double*)(double[9]) {    1, 0,-2,-2, 0,     -10,    0,       0,    0 };
    C[100] = (double*)(double[9]) {    0,-1, 2, 0, 1,     -10,    0,       0,    0 };
    C[101] = (double*)(double[9]) {    1, 1, 0,-2, 1,     -10,    0,       0,    0 };
    C[102] = (double*)(double[9]) {    1, 0,-2, 2, 0,     -10,    0,       0,    0 };
    C[103] = (double*)(double[9]) {    2, 0, 0, 2, 0,      10,    0,       0,    0 };
    C[104] = (double*)(double[9]) {    0, 0, 2, 4, 2,     -10,    0,       0,    0 };
    C[105] = (double*)(double[9]) {    0, 1, 0, 1, 0,      10,    0,       0,    0 };

    double l = mod(485866.733 + (1325.0 * rev + 715922.633)*T + 31.310*T2 + 0.064*T3, rev);
    double lp = mod(1287099.804 + (99.0 * rev + 1292581.224)*T - 0.577*T2 - 0.012*T3, rev);
    double F = mod(335778.877 + (1342.0 * rev + 295263.137)*T - 13.257*T2 + 0.011*T3, rev);
    double D = mod(1072261.307 + (1236.0 * rev + 1105601.328)*T - 6.891*T2 + 0.019*T3, rev);
    double Om = mod(450160.280 - (5.0 * rev + 482890.539)*T + 7.455*T2 + 0.008*T3, rev);

    double dpsi_aux = 0;
    double deps_aux = 0;
    double arg = 0;

    for(int i=0; i<N_coeff; i++)
    {
        arg = ((C[i][0]*l)+(C[i][1]*lp)+(C[i][2]*F)+(C[i][3]*D)+(C[i][4]*Om)) / Arcs;
        dpsi_aux = dpsi_aux + ((C[i][5] + (C[i][6] * T)) * sin(arg));
        deps_aux = deps_aux + ((C[i][7] + (C[i][8] * T)) * cos(arg));
    }

    dpsi_aux = 1e-5 * (dpsi_aux/(Arcs));
    deps_aux = 1e-5 * (deps_aux/(Arcs));

    *dpsi = dpsi_aux;
    *deps = deps_aux;

    free(C);
}
