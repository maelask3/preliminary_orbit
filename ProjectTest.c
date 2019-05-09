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
#include <stdlib.h>
#include <stdio.h>

void test_Position()
{
    double lon = -2.11796961366573;
    double lat = 0.683053277790977;
    double h = 99.81638;

    double *expected = (double*)(double[3]) {-2577383.6395731382, -4230610.2467987221, 4004108.3320587045};
    double *actual = Position(lon, lat, h);
    array_test_delta("Position() 1", expected, actual, 3, 1e-7);
    free(actual);

    lon = -1.50472339730215;
    lat = 0.533589040236714;
    h = 0.;

    expected = (double*)(double[3]) {362889.51475075335,-5484262.3610134749,3225167.7284776145};
    actual = Position(lon, lat, h);
    array_test_delta("Position() 2", expected, actual, 3, 1e-7);
    free(actual);

    lon = -1.91986217719376;
    lat = 0.698131700797732;
    h = 2000.;

    expected = (double*)(double[3]) {-1673928.5598879098,-4599080.9200722268,4079271.1474197493};
    actual = Position(lon, lat, h);
    array_test_delta("Position() 7", expected, actual, 3, 1e-7);
    free(actual);
}

void test_Mjday()
{
    int year, month, day, hour, min;
    double sec;
    double expected, actual;

    year = 2009;
    month = 5;
    day = 26;
    hour = 16;
    min = 0;
    sec = 20.475;
    expected = 54977.666903645732;
    actual = Mjday(year, month, day, hour, min, sec);
    double_test("Mjday() 1", expected, actual);

    year = 2011;
    month = 1;
    day = 4;
    hour = 13;
    min = 0;
    sec = 46.5;
    expected = 55565.542204861064;
    actual = Mjday(year, month, day, hour, min, sec);
    double_test("Mjday() 2", expected, actual);

    year = 2006;
    month = 9;
    day = 11;
    hour = 4;
    min = 45;
    sec = 44.073;
    expected = 53989.198426770978;
    actual = Mjday(year, month, day, hour, min, sec);
    double_test("Mjday() 3", expected, actual);
}

void test_MeanObliquity()
{
    double Mjd_TT;
    double expected;
    double actual;

    Mjd_TT = 54977.667669664253;
    expected = 0.409071470559;
    actual = MeanObliquity(Mjd_TT);

    double_test("MeanObliquity() 1", expected, actual);

    Mjd_TT = 55565.542970879585;
    expected = 0.409067817510;
    actual = MeanObliquity(Mjd_TT);

    double_test("MeanObliquity() 2", expected, actual);

    Mjd_TT = 53989.199181215423;
    expected = 0.409077612887;
    actual = MeanObliquity(Mjd_TT);

    double_test("MeanObliquity() 3", expected, actual);
}

void test_NutAngles()
{
    double Mjd_TT;
    double exp_dpsi, exp_deps;
    double dpsi = 0;
    double deps = 0;

    Mjd_TT = 54977.667669664253;
    exp_dpsi = 0.000064869339;
    exp_deps = 0.000022305134;

    NutAngles(Mjd_TT, &dpsi, &deps);

    double_test("NutAngles() 1, dpsi", exp_dpsi, dpsi);
    double_test("NutAngles() 1, deps", exp_deps, deps);

    Mjd_TT = 55565.542970879585;
    exp_dpsi = 0.000087228574;
    exp_deps = -0.000000813487;

    NutAngles(Mjd_TT, &dpsi, &deps);

    double_test("NutAngles() 2, dpsi", exp_dpsi, dpsi);
    double_test("NutAngles() 2, deps", exp_deps, deps);

    Mjd_TT = 53989.199181215423;
    exp_dpsi = 0.000007098336;
    exp_deps = 0.000046734356;

    NutAngles(Mjd_TT, &dpsi, &deps);

    double_test("NutAngles() 2, dpsi", exp_dpsi, dpsi);
    double_test("NutAngles() 2, deps", exp_deps, deps);
}

void test_timediff()
{
    double UT1_UTC = 0.258022690875596;
    double TAI_UTC = 34;

    double UT1_TAI_e = 0.;
    double UTC_GPS_e = 0.;
    double UT1_GPS_e = 0.;
    double TT_UTC_e = 0.;
    double GPS_UTC_e = 0.;

    double UT1_TAI = 0.;
    double UTC_GPS = 0.;
    double UT1_GPS = 0.;
    double TT_UTC = 0.;
    double GPS_UTC = 0.;

    UT1_UTC = 0.258022690875596;
    TAI_UTC = 34;

    UT1_TAI_e = -33.7419773091244;
    UTC_GPS_e = -15;
    UT1_GPS_e = -14.7419773091244;
    TT_UTC_e = 66.184;
    GPS_UTC_e = 15;

    timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
    double_test("timediff(), UT1_TAI", UT1_TAI_e, UT1_TAI);
    double_test("timediff(), UTC_GPS", UTC_GPS_e, UTC_GPS);
    double_test("timediff(), UT1_GPS", UT1_GPS_e, UT1_GPS);
    double_test("timediff(), TT_UTC", TT_UTC_e, TT_UTC);
    double_test("timediff(), GPS_UTC", GPS_UTC_e, GPS_UTC);

}

void test_Frac()
{
    double in;
    double expected;
    double actual;

    in = 10.3456636914056;
    expected = 0.3456636914056;
    actual = Frac(in);
    double_test("Frac() 1", expected, actual);

    in = 10.352627149095715;
    expected = 0.352627149095715;
    actual = Frac(in);
    double_test("Frac() 2", expected, actual);
}

void test_R_x()
{
    double in;
    double **expected = malloc(3 * sizeof(double*));
    double **actual;

    in = -0.409093775692299;
    expected[0] = (double*)(double[3]) {1., 0., 0.};
    expected[1] = (double*)(double[3])
    {0., 0.917481675640187, -0.397778047237997};
    expected[2] = (double*)(double[3]) {0., 0.397778047237997, 0.917481675640187};

    actual = R_x(in);

    matrix_test("R_x() 1", expected, actual, 3, 3);

    free(actual);

    in = 0.409071470558628;
    expected[0] = (double*)(double[3]) {1., 0., 0.};
    expected[1] = (double*)(double[3]) {0., 0.917490547904469, 0.397757582587632};
    expected[2] = (double*)(double[3]) {0., -0.397757582587632, 0.917490547904469};

    actual = R_x(in);

    matrix_test("R_x() 2", expected, actual, 3, 3);

    free(actual);

    in = -2.56777193042581e-6;
    expected[0] = (double*)(double[3]) {1., 0., 0.};
    expected[1] = (double*)(double[3]) {0., 0.999999999996703, -2.56777193042298e-6};
    expected[2] = (double*)(double[3]) {0., 2.56777193042298e-6, 0.999999999996703};

    actual = R_x(in);

    matrix_test("R_x() 3", expected, actual, 3, 3);

    free(actual);
    free(expected);
}

void test_R_y()
{
    double in;
    double **expected = malloc(3 * sizeof(double*));
    double **actual;

    in = 0.000913347353936069;
    expected[0] = (double*)(double[3]) {0.999999582898335, 0., -0.000913347226949831};
    expected[1] = (double*)(double[3]) {0., 1., 0.};
    expected[2] = (double*)(double[3]) {0.000913347226949831, 0. , 0.999999582898335};

    actual = R_y(in);

    matrix_test("R_y() 1", expected, actual, 3, 3);

    free(actual);

    in = -7.57892008067929e-8;
    expected[0] = (double*)(double[3]) {0.999999999999997, 0., 7.57892008067929e-8};
    expected[1] = (double*)(double[3]) {0., 1., 0.};
    expected[2] = (double*)(double[3]) {-7.57892008067929e-8, 0., 0.999999999999997};

    actual = R_y(in);

    matrix_test("R_y() 2", expected, actual, 3, 3);

    free(actual);

    in = 0.000913349201373003;
    expected[0] = (double*)(double[3]) {0.999999582896647, 0., -0.000913349074385995};
    expected[1] = (double*)(double[3]) {0., 1., 0.};
    expected[2] = (double*)(double[3]) {0.000913349074385995, 0., 0.999999582896647};

    actual = R_y(in);

    matrix_test("R_y() 3", expected, actual, 3, 3);

    free(actual);
    free(expected);
}

void test_R_z()
{
    double in;
    double **expected = malloc(3 * sizeof(double*));
    double **actual;

    in = -0.001050992069582;
    expected[0] = (double*)(double[3]) {0.999999447707886, -0.00105099187608332, 0.};
    expected[1] = (double*)(double[3]) {0.00105099187608332, 0.999999447707886, 0.};
    expected[2] = (double*)(double[3]) {0., 0., 1.};

    actual = R_z(in);

    matrix_test("R_z() 1", expected, actual, 3, 3);

    free(actual);

    in = 2.17192854407046;
    expected[0] = (double*)(double[3]) {-0.565576570524664, 0.824695788077977, 0.};
    expected[1] = (double*)(double[3]) {-0.824695788077977, -0.565576570524664, 0.};
    expected[2] = (double*)(double[3]) {0., 0., 1.};

    actual = R_z(in);

    matrix_test("R_z() 2", expected, actual, 3, 3);

    free(actual);

    in = 2.50337730618796;
    expected[0] = (double*)(double[3]) {-0.803160266383483, 0.595763028814992, 0.};
    expected[1] = (double*)(double[3]) {-0.595763028814992, -0.803160266383483, 0.};
    expected[2] = (double*)(double[3]) {0., 0., 1.};

    actual = R_z(in);

    matrix_test("R_z() 3", expected, actual, 3, 3);

    free(actual);
    free(expected);
}

void test_unit()
{
    double *vec = (double*)(double[3]){-9341115904217.03, 16158801718408.7, 29720729511155.2};
    double *actual = unit(vec);
    double *expected = (double*)(double[3]){-0.266163758239956, 0.460425439329396, 0.846856108567397};

    array_test("unit() 1", expected, actual, 3);

    vec = (double*)(double[3]){2644218480744.88, -9104674845836.32, 19019299944092.7};
    actual = unit(vec);
    expected = (double*)(double[3]) {0.124425783694902, -0.428427647423945, 0.894968142044561};

    array_test("unit() 2", expected, actual, 3);

    vec = (double*)(double[3]){218961962953.224, 11038749223.2812, 38367764388851.1};
    actual = unit(vec);
    expected = (double*)(double[3]) {0.00570683207063089, 0.000287704227394504, 0.999983674513737};

    array_test_delta("unit() 3", expected, actual, 3, 1e-10);
    // La norma sale chunga, perdemos precision
}

void test_EqnEquinox()
{
    double in;
    double expected;
    double actual;

    in = 54977.6676696643;
    expected = 5.95170051422054e-5;
    actual = EqnEquinox(in);

    double_test("EqnEquinox() 1", expected, actual);

    in = 55565.9051733796;
    expected = 8.02104092023363e-5;
    actual = EqnEquinox(in);

    double_test("EqnEquinox() 2", expected, actual);

    in = 54332.4868655555;
    expected = 3.23022519949022e-5;
    actual = EqnEquinox(in);

    double_test("EqnEquinox() 3", expected, actual);
}

void test_gmst()
{
    double in;
    double expected;
    double actual;

    in = 54977.6669066321;
    expected = 2.17186902706532;
    actual = gmst(in);

    double_test_delta("gmst() 1", expected, actual, 1e-10);

    in = 53989.1984286448;
    expected = 1.07347347929379;
    actual = gmst(in);

    double_test_delta("gmst() 2", expected, actual, 1e-9);

    in = 55565.9044057253;
    expected = 1.21707647675569;
    actual = gmst(in);

    double_test_delta("gmst() 3", expected, actual, 1e-10);
}

void test_NutMatrix()
{
    double in;
    double **expected = malloc(3 * sizeof(double*));
    double **actual;

    in = 54977.6676696643;
    expected[0] = (double*)(double[3]) {0.999999997896984, -5.95170051004573e-5, -2.58022713429399e-5};
    expected[1] = (double*)(double[3]) {5.95164295625406e-5, 0.999999997980121, -2.23059014984317e-5};
    expected[2] = (double*)(double[3]) {2.58035988712757e-5, 2.23043657924249e-5, 0.999999999418345};

    actual = NutMatrix(in);

    matrix_test("NutMatrix() 1", expected, actual, 3, 3);

    free(actual);

    in = 53989.1991812154;
    expected[0] = (double*)(double[3]) {0.999999999974807, -6.51263906813664e-6, -2.82345706801941e-6};
    expected[1] = (double*)(double[3]) {6.51250710857693e-6, 0.999999998886743, -4.67343651293639e-5};
    expected[2] = (double*)(double[3]) {2.82376142892835e-6, 4.67343467404624e-5, 0.999999998903964};

    actual = NutMatrix(in);

    matrix_test("NutMatrix() 2", expected, actual, 3, 3);

    free(actual);

    in = 55565.9051733796;
    expected[0] = (double*)(double[3]) {0.999999996178561, -8.02104091001633e-5, -3.47730872384987e-5};
    expected[1] = (double*)(double[3]) {8.02104352944664e-5, 0.99999999678286, 7.51898494688596e-7};
    expected[2] = (double*)(double[3]) {3.47730268165429e-5, -7.54687656367992e-07, 0.999999999395134};

    actual = NutMatrix(in);

    matrix_test("NutMatrix() 3", expected, actual, 3, 3);

    free(actual);
    free(expected);
}

void test_PrecMatrix()
{
    double Mjd_1;
    double Mjd_2;
    double **expected = malloc(3 * sizeof(double*));
    double **actual;

    Mjd_1 = 51544.5;
    Mjd_2 = 54977.66766966443;
    expected[0] = (double*)(double[3]) {0.999997373802329, -0.00210194819368335, -0.00091334672251536};
    expected[1] = (double*)(double[3]) {0.00210194819366918, 0.999997790903995, -9.59920515567498e-7};
    expected[2] = (double*)(double[3]) {0.000913346722547957, -9.59889498958352e-7,0.999999582898335};

    actual = PrecMatrix(Mjd_1, Mjd_2);

    matrix_test_delta("PrecMatrix() 1", expected, actual, 3, 3, 1e-11);

    free(actual);

    Mjd_1 = 51544.5;
    Mjd_2 = 53989.1991812154;
    expected[0] = (double*)(double[3]) {0.999998668371322, -0.00149674925171971, -0.000650382395108268};
    expected[1] = (double*)(double[3]) {0.00149674925171607, 0.999998879870093, -4.86735605256216e-7};
    expected[2] = (double*)(double[3]) {0.000650382395116649, -4.86724406102476e-7, 0.999999788501229};

    actual = PrecMatrix(Mjd_1, Mjd_2);

    matrix_test("PrecMatrix() 2", expected, actual, 3, 3);

    free(actual);


    Mjd_1 = 51544.5;
    Mjd_2 = 55565.9051733796;
    expected[0] = (double*)(double[3]) {0.999996396736129, -0.00246210631618179, -0.00106983514929105};
    expected[1] = (double*)(double[3]) {0.00246210631615512, 0.999996969010783, -1.3170512357239e-6};
    expected[2] = (double*)(double[3]) {0.00106983514935241, -1.31700138827348e-6, 0.999999427725346};

    actual = PrecMatrix(Mjd_1, Mjd_2);

    matrix_test_delta("PrecMatrix() 3", expected, actual, 3, 3, 1e-11);

    free(actual);
    free(expected);
}

void test_PoleMatrix()
{
    double xp;
    double yp;
    double **expected = malloc(3 * sizeof(double*));
    double **actual;

    xp = 7.5789200806793e-8;
    yp = 2.56777193042581e-6;
    expected[0] = (double*)(double[3]) {0.999999999999997, 1.94609382460874e-13, 7.57892008065431e-08};
    expected[1] = (double*)(double[3]) {0., 0.999999999996703, -2.56777193042298e-6};
    expected[2] = (double*)(double[3]) {-7.57892008067929e-8, 2.56777193042298e-6, 0.9999999999967};

    actual = PoleMatrix(xp, yp);

    matrix_test("PoleMatrix() 1", expected, actual, 3, 3);

    free(actual);

    xp = 3.41489253490265e-7;
    yp = 1.23201492955798e-6;
    expected[0] = (double*)(double[3]) {0.999999999999942, 4.20719858583501e-13, 3.41489253489999e-7};
    expected[1] = (double*)(double[3]) {0., 0.999999999999241, -1.23201492955767e-6};
    expected[2] = (double*)(double[3]) {-3.41489253490258e-7, 1.23201492955759e-6, 0.999999999999183};

    actual = PoleMatrix(xp, yp);

    matrix_test("PoleMatrix() 2", expected, actual, 3, 3);

    free(actual);

    xp = 5.72229229392844e-7;
    yp = 9.71655423937364e-7;
    expected[0] = (double*)(double[3]) {0.999999999999836, 5.56009634474937e-13, 5.72229229392543e-7};
    expected[1] = (double*)(double[3]) {0., 0.999999999999528, -9.71655423937211e-7};
    expected[2] = (double*)(double[3]) {-5.72229229392813e-7, 9.71655423937052e-7, 0.999999999999364};

    actual = PoleMatrix(xp, yp);

    matrix_test("PoleMatrix() 3", expected, actual, 3, 3);

    free(actual);
    free(expected);
}

void test_newtonnu()
{
    double ecc;
    double nu;
    double e0_e;
    double e0_a = 0.;
    double m_e;
    double m_a =0.;

    ecc = 0.0825331061735532;
    nu = 0.181200311058077;
    e0_e = 0.166883937745775;
    m_e = 0.153174331359043;

    newtonnu(ecc, nu, &e0_a, &m_a);
    double_test("newtonnu() 1, e0", e0_e, e0_a);
    double_test("newtonnu() 1, m", m_e, m_a);

    ecc = 0.0866943121616373;
    nu = 0.171305852271464;
    e0_e = 0.157107124096794;
    m_e = 0.143542791750356;

    newtonnu(ecc, nu, &e0_a, &m_a);
    double_test("newtonnu() 2, e0", e0_e, e0_a);
    double_test("newtonnu() 2, m", m_e, m_a);

    ecc = 0.0791291077778817;
    nu = 0.190267237462258;
    e0_e = 0.175840369293131;
    m_e = 0.161997870558239;

    newtonnu(ecc, nu, &e0_a, &m_a);
    double_test("newtonnu() 3, e0", e0_e, e0_a);
    double_test("newtonnu() 3, m", m_e, m_a);
}


void test_rv2coe()
{
    double *r;
    double *v;

    double p_e;
    double p_a = 0.;
    double a_e;
    double a_a = 0.;
    double ecc_e;
    double ecc_a = 0.;
    double incl_e;
    double incl_a = 0.;
    double omega_e;
    double omega_a = 0.;
    double argp_e;
    double argp_a = 0.;
    double nu_e;
    double nu_a = 0.;
    double m_e;
    double m_a = 0.;
    double arglat_e;
    double arglat_a = 0.;
    double truelon_e;
    double truelon_a = 0.;
    double lonper_e;
    double lonper_a = 0.;

    r = (double*)(double[3]) {20435422.3521528, 1070699.44671798, 1012905.49143388};
    v = (double*)(double[3]) {17.1964697822141, -2657.51027611546, 3738.38685080841};

    p_e = 22151801.03597;
    a_e = 22303727.6412011;
    ecc_e = 0.0825331061735532;
    incl_e = 2.18724254265762;
    omega_e = 0.0874259916432499;
    argp_e = 6.16261219029842;
    nu_e = 0.181200311058077;
    m_e = 0.153174331359043;
    arglat_e = 999999.1;
    truelon_e = 999999.1;
    lonper_e = 999999.1;

    rv2coe(r, v, &p_a, &a_a, &ecc_a, &incl_a, &omega_a, &argp_a, &nu_a, &m_a, &arglat_a, &truelon_a, &lonper_a);
    double_test_delta("rv2coe() 1, p", p_e, p_a, 1e-7);
    double_test_delta("rv2coe() 1, a", a_e, a_a, 1e-7);
    double_test_delta("rv2coe() 1, ecc", ecc_e, ecc_a, 1e-7);
    double_test_delta("rv2coe() 1, incl", incl_e, incl_a, 1e-7);
    double_test_delta("rv2coe() 1, omega", omega_e, omega_a, 1e-7);
    double_test_delta("rv2coe() 1, argp", argp_e, argp_a, 1e-7);
    double_test_delta("rv2coe() 1, nu", nu_e, nu_a, 1e-7);
    double_test_delta("rv2coe() 1, m", m_e, m_a, 1e-7);
    double_test_delta("rv2coe() 1, arglat", arglat_e, arglat_a, 1e-7);
    double_test_delta("rv2coe() 1, truelon", truelon_e, truelon_a, 1e-7);
    double_test_delta("rv2coe() 1, lonper", lonper_e, lonper_a, 1e-7);
}

void test_gibbs()
{
    double *r1;
    double *r2;
    double *r3;

    double *v2_e;
    double **v2_a = malloc(sizeof(double*));
    double theta_e;
    double theta_a = 0.;
    double theta1_e;
    double theta1_a = 0.;
    double copa_e;
    double copa_a = 0.;
    char *error_e = "          ok";
    char **error_a = malloc(sizeof(char*));

    r1 = (double*)(double[3]) {20387627.0717525, 1865163.69633389, -109943.688555792};
    r2 = (double*)(double[3]) {20435422.3521528, 1070699.44671798, 1012905.49143388};
    r3 = (double*)(double[3]) {20398157.066627, 271778.615869933, 2131538.39542067};

    v2_e = (double*)(double[3]) {17.4448460122937, -2659.68695026789, 3741.47770474848};
    theta_e = 0.0672088229314452;
    theta1_e = 0.0670842334897569;
    copa_e = 3.7999117741272e-15;

    gibbs(r1, r2, r3, v2_a, &theta_a, &theta1_a, &copa_a, error_a);
    array_test_delta("gibbs() 1, v2", v2_e, *v2_a, 3, 1e-8);
    double_test("gibbs() 1, theta", theta_e, theta_a);
    double_test("gibbs() 1, theta_1", theta1_e, theta1_a);
    double_test("gibbs() 1, copa", copa_e, copa_a);
    printf("gibbs() 1, error: %s\n", *error_a);

    free(error_a);
    free(v2_a);
}

void test_hgibbs()
{
    double *r1;
    double *r2;
    double *r3;
    double MJD1;
    double MJD2;
    double MJD3;

    double *v2_e;
    double **v2_a = malloc(sizeof(double*));
    double theta_e;
    double theta_a = 0.;
    double theta1_e;
    double theta1_a = 0.;
    double copa_e;
    double copa_a = 0.;
    char *error_e = "          ok";
    char **error_a = malloc(sizeof(char*));

    r1 = (double*)(double[3]) {20387627.0717525, 1865163.69633389, -109943.688555792};
    r2 = (double*)(double[3]) {20435422.3521528, 1070699.44671798, 1012905.49143388};
    r3 = (double*)(double[3]) {20398157.066627, 271778.615869933, 2131538.39542067};
    MJD1 = 55565.9044073611;
    MJD2 = 55565.9078795835;
    MJD3 = 55565.9113518056;

    v2_e = (double*)(double[3]) {17.4309174158407, -2657.49386520115, 3738.39266893676};
    theta_e = 0.0672088229314452;
    theta1_e = 0.0670842334897569;
    copa_e = 3.7999117741272e-15;

    hgibbs(r1, r2, r3, MJD1, MJD2, MJD3, v2_a, &theta_a, &theta1_a, &copa_a, error_a);
    array_test_delta("hgibbs() 1, v2", v2_e, *v2_a, 3, 1e-5);
    double_test_delta("hgibbs() 1, theta", theta_e, theta_a, 1e-5);
    double_test_delta("hgibbs() 1, theta_1", theta1_e, theta1_a, 1e-5);
    double_test_delta("hgibbs() 1, copa", copa_e, copa_a, 1e-5);
    printf("hgibbs() 1, error: %s\n", *error_a);

    free(error_a);
    free(v2_a);
}

void test_IERS()
{
    FILE *fp = fopen("/home/AD.MTHREE.ES/meoberto/preliminary_orbit/eop19620101.txt","r");
    if(!fp)
    {
        fprintf(stderr, "ERROR: No se ha podido abrir eop19620101.txt\n");
        exit(1);
    }
    double **eopdata = malloc(20026 * sizeof(double*));
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

    double Mjd_UTC = 55565.9044073611;
    char interp = 'l';

    double UT1_UTC_e = -0.141329177689613;
    double UT1_UTC_a = 0.;
    double TAI_UTC_e = 34;
    double TAI_UTC_a = 0.;
    double x_pole_e = 5.72229229392844e-7;
    double x_pole_a = 0.;
    double y_pole_e = 9.71655423937364e-7;
    double y_pole_a = 0.;
    double ddpsi_e = -3.22818789921021e-7;
    double ddpsi_a = 0.;
    double ddeps_e = -3.03101990895201e-8;
    double ddeps_a = 0.;

    Mjd_UTC = 55565.9044073611;
    interp = 'l';

    UT1_UTC_e = -0.141329177689613;
    TAI_UTC_e = 34.;
    x_pole_e = 5.72229229392844e-7;
    y_pole_e = 9.71655423937364e-7;
    ddpsi_e = -3.22818789921021e-7;
    ddeps_e = -3.03101990895201e-8;

    IERS(eopdata, 20026, Mjd_UTC, interp, &UT1_UTC_a, &TAI_UTC_a, &x_pole_a, &y_pole_a, &ddpsi_a, &ddeps_a);

    double_test_delta("IERS() 1, UT1_UTC", UT1_UTC_e, UT1_UTC_a, 1e-8);
    double_test_delta("IERS() 1, TAI_UTC", TAI_UTC_e, TAI_UTC_a, 1e-8);
    double_test_delta("IERS() 1, x_pole", x_pole_e, x_pole_a, 1e-8);
    double_test_delta("IERS() 1, y_pole", y_pole_e, y_pole_a, 1e-8);
    double_test_delta("IERS() 1, ddpsi", ddpsi_e, ddpsi_a, 1e-8);
    double_test_delta("IERS() 1, ddeps", ddeps_e, ddeps_a, 1e-8);

    for(int i=0; i<20026; i++)
            free(eopdata[i]);
    free(eopdata);
}


int main()
{
    /*MatlabUtilsTest();
    test_Position();
    test_Mjday();
    test_MeanObliquity();
    test_NutAngles();
    test_timediff();
    test_Frac();
    test_R_x();
    test_R_y();
    test_R_z();
    test_unit();
    test_EqnEquinox();
    test_gmst();
    test_NutMatrix(); */
    test_IERS();
    //test_gast();
    //test_GHAMatrix();
    /*test_PrecMatrix();
    test_PoleMatrix();
    test_newtonnu();
    test_rv2coe();
    test_gibbs();
    test_hgibbs();*/
    return 0;
}
