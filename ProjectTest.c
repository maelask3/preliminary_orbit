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
#include <stdlib.h>
#include <stdio.h>
#include "lambert_gooding.h"

void test_Position()
{
    double lon = -2.11796961366573;
    double lat = 0.683053277790977;
    double h = 99.81638;

    double *expected = (double*)(double[3]) {-2577383.6395731382, -4230610.2467987221, 4004108.3320587045};
    double *actual = Position(lon, lat, h);
    array_test_delta("Position() 1", expected, actual, 3, 10e-7);
    free(actual);

    lon = -1.50472339730215;
    lat = 0.533589040236714;
    h = 0.;

    expected = (double*)(double[3]) {362889.51475075335,-5484262.3610134749,3225167.7284776145};
    actual = Position(lon, lat, h);
    array_test_delta("Position() 2", expected, actual, 3, 10e-7);
    free(actual);

    lon = -1.91986217719376;
    lat = 0.698131700797732;
    h = 2000.;

    expected = (double*)(double[3]) {-1673928.5598879098,-4599080.9200722268,4079271.1474197493};
    actual = Position(lon, lat, h);
    array_test_delta("Position() 7", expected, actual, 3, 10e-7);
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
    expected[1] = (double*)(double[3]) {0., 0.917481675640187, -0.397778047237997};
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

int main()
{
    MatlabUtilsTest();
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
    test_NutMatrix();

    //Test lambert_gooding
/*
	double v1[3], v2[3];
	double r1[] = {40706177.3896682, 9894896.081818309, -235154.0576340267};
	double r2[] = {40481024.54943731, 10782719.6281578, -234124.5674719382};
	double tof = 300.0000223517418;
	double mu = 398600441800000.0;
	double long_way = 0.0;
	double multi_revs = 1.0;

	double sol_v1[] = {-717.4652239083363, 2967.699531994756, 3.240668317805917};
	double sol_v2[] = {-783.4918306318922, 2950.883185650744, 3.622315600677195};

	lambert_gooding(r1, r2, tof, mu, long_way, multi_revs, v1, v2);

	printf("Actual  -------------------  Esperado \n");
	for(int i = 0; i<3; i++)
	{
		printf("%lf  ---------v1----------  %lf \n", v1[i], sol_v1[i]);
		printf("%lf  ---------v2----------  %lf \n", v2[i], sol_v2[i]);
	}*/
	//Test vlamb
	/*double gm = 398600441800000.0;
	double r1 = 41892208.62618656;
	double r2 = 41893140.31514826;
	double th = 0.02186412186035686;
	double tdelt = 300.0000223517418;
	double vri[2], vti[2], vrf[2], vtf[2];

	double sol_n = 1.0;
	double sol_vri[] = {3.796642495817074, 0.0};
	double sol_vti[] = {3053.19389282054, 0.0};
	double sol_vrf[] = {2.414379595078204, 0.0};
	double sol_vtf[] = {3053.125990843594, 0.0};

	double n = vlamb(gm, r1, r2, th, tdelt, vri, vti, vrf, vtf);

	printf("Actual  -------------------  Esperado \n");
	printf("%lf  ---------n----------  %lf \n", n, sol_n);
    	printf("%lf  ---------vri----------  %lf \n", vri[0], sol_vri[0]);
    	printf("%lf  ---------vti----------  %lf \n", vti[0], sol_vti[0]);
    	printf("%lf  ---------vrf----------  %lf \n", vrf[0], sol_vrf[0]);
    	printf("%lf  ---------vtf----------  %lf \n", vtf[0], sol_vtf[0]);*/

    //Test tlamb
/*
	double m = 0.0;
    	double q = 0.9891272559119239;
    	double qsqfm1 = 0.02162727161214778;
    	double x = 0.0;
    	double n = 0.0;

    	double sol[] = {0.5861212399757942, 0.0, 0.0, 0.0};

    	double *aux = tlamb(m, q, qsqfm1, x, n);

    	printf("Actual  -------------------  Esperado \n");
    	printf("%lf  ---------t----------  %lf \n", aux[0], sol[0]);
    	printf("%lf  ---------dt----------  %lf \n", aux[1], sol[1]);
    	printf("%lf  ---------d2t----------  %lf \n", aux[2], sol[2]);
    	printf("%lf  ---------d2t----------  %lf \n", aux[3], sol[3]);*/
    //Test xlamb
    	double m = 0.0;
	double q = 0.9891272559119239;
	double qsqfm1 = 0.02162727161214778;
	double tin = 0.06146745322679799;

	double sol[] = {1.0, 0.6959029065435739, 0.0};
	double *aux = xlamb(m, q, qsqfm1, tin);
	printf("Actual  -------------------  Esperado \n");
        printf("%lf  ---------n----------  %lf \n", aux[0], sol[0]);
        printf("%lf  ---------x----------  %lf \n", aux[1], sol[1]);
        printf("%lf  ---------xpl----------  %lf \n", aux[2], sol[2]);

    return 0;
}
