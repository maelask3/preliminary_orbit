#include "TestUtils.h"
#include "MatlabUtilsTest.h"
#include "Position.h"
#include "Mjday.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include <stdlib.h>

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
int main()
{
    test_Position();
    test_Mjday();
    test_MeanObliquity();
    test_NutAngles();
    return MatlabUtilsTest();
}
