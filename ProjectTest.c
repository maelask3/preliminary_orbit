#include "TestUtils.h"
#include "MatlabUtilsTest.h"
#include "Position.h"
#include "Mjday.h"
#include <stdlib.h>

void test_Position()
{
    double lon = -2.117969613666;
    double lat = 0.683053277791;
    double h = 99.81638;

    double *expected = (double*)(double[3]) {-2577383.6395731382, -4230610.2467987221, 4004108.3320587045};
    double *actual = Position(lon, lat, h);
    array_test_delta("Position() 1", expected, actual, 3, 10e-6);
    free(actual);

    lon = -1.504723397302;
    lat = 0.533589040237;
    h = 0.;

    expected = (double*)(double[3]) {362889.51475075335,-5484262.3610134749,3225167.7284776145};
    actual = Position(lon, lat, h);
    array_test_delta("Position() 2", expected, actual, 3, 10e-6);
    free(actual);

    lon = -1.919862177194;
    lat = 0.698131700798;
    h = 2000.;

    expected = (double*)(double[3]) {-1673928.5598879098,-4599080.9200722268,4079271.1474197493};
    actual = Position(lon, lat, h);
    array_test_delta("Position() 7", expected, actual, 3, 10e-6);
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
int main()
{
    test_Position();
    test_Mjday();
    return MatlabUtilsTest();
}
