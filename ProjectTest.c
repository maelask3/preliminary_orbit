#include "TestUtils.h"
#include "MatlabUtilsTest.h"
#include "Position.h"
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
}
int main()
{
    test_Position();
    return MatlabUtilsTest();
}
