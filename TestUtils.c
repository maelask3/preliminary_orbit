#include "TestUtils.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

void array_pretty_print(double *v, int sz)
{
    printf("[");
    for(int i=0; i<sz; i++)
    {
        printf("%lf", v[i]);
        if(i != sz-1)
        {
            printf(", ");
        } else {
            printf("]\n");
        }
    }

}

void matrix_pretty_print(double **m, int rows, int cols)
{
    for(int i=0; i<rows; i++)
        array_pretty_print(m[i], cols);
}

void matrix_test(char *test_name, double **expected, double **actual, int rows, int cols)
{
    printf("================================================================================\n");
    printf("Now running test %s:\n",test_name);
    printf("Actual = \n");
    matrix_pretty_print(actual, rows, cols);
    for(int i=0; i<rows; i++)
        for(int j=0; j<cols; j++)
            assert(fabs(expected[i][j] - actual[i][j]) < DELTA);
    printf("Expected = \n");
    matrix_pretty_print(expected, rows, cols);
    printf("PASS\n");
}

void double_test(char *test_name, double expected, double actual)
{
    printf("================================================================================\n");
    printf("Now running test %s:\n",test_name);
    printf("Actual = %lf\n", actual);
    assert(fabs(expected - actual) < DELTA);
    printf("Expected = %lf\n", expected);
    printf("PASS\n");
}

void array_test_delta(char *test_name, double *expected, double *actual, int sz, double delta)
{
    printf("================================================================================\n");
    printf("Now running test %s:\n",test_name);
    printf("Actual = ");
    array_pretty_print(actual, sz);
    printf("Expected = ");
    array_pretty_print(expected, sz);
    for(int i=0; i<sz; i++)
    {
        double diff = fabs(expected[i] - actual[i]);
        if(diff > delta)
            fprintf(stderr,"Divergence = %lf\n",diff);
        assert(fabs(expected[i] - actual[i]) < delta);
    }
    printf("PASS\n");
}

void array_test(char *test_name, double *expected, double *actual, int sz)
{
    array_test_delta(test_name, expected, actual, sz, DELTA);
}
