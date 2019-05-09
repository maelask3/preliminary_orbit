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

void matrix_test_delta(char *test_name, double **expected, double **actual, int rows, int cols, double delta)
{
    printf("================================================================================\n");
    printf("Now running test %s:\n",test_name);
    printf("Actual = \n");
    matrix_pretty_print(actual, rows, cols);
    double diff = 0.;
    for(int i=0; i<rows; i++)
        for(int j=0; j<cols; j++)
        {
            diff = fabs(expected[i][j] - actual[i][j]);
            if(diff > delta)
                fprintf(stderr, "Divergence = %lf\n", diff);
            assert(fabs(expected[i][j] - actual[i][j]) < delta);
        }
    printf("Expected = \n");
    matrix_pretty_print(expected, rows, cols);
    printf("PASS\n");
}

void matrix_test(char *test_name, double **expected, double **actual, int rows, int cols)
{
    matrix_test_delta(test_name, expected, actual, rows, cols, DELTA);
}

void double_test_delta(char *test_name, double expected, double actual, double delta)
{
    printf("================================================================================\n");
    printf("Now running test %s:\n",test_name);
    FILE *os;
    double diff = fabs(expected - actual);
    if(diff > delta)
        os = stderr;
    else
        os = stdout;
    fprintf(os, "Actual = %lf\n", actual);
    fprintf(os, "Expected = %lf\n", expected);
    fprintf(os,"Divergence = %lf\n",diff);
    assert(diff < delta);
    printf("PASS\n");
}

void double_test(char *test_name, double expected, double actual)
{
    double_test_delta(test_name, expected, actual, DELTA);
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
