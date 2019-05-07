#include <assert.h>
#include "MatlabUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#define DELTA 10^(-12)

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
	printf("Actual = ");
	matrix_pretty_print(actual, rows, cols);
	for(int i=0; i<rows; i++)
        for(int j=0; j<cols; j++)
			assert(fabs(expected[i][j] - actual[i][j]) < DELTA);
	printf("Expected = ");
	matrix_pretty_print(expected, rows, cols);
}

void double_test(char *test_name, double expected, double actual)
{
	printf("================================================================================\n");
	printf("Now running test %s:\n",test_name);
	printf("Actual = %lf\n", actual);
	assert(fabs(expected - actual) < DELTA);
	printf("Expected = %lf\n", expected);
}

void array_test(char *test_name, double *expected, double *actual, int sz)
{
	printf("================================================================================\n");
	printf("Now running test %s:\n",test_name);
	printf("Actual = ");
	array_pretty_print(actual, sz);
	for(int i=0; i<sz; i++){
		assert(fabs(actual[i] - expected[i]) < DELTA);
	}
	printf("Expected = ");
	array_pretty_print(expected, sz);
}

void test_norm()
{
	double v[3] = {1., 1., 1.};
	double_test("norm() 1", sqrt(3.), norm(v));
	double w[3] = {1., 2., 3.};
	double_test("norm() 2", sqrt(14.), norm(w));
}

void test_dot()
{
	double v[3] = {1., 5., 8.};
	double w[3] = {5., 9., 1.};
	double_test("dot()", 58., dot(v, w));
}

void test_sign()
{
	double_test("sign() 1", 1., sign(45));
	double_test("sign() 2", 0., sign(0));
	double_test("sign() 3", -1., sign(-666));
}

void test_zeros()
{
	double **v = (double**) calloc(3, sizeof(double*));
	for(int i=0; i<3; i++)
	{
		v[i] = (double*) calloc(3, sizeof(double));
	}
	double **m = zeros(3, 3);
	matrix_test("zeros()", v, m, 3, 3);
}

void test_det2x2()
{
	double m[][2] = {{1., 2.}, {3., 4.}};
	double_test("det2x2() 1", -2., det2x2(m));
	double n[][2] = {{5., 9.}, {74., 25.}};
	double_test("det2x2() 2", -541., det2x2(n));
}

void test_det()
{
	double m[][3] = {{1.,2.,3.},{4.,5.,6.},{7.,8.,9.}};
	double_test("det() 1", 0., det(m));
	double n[][3] = {{9.,1.,1.},{1.,1.,2.},{5.,4.,3.}};
	double_test("det() 2", -39., det(n));
}

void test_roots()
{
	double m[] = {1., 0., -73120740632072.03, 0., 0., -1.587936795676375e+36, 0., 0., -1.198538485369091e+58, 0., 0., 0., 0., 0., 0., 0.};
	double *v;
	v = (double*) malloc(15*sizeof(double));
	int num = roots(m, v);
	double sol[] = {0., 0., 0., 0., 0., 0., 0., 20488505.59583733, -16734286.96763425};
	array_test("roots()", sol, v, num);

}

void test_cross()
{
	double v[3] = {1.,5.,8.};
	double w[3] = {5.,9.,1.};
	double vwce[3] = {-67., 39., -16.};
	double vwca[3];
       	memcpy(vwca, cross(v, w), sizeof(double[3]));
	array_test("cross()", vwce, vwca, 3);
}

void test_productMatrix(){
        double m[][3] = {{1.,2.,3.},{4.,5.,6.},{7.,8.,9.}};
        double n[][3] = {{9.,1.,1.},{1.,1.,2.},{5.,4.,3.}};
	double sol[][3] = {{26.,15.,14.}, {71.,33.,32.},{116.,51.,50.}};
        matrix_test("productMatrix()", sol, productMatrix(m,n), 3,3);
}

int main()
{
	test_norm();
	test_dot();
	test_sign();
	test_zeros();
	test_det2x2();
	test_det();
	test_roots();
	test_cross();
//	test_productMatrix();
	
	return 0;
}
