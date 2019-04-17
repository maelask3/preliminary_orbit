#include <assert.h>
#include "MatlabUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

void array_pretty_print(double *v)
{
	printf("[%lf, %lf, %lf]", v[0], v[1], v[2]);
}

void test_norm()
{
	printf("TESTING: norm()\n");
	double v[3] = {1., 1., 1.};
	assert(fabs(norm(v) - sqrt(3.)) < 10e-12);
	double w[3] = {1., 2., 3.};
	assert(fabs(norm(w) - sqrt(14.)) < 10e-12);
}

void test_dot()
{
	printf("TESTING: dot()\n");
	double v[3] = {1., 5., 8.};
	double w[3] = {5., 9., 1.};
	assert(fabs(dot(v, w) - 58.) < 10e-12);
}

void test_sign()
{
	printf("TESTING: sign()\n");
	assert(fabs(sign(45) - 1.) < 10e-12);
	assert(sign(0) == 0.);
	assert(sign(-666) == -1.);
}

void test_zeros()
{

}

void test_det2x2()
{
	printf("TESTING: det2x2()\n");
	double m[][2] = {{1., 2.}, {3., 4.}};
	
	assert(fabs(det2x2(m) - -2.) < 10e-12);
	
	double n[][2] = {{5., 9.}, {74., 25.}};
	assert(fabs(det2x2(n) - -541.) < 10e-12);
}

void test_det()
{
	printf("TESTING: det()\n");
	double m[][3] = {{1.,2.,3.},{4.,5.,6.},{7.,8.,9.}};
	assert(det(m) == 0.);
	double n[][3] = {{9.,1.,1.},{1.,1.,2.},{5.,4.,3.}};
	assert(fabs(det(n) - -39.) < 10e-12);
}

void test_roots()
{
}

void test_isreal()
{
	printf("TESTING: isreal()\n");
	double complex z = -I;
	double complex q = 25;
	double r = 88;

	assert(isreal(z) == false);
	assert(isreal(q) == true);
	assert(isreal(r) == true);
}

void test_cross()
{
	printf("TESTING: cross()\n");
	double v[3] = {1.,5.,8.};
	double w[3] = {5.,9.,1.};
	printf("Vector v: ");
	array_pretty_print(v);
	printf("\nVector w: ");
	array_pretty_print(w);
	printf("\nProducto vectorial: ");
	double vwc[3];
       	memcpy(vwc, cross(v, w), sizeof(double[3]));
	array_pretty_print(vwc);
	printf("\n");
	assert(vwc[0] == -67.0);
	assert(vwc[1] == 39.0);
	assert(vwc[2] == -16.0);

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
	test_isreal();
	test_cross();
}
