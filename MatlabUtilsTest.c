/**
 * @file MatlabUtilsTest.c
 * @author Davide Pérez y Millán Santamaría
 * @brief Tests unitarios sobre MatlabUtils.c
 */
#include <assert.h>
#include "MatlabUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "TestUtils.h"

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

    for(int i=0; i<3; i++)
    {
        free(m[i]);
    }
    free(m);

    for(int i=0; i<3; i++)
    {
        free(v[i]);
    }
    free(v);
}

void test_det2x2()
{
    double **m = malloc(2 * sizeof(double*));
    m[0] = (double*)(double[2]){1., 2.};
    m[1] = (double*)(double[2]){3., 4.};
    double_test("det2x2() 1", -2., det2x2(m));
    free(m);

    double **n = malloc(2 * sizeof(double*));
    n[0] = (double*)(double[2]){5., 9.};
    n[1] = (double*)(double[2]){74., 25.};
    double_test("det2x2() 2", -541., det2x2(n));
    free(n);
}

void test_det()
{
    double **mm = malloc(3 * sizeof(double*));
    mm[0] = (double*)(double[3]){1.,2.,3.};
    mm[1] = (double*)(double[3]){4.,5.,6.};
    mm[2] = (double*)(double[3]){7.,8.,9.};
    double_test("det() 1", 0., det(mm));
    free(mm);

    double **nn = malloc(3 * sizeof(double*));
    nn[0] = (double*)(double[3]) {9.,1.,1.};
    nn[1] = (double*)(double[3]) {1.,1.,2.};
    nn[2] = (double*)(double[3]) {5.,4.,3.};
    double_test("det() 2", -39., det(nn));
    free(nn);
}

void test_roots()
{
    double *v = calloc(16, sizeof(double));
    double *m = (double*)(double[16]) {1.0, 0.0, -73120740632127.34375, 0., 0., -1587936795677189147685214247486226432.0, 0., 0.,  -11985384853690594734217583339479868539727097108863410765824.0, 0., 0., 0., 0., 0., 0., 0.};
    int num = roots(m, &v);
    printf("%d\n", num);
    double sol[] = {20488505.5958389, -16734286.9676338};
    array_test_delta("roots()", sol, v, num, 10e-7);
    free(v);
}

void test_cross()
{
    double v[3] = {1.,5.,8.};
    double w[3] = {5.,9.,1.};
    double vwce[3] = {-67., 39., -16.};
    double *vwcap = cross(v, w);
    double vwca[3];
    memcpy(vwca, vwcap, sizeof(double[3]));
    free(vwcap);
    array_test("cross()", vwce, vwca, 3);
}

void test_producto()
{
    double **mm = malloc(3 * sizeof(double*));
    mm[0] = (double*)(double[3]){1.,2.,3.};
    mm[1] = (double*)(double[3]){4.,5.,6.};
    mm[2] = (double*)(double[3]){7.,8.,9.};
    double **nn = malloc(3 * sizeof(double*));
    nn[0] = (double*)(double[3]) {9.,1.,1.};
    nn[1] = (double*)(double[3]) {1.,1.,2.};
    nn[2] = (double*)(double[3]) {5.,4.,3.};
    double **pp = productMatrix(mm, nn);
    free(mm);
    free(nn);

    double **sol = malloc(3 * sizeof(double*));
    sol[0] = (double*)(double[3]) {26., 15., 14.};
    sol[1] = (double*)(double[3]) {71., 33., 32.};
    sol[2] = (double*)(double[3]) {116., 51., 50.};
    matrix_test("productMatrix()", sol, pp, 3, 3);
    free(pp);
    free(sol);
}

void test_suma()
{
    double **mm = malloc(3 * sizeof(double*));
    mm[0] = (double*)(double[3]){1.,2.,3.};
    mm[1] = (double*)(double[3]){4.,5.,6.};
    mm[2] = (double*)(double[3]){7.,8.,9.};
    double **nn = malloc(3 * sizeof(double*));
    nn[0] = (double*)(double[3]) {9.,1.,1.};
    nn[1] = (double*)(double[3]) {1.,1.,2.};
    nn[2] = (double*)(double[3]) {5.,4.,3.};
    double **pp = sumMatrix(mm, nn);
    free(mm);
    free(nn);

    double **sol = malloc(3 * sizeof(double*));
    sol[0] = (double*)(double[3]) {10., 3., 4.};
    sol[1] = (double*)(double[3]) {5., 6., 8.};
    sol[2] = (double*)(double[3]) {12., 12., 12.};

    matrix_test("sumMatrix()", sol, pp, 3, 3);
    free(sol);
    free(pp);
}

void test_transpuesta()
{
    double **mm = malloc(3 * sizeof(double*));
    mm[0] = (double*)(double[3]){1.,2.,3.};
    mm[1] = (double*)(double[3]){4.,5.,6.};
    mm[2] = (double*)(double[3]){7.,8.,9.};
    double **tt = transposeMatrix(mm);
    free(mm);

    double **sol = malloc(3 * sizeof(double*));
    sol[0] = (double*)(double[3]) {1., 4., 7.};
    sol[1] = (double*)(double[3]) {2., 5., 8.};
    sol[2] = (double*)(double[3]) {3., 6., 9.};
    matrix_test("transposeMatrix()", sol, tt, 3, 3);
    free(sol);
    free(tt);
}

void test_fix()
{
    double a = -1.;
    double b = -69.420;
    double c = 13.69;
    double d = 2.;

    double_test("fix() 1", fix(a), -1.);
    double_test("fix() 2", fix(b), -69.);
    double_test("fix() 3", fix(c), 13.);
    double_test("fix() 4", fix(d), 2.);
}
int MatlabUtilsTest()
{
    test_norm();
    test_dot();
    test_sign();
    test_zeros();
    test_det2x2();
    test_det();
    test_roots();
    test_cross();
    test_producto();
    test_suma();
    test_transpuesta();
    test_fix();
    return 0;
}
