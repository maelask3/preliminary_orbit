#ifndef TESTUTILS_H
#define TESTUTILS_H
#define DELTA 10e-12

void array_pretty_print(double *v, int sz);
void matrix_pretty_print(double **m, int rows, int cols);
void matrix_test(char *test_name, double **expected, double **actual, int rows, int cols);
void double_test(char *test_name, double expected, double actual);
void array_test_delta(char *test_name, double *expected, double *actual, int sz, double delta);
void array_test(char *test_name, double *expected, double *actual, int sz);

#endif // TESTUTILS_H
