/**
 * @file MatlabUtils.c
 * @author Davide Pérez y Millán Santamaría
 * @brief Es una libreria de funciones presentes en Matlab y no en C.
 */
#include "MatlabUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "rpoly.h"

/**
 * @brief Calcula la norma de un vector de dimensión 3
 * @param v es un vector de doubles de dimensión 3
 * @return Devuelve la norma de v
 */
double norm(double *v)
{
        return(sqrt((v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2])));
}

/**
 * @brief Calcula el producto escalar de dos vectores de tamaño 3
 * @param v es un vector de doubles de tamaño 3
 * @param t es un vector de doubles de tamaño 3
 * @return Devuelve el producto escalar de v y t
 */
double dot(double *v, double *t)
{
	return(v[0]*t[0]+v[1]*t[1]+v[2]*t[2]);
}

/**
 * @brief Determina el signo de un número.
 * @param n es un número de tipo double con signo.
 * @return Devuelve 0.0 si n es cero, -1.0 sin n es menor que cero y 1.0 si n es mayor que cero.
 */
double sign(double n)
{
    if(fabs(n) < 10e-12)
		return 0.0;
	else
		return ((n>0)? 1.0 : -1.0);
}

/**
 * @brief Inicializa matrices de tamano rows x cols.
 * @param rows es un int sin signo que refiere al número de filas.
 * @param cols es u int sin signo que refiere al número de columnas.
 * @return Devuelve una matriz de ceros de tamano rows x cols.
 */
double **zeros(unsigned int rows, unsigned int cols)
{
	double **matrix;
    matrix = malloc (rows*sizeof(double *));
    for(unsigned int i=0; i<rows; i++)
	{
            matrix[i] = calloc(cols, sizeof(double));
  	}
	return(matrix);
}

/**
 * @brief Calcula el determinante de una matrix 3x3.
 * @param m es una matriz de doubles de tamaño 3x3.
 * @return Devuelve el determinante de m.
 */
double det(double **m)
{
	return m[0][0]*m[1][1]*m[2][2]+m[0][1]*m[1][2]*m[2][0]+m[0][2]*m[1][0]*m[2][1]-m[2][0]*m[1][1]*m[0][2]-m[2][1]*m[1][2]*m[0][0]-m[2][2]*m[1][0]*m[0][1];
}

/**
 * @brief Calcula el determinante de una matrix 2x2.
 * @param m es una matriz de doubles de tamaño 2x2.
 * @return Devuelve el determinante de m.
 */
double det2x2(double **m)
{
	return m[0][0]*m[1][1]-m[1][0]*m[0][1];
}

/**
 * @brief Calcula las raices reales de un polinomio.
 * @param coef es un vector de doubles de tamaño 16, sus componentes son los coeficientes ordenamos de de mayor a menor grado
 * @param sols_reales es un vector de doubles vacio de tamaño 15 . Es un parametro de entrada/salida, en el se devuelven las raices encontradas.
 * @return Devuelve un int que es el número de raices encontradas, es decir hasta que indice de sols_reales debemos iterar.
 */
int roots(double *coef, double **sols_reales)
{
	double *real, *imagin, *aux;
    real = malloc(15*sizeof(double));
    imagin = malloc(15*sizeof(double));
    int solutions = real_poly_roots(coef, 15, real, imagin);
	
	if(solutions < 0)
		return solutions;

    aux = calloc(solutions, sizeof(double));
    int j = 0;
    for(int i = 0; i< solutions; i++)
    {

            if((fabs(imagin[i]) < 10e-12) && (fabs(real[i]) > 10e-12))
            {
                    aux[j] = real[i];
                    j++;
            }
    }
    *sols_reales = aux;

    free(real);
    free(imagin);
    return j;
}

/**
 * @brief Redondea un número al entero más cercano a 0
 * @param in Número
 * @return Redondea hacia arriba si in es negativo, y hacia abajo si in es positivo.
 */

double fix(double in)
{
    return (in>0) ? floor(in) : ceil(in);
}

/**
 * @brief Devuelve el módulo de la división entre a y m
 * @param a Dividendo
 * @param m Divisor
 * @return Devuelve el módulo b tal que b = a - m*floor(a/m)
 */
double mod(double a, double m)
{
    return a - m*floor(a/m);
}

/**
 * @brief Calcula el producto cruzado de dos vectores de tamaño 3 .
 * @param v1 es un vector de doubles de tamaño 3 .
 * @param v2 es un vector de doubles de tamaño 3 .
 * @return Devuelve el producto cruzado de v1 y v2.
 */
double *cross(double *v1, double *v2)
{
    double *res = calloc(3, sizeof(double));

    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return res;
}

//O(n^3), fuente: www.kkhaydarov.com/matrux-multiplication-algorithms/, naive matrix
/**
 * @brief Calcula el producto de dos matrices de tamaño 3x3.
 * @param m1 es una matriz de doubles de tamaño 3x3.
 * @param m2 es una matriz de doubles de tamaño 3x3.
 * @return Devuelve el producto de m1 y m2, una matrix 3x3 de doubles.
 */
double **productMatrix(double **m1, double **m2){
	double **product = zeros(3, 3);
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			product[i][j] = 0;
			for(int k=0; k<3; k++)
			{
				product[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}

	return product;
}

/**
 * @brief Calcula la suma de dos matrices 3x3.
 * @param m1 es una matriz de doubles de tamaño 3x3.
 * @param m2 es una matriz de doubles de tamaño 3x3.
 * @return Devuelve la suma de m1 y m2, una matrix 3x3 de doubles.
 */
double **sumMatrix(double **m1, double **m2){
	double **sum = zeros(3, 3);
	for(int i = 0; i<3; i++){
		for(int j = 0; j<3; j++){
			sum[i][j] = m1[i][j] + m2[i][j];
		}
	}
	return sum;
}

/**
 * @brief Calcula la matriz transpuesta de una matriz de tamaño 3x3.
 * @param m es una matriz de doubles de tamaño 3x3.
 * @return Devuelve la matriz transpuesta de m.
 */
double **transposeMatrix(double **m){
	double **transpose = zeros(3, 3);
	for(int i = 0; i<3; i++){
		for(int j = 0; j<3; j++){
			transpose[i][j] = m[j][i];
		}
	}
	return transpose;
}

/**
 * @brief Calcula el producto de un vector de tamaño 3 por un double
 * @param v es un vector de tamaño 3
 * @param d es un double
 * @return Devuelve el producto del vector por el double.
 */

double *vectorProductDouble(double *v, double d)
{
	double *res = calloc(3, sizeof(double));
	for(int i=0; i<3; i++)
	{
        res[i] = d * v[i];
	}

	return res;
}

/**
 * @brief Calcula la suma de dos vectores de tamaño 3
 * @param v1 es un vector de tamaño 3
 * @param v1 es un vector de tamaño 3
 * @return Devuelve la suma de los dos vectores.
 */

double *sumVector(double *v1, double *v2)
{
	double *res = calloc(3, sizeof(double));
	for(int i=0; i<3; i++)
	{
		res[i] = v1[i] + v2[i];
	}

	return res;
}
