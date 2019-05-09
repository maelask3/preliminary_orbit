/**
 * @file angl.c
 * @Autor Davide Pérez y Millán Santamaría
 */

#include "MatlabUtils.h"
#include <stdio.h>
#include <math.h>

double angl(double *vec1, double *vec2){
	double small = 0.00000001;
	double undefined = 999999.1;

	double magv1 = norm(vec1);
	double magv2 = norm(vec2);
	
	double theta;
    if(magv1*magv2 > (small*small))
    {
		double temp = dot(vec1, vec2)/(magv1*magv2);
		if(fabs(temp) > 1){
			temp = sign(temp);
		}
		theta = acos(temp);
	}else{
		theta = undefined;
	}

    return theta;
}
