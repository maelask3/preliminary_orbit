#include "R_x.h"
#include "R_y.h"
#include "PoleMatrix.h"
#include "MatlabUtils.h"

double **PoleMatrix(double xp, double yp)
{
    double **PoleMat = productMatrix(R_y(-xp), R_x(-yp));

    return PoleMat;
}
