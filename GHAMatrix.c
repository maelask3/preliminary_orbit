#include "GHAMatrix.h"
#include "R_z.h"
#include "gast.h"

double **GHAMatrix(double Mjd_UT1)
{
    double **GHAmat = R_z(gast(Mjd_UT1));
    return GHAmat;
}
