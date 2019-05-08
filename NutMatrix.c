#include "NutMatrix.h"
#include "R_x.h"
#include "R_z.h"
#include "NutAngles.h"
#include "MeanObliquity.h"
#include "MatlabUtils.h"

double **NutMatrix(double Mjd_TT)
{
    double ep = MeanObliquity(Mjd_TT);
    double dpsi = 0.;
    double deps = 0.;
    NutAngles(Mjd_TT, &dpsi, &deps);

    double **NutMat = productMatrix(R_x(-ep-deps),productMatrix(R_z(-dpsi), R_x(ep)));

    return NutMat;
}
