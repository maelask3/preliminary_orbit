#include "EqnEquinox.h"
#include "NutAngles.h"
#include "MeanObliquity.h"
#include <math.h>

double EqnEquinox(double Mjd_TT)
{
    double dpsi = 0.;
    double deps = 0.;
    NutAngles(Mjd_TT, &dpsi, &deps);

    return dpsi * cos(MeanObliquity(Mjd_TT));
}
