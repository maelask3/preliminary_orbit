#include "hgibbs.h"
#include "SAT_Const.h"
#include "MatlabUtils.h"
#include "unit.h"
#include "angl.h"
#include <math.h>

void hgibbs(double *r1, double *r2, double *r3, double MJD1, double MJD2, double MJD3, double **v2, double *theta, double *theta1, double *copa, char **error)
{
    *error =  "          ok";

    *theta = 0.;
    *theta1 = 0.;
    double magr1 = norm(r1);
    double magr2 = norm(r2);
    double magr3 = norm(r3);

    double tolangle = 0.01745329251994;
    double dt21 = (MJD2-MJD1)*86400.;
    double dt31 = (MJD3-MJD1)*86400.;
    double dt32 = (MJD3-MJD2)*86400.;

    double *p = cross(r2, r3);
    double *pn = unit(p);
    double *r1n = unit(r1);
    *copa = asin(dot(pn, r1n));

    if(fabs(dot(r1n, pn)) > 0.017452406)
        *error = "not coplanar";

    *theta = angl(r1, r2);
    *theta1 = angl(r2, r3);

    if((*theta > tolangle) || (*theta1 > tolangle))
        *error = "   angl > 1�";

    double term1 = -dt32*((1/(dt21*dt31)) + GM_Earth/(12*magr1*magr1*magr1));
    double term2 = (dt32-dt21)*(1/(dt21*dt32) + GM_Earth/(12*magr2*magr2*magr2));
    double term3 = dt21*(1/(dt32*dt31) + GM_Earth/(12*magr3*magr3*magr3));

    *v2 = sumVector(vectorProductDouble(r1, term1), sumVector(vectorProductDouble(r2, term2), vectorProductDouble(r3, term3)));
}
