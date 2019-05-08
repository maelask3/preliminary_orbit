#include "unit.h"
#include "MatlabUtils.h"

double *unit(double *vec)
{
    double small = 0.000001;
    double magv = norm(vec);

    double *outvec = calloc(3, sizeof(double));
    if(magv > small)
    {
        for(int i=0; i<3; i++)
            outvec[i] = vec[i]/magv;
    }

    return outvec;
}
