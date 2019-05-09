#include "newtonnu.h"
#include "SAT_Const.h"
#include <math.h>
void newtonnu(double ecc, double nu, double *e0, double *m)
{
    double e0_aux = 999999.9;
    double m_aux = 999999.9;
    double small = 0.00000001;
    double sine = 0.;
    double cose = 0.;

    if(fabs(ecc) < small)
    {
        m_aux = nu;
        e0_aux = nu;
    } else {
        if(ecc < (1.0-small))
        {
            sine = ((sqrt(1.0-(ecc*ecc)) * sin(nu)) / (1.0+(ecc*cos(nu))));
            cose = (ecc + cos(nu)) / (1.0+ (ecc*cos(nu)));
            e0_aux = atan2(sine, cose);
            m_aux = e0_aux - (ecc * sin(e0_aux));
        } else {
            if(ecc > (1.0+small))
            {
                if((ecc > 1.0) && (fabs(nu)+0.00001 < (pi-acos(1.0/ecc))))
                {
                    sine = ((sqrt(ecc*ecc-1.0) * sin(nu)) / (1.0 + (ecc*cos(nu))));
                    e0_aux = asinh(sine);
                    m_aux = ecc*sinh(e0_aux) - e0_aux;
                }
            } else {
                if(fabs(nu) < (168.0*pi/180.0))
                {
                    e0_aux = tan(nu*0.5);
                    m_aux = e0_aux + (e0_aux*e0_aux*e0_aux)/3.;
                }
            }
        }
    }

    if(ecc < 1.0)
    {
        m_aux = remainder(m_aux, pi2);
        if(m_aux < 0.)
            m_aux += pi2;
        e0_aux = remainder(e0_aux, pi2);
    }

    *m = m_aux;
    *e0 = e0_aux;
}
