#include "IERS.h"
#include "SAT_Const.h"
#include <math.h>
#include <stdlib.h>

void IERS(double **eop, int eop_length, double Mjd_UTC, char interp, double *UT1_UTC, double *TAI_UTC, double *x_pole, double *y_pole, double *ddpsi, double *ddeps)
{
    if(interp == 'l')
    {
        double *preeop = NULL;
        double *nexteop = NULL;
        double mj = floor(Mjd_UTC);

        for(int i=0; i<eop_length; i++)
        {
            if(fabs(eop[i][3] - mj) < DELTA)
            {
                preeop = eop[i];
                nexteop = eop[i+1];
                break;
            }
        }

        double mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
        double fixf = mfme/1440;

        if(!preeop || !nexteop)
            return;

        *UT1_UTC = preeop[6]+(nexteop[6]-preeop[6])*fixf;
        *TAI_UTC = preeop[12];
        *x_pole = preeop[4]+(nexteop[4]-preeop[4])*fixf;
        *y_pole = preeop[5]+(nexteop[5]-preeop[5])*fixf;
        *ddpsi = preeop[8]+(nexteop[8]-preeop[8])*fixf;
        *ddeps = preeop[9]+(nexteop[9]-preeop[9])*fixf;

        *x_pole /= Arcs;
        *y_pole /= Arcs;
        *ddpsi /= Arcs;
        *ddeps /= Arcs;

    } else if(interp == 'n') {

        double mj = floor(Mjd_UTC);
        double *cureop = NULL;
        for(int i=0; i<eop_length; i++)
        {
            if(fabs(eop[i][3] - mj) < DELTA)
            {
                cureop = eop[i];
                break;
            }
        }

        if(!cureop)
            return;

        *UT1_UTC = cureop[6];
        *TAI_UTC = cureop[12];
        *x_pole = cureop[4]/Arcs;
        *y_pole = cureop[5]/Arcs;
        *ddpsi = cureop[8]/Arcs;
        *ddeps = cureop[9]/Arcs;
    }
}
