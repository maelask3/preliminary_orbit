/**
 * @file timediff.c
 * @authors Davide Pérez y Millán Santamaría
 */
#include "timediff.h"
/**
 * @brief Diferencias entre tiempo (s)
 * @param UT1_UTC Diferencia UT1-UTC (s)
 * @param TAI_UTC Diferencia TAI-UTC (s)
 * @param UT1_TAI (Salida) Diferencia UT1-TAI (s)
 * @param UTC_GPS (Salida) Diferencia UTC-GPS (s)
 * @param UT1_GPS (Salida) Diferencia UT1-GPS (s)
 * @param TT_UTC (Salida) Diferencia TT-UTC (s)
 * @param GPS_UTC (Salida) Diferencia GPS-UTC (s)
 */
void timediff(double UT1_UTC, double TAI_UTC, double *UT1_TAI, double *UTC_GPS, double *UT1_GPS, double *TT_UTC, double *GPS_UTC)
{
    double TT_TAI = 32.184;
    double GPS_TAI = -19.0;
    /*double TT_GPS = TT_TAI - GPS_TAI;
    double TAI_GPS = -GPS_TAI;*/
    *UT1_TAI = UT1_UTC-TAI_UTC;
    double UTC_TAI = -TAI_UTC;
    *UTC_GPS = UTC_TAI-GPS_TAI;
    *UT1_GPS = *UT1_TAI-GPS_TAI;
    *TT_UTC = TT_TAI - UTC_TAI;
    *GPS_UTC = GPS_TAI-UTC_TAI;
}
