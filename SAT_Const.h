/**
 * @file SAT_Const.h
 * @Autor Davide Pérez y Millán Santamaría
 * @brief Es el fichero de constantes
 */
#ifndef SAT_CONST_H
#define SAT_CONST_H

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define pi2 2.0*M_PI
#define Rad M_PI/180.0
#define Deg 180.0/M_PI
#define Arcs  3600.0*180.0/M_PI

#define MJD_J2000 51544.5
#define T_B1950 -0.500002018
#define c_light 299792457.999999984
#define AU 149597870659.999996

#define R_Earth 6378.137e3
#define f_Earth 1.0/298.257223563
#define R_Sun 696000.0e3
#define R_Moon 1738.0e3

#define omega_Earth 7.2921158553e-5

#define GM_Earth 398600.4418e9
#define GM_Sun 1.327124399354841e20
#define GM_Moon GM_Earth/81.3005869999999931
#define GM_Mercury 22032.08047272131e9
#define GM_Venus 324858.7656168717e9
#define GM_Mars 42828.28658876890e9
#define GM_Jupiter 126712597.0817946e9
#define GM_Saturn 37939519.70882996e9
#define GM_Uranus 5780158.533597719e9
#define GM_Neptune 6871307.771479524e9
#define GM_Pluto 1020.864920706286e9

#define P_Sol 4.560e-6
#endif // SAT_CONST_H
