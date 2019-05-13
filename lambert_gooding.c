/**
 * @file lambert_gooding.c
 * @Autor Davide Pérez y Millán Santamaría
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "SAT_Const.h"
#include "MatlabUtils.h"
#include "unit.h"
#include "lambert_gooding.h"
void lambert_gooding(double *r1, double *r2, double tof, double mu, bool long_way, double multi_revs, double *v1, double *v2)
{
    double r1mag = norm(r1);
    double r2mag = norm(r2);

    if ( r1mag==0.0 || r2mag==0.0 || mu<=0.0 || tof<=0.0 )
    {
        printf("Error in solve_lambert_gooding: invalid input\n");
        return;
    }

    bool solution_exists[10];
    for(int i=0; i<10; i++)
        solution_exists[i] = false;
    double dr = r1mag - r2mag;
    double r1r2 = r1mag*r2mag;
    double *r1hat = vectorProductDouble(r1, (1.0/r1mag));
    double *r2hat = vectorProductDouble(r2,(1.0/r2mag));
    double *r1xr2 = cross(r1,r2);

    if(r1xr2[0]==0.0 && r1xr2[1]==0.0 && r1xr2[2]==0.0)
    {
        r1xr2[0] = 0.0;
        r1xr2[1] = 0.0;
        r1xr2[2] = 1.0;
    }
    double *r1xr2_hat = unit(r1xr2);

    double pa = acos(fmax(-1.0,fmin(1.0,dot(r1hat,r2hat))));

    int n_solutions = 0;
    double all_vt1[3][10], all_vt2[3][10];

    for(int i = 0; i<=multi_revs; i++)
    {
        int num_revs = i;
        double ta;
        double *rho;
        if(long_way)
        {
            ta = num_revs*2.0*pi + (2.0*pi-pa);
            rho   = vectorProductDouble(r1xr2_hat, -1.0);
        }
        else
        {
            ta = num_revs*2.0*pi + pa;
            rho = r1xr2_hat;
        }

        double *etai = cross(rho,r1hat);
        double *etaf = cross(rho,r2hat);

        double vri[2] = {0., 0.}, vti[2]= {0., 0.}, vrf[2]= {0., 0.}, vtf[2]= {0., 0.};
        int n = vlamb(mu,r1mag,r2mag,ta,tof,vri,vti,vrf,vtf);

        double vt1[3][2], vt2[3][2];
        double *aux;
        switch(n)
        {
        case 1 :
            aux = sumVector(vectorProductDouble(r1hat, vri[0]), vectorProductDouble(etai, vti[0]));
            vt1[0][0] = aux[0];
            vt1[1][0] = aux[1];
            vt1[2][0] = aux[2];
            aux = sumVector(vectorProductDouble(r2hat, vrf[0]), vectorProductDouble(etaf, vtf[0]));
            vt2[0][0] = aux[0];
            vt2[1][0] = aux[1];
            vt2[2][0] = aux[2];
            break;
        case 2 :
            aux = sumVector(vectorProductDouble(r1hat, vri[0]), vectorProductDouble(etai, vti[0]));
            vt1[0][0] = aux[0];
            vt1[1][0] = aux[1];
            vt1[2][0] = aux[2];
            aux = sumVector(vectorProductDouble(r2hat, vrf[0]), vectorProductDouble(etaf, vtf[0]));
            vt2[0][0] = aux[0];
            vt2[1][0] = aux[1];
            vt2[2][0] = aux[2];
            aux = sumVector(vectorProductDouble(r1hat, vri[1]), vectorProductDouble(etai, vti[1]));
            vt1[0][1] = aux[0];
            vt1[1][1] = aux[1];
            vt1[2][1] = aux[2];
            aux = sumVector(vectorProductDouble(r2hat, vrf[1]), vectorProductDouble(etaf, vtf[1]));
            vt2[0][1] = aux[0];
            vt2[1][1] = aux[1];
            vt2[2][1] = aux[2];
            break;
        }

        if(i == 0 && n==1)
        {
            all_vt1[0][0] = vt1[0][0];
            all_vt1[1][0] = vt1[1][0];
            all_vt1[2][0] = vt1[2][0];

            all_vt2[0][0] = vt2[0][0];
            all_vt2[1][0] = vt2[1][0];
            all_vt2[2][0] = vt2[2][0];

            solution_exists[0] = true;
            n_solutions ++;
        }
        else
        {
            switch(n)
            {
            case 1 :
                all_vt1[0][2*i-1] = vt1[0][0];
                all_vt1[1][2*i-1] = vt1[1][0];
                all_vt1[2][2*i-1] = vt1[2][0];

                all_vt2[0][2*i-1] = vt2[0][0];
                all_vt2[1][2*i-1] = vt2[1][0];
                all_vt2[2][2*i-1] = vt2[2][0];

                solution_exists[2*i-1] = true;///Igual esta mal el -1
                n_solutions ++;
                break;
            case 2:
                all_vt1[0][2*i-1] = vt1[0][0];
                all_vt1[1][2*i-1] = vt1[1][0];
                all_vt1[2][2*i-1] = vt1[2][0];

                all_vt2[0][2*i-1] = vt2[0][0];
                all_vt2[1][2*i-1] = vt2[1][0];
                all_vt2[2][2*i-1] = vt2[2][0];

                solution_exists[2*i-1] = true;

                all_vt1[0][2*i] = vt1[0][1];
                all_vt1[1][2*i] = vt1[1][1];
                all_vt1[2][2*i] = vt1[2][1];

                all_vt2[0][2*i] = vt2[0][1];
                all_vt2[1][2*i] = vt2[1][1];
                all_vt2[2][2*i] = vt2[2][1];

                solution_exists[2*i] = true;
                n_solutions ++;
            }
        }

    }

    int k = 0;
    for(int i = 0; i<n_solutions; i++)
    {
        if(solution_exists[i])
        {
            k = k+1;
            v1[0] = all_vt1[0][i];
            v1[1] = all_vt1[1][i];
            v1[2] = all_vt1[2][i];

            v2[0] = all_vt2[0][i];
            v2[1] = all_vt2[1][i];
            v2[2] = all_vt2[2][i];
        }
    }
}



double vlamb(double gm, double r1, double r2, double th, double tdelt, double *vri, double *vti, double *vrf, double *vtf)
{
    double thr2 = th;
    int m = 0;
    while(thr2 > 2*pi)
    {
        thr2 = thr2 - 2*pi;
        m = m + 1;
    }
    thr2 = thr2/2;

    double r1mag = r1;
    double r2mag = r2;
    double dr = r1mag-r2mag;
    dr = r1mag-r2mag;
    double r1r2 = r1mag*r2mag;
    double r1r2th = 4.0*r1r2*sin(thr2)*sin(thr2);
    double csq = dr*dr + r1r2th;
    double c = sqrt(csq);
    double s = (r1 + r2 + c)/2.0;
    double gms = sqrt(gm*s/2.0);
    double qsqfm1 = c/s;
    double q = sqrt(r1r2)*cos(thr2)/s;

    double rho, sig, n;
    if(c != 0.0)
    {
        rho = dr/c;
        sig = r1r2th/csq;
    }
    else
    {
        rho = 0.0;
        sig = 1.0;
    }

    double t = 4.0*gms*tdelt/(s*s);

    double *aux = NULL;
    aux = xlamb(m,q,qsqfm1,t);
    n = aux[0];
    double x1 = aux[1];
    double x2 = aux[2];
    double x;
    free(aux);
    for(int i = 0; i<n; i++)
    {
        if(i==0)
        {
            x = x1;
        }
        else
        {
            x = x2;
        }
        aux = tlamb(m,q,qsqfm1,x,-1);
        double qzminx = aux[1];
        double qzplx = aux[2];
        double zplqx = aux[3];
        double vt2 = gms*zplqx*sqrt(sig);
        double vr1 = gms*(qzminx - qzplx*rho)/r1;
        double vt1 = vt2/r1;
        double vr2 = -gms*(qzminx + qzplx*rho)/r2;
        vt2 = vt2/r2;
        vri[i] = vr1;
        vti[i] = vt1;
        vrf[i] = vr2;
        vtf[i] = vt2;
        free(aux);
    }
    /* esto que he hecho es muy malo */
    vri[1] = 0.;
    vti[1] = 0.;
    vrf[1] = 0.;
    vtf[1] = 0.;
    /* hay que buscar una solucion alternativa */
    return n;
}


double *tlamb(double m, double q, double qsqfm1, double x, double n)
{
    double sw = 0.4;
    double t = 0.0;

    bool lm1 = fabs(n - (-1)) < 1e-12;
    bool l1 = (n>=1);
    bool l2 = (n>=2);
    bool l3 = fabs(n - (3)) < 1e-12;
    double qsq = q*q;
    double xsq = x*x;
    double u = (1.0 - x)*(1.0 + x);

    double dt,d2t,d3t;
    if(!lm1)
    {
        dt=0.0;
        d2t=0.0;
        d3t=0.0;
    }
    double y=0., z=0., qx=0., aa=0., bb=0., a=0., b=0., g=0., f=0., fg1=0., term=0., fg1sq=0., towi1=0., told=0., qz=0.,qz2=0., u0i=0., u1i=0.,u2i=0.,u3i=0., twoi1=0., tq=0., tqsum=0., ttmold=0., p=0., v=0., tterm=0., tqterm=0.;
    if(lm1 || m>0.0 || x<0.0 || fabs(u)>sw)
    {
        y = sqrt(fabs(u));
        z = sqrt(qsqfm1 + qsq*xsq);
        qx = q*x;
        if(qx<=0.0)
        {
            a = z-qx;
            b = q*z -x;
        }
        if(qx<0.0 && lm1)
        {
            aa = qsqfm1/a;
            bb = qsqfm1*(qsq*u -xsq)/b;
        }
        if((qx==0.0 && lm1) ||qx>0.0)
        {
            aa = z + qx;
            bb = q*z +x;
        }
        if(qx>0.0)
        {
            a = qsqfm1/aa;
            b = qsqfm1*(qsq*u - xsq)/bb;
        }
        if(!lm1)
        {
            if(qx*u>=0.0)
            {
                g = x*z+q*u;
            }
            else
            {
                g = (xsq - qsq*u)/(x*z - q*u);
            }
            f = a*y;
            if(x<=1.0)
            {
                t = m*pi + atan2(f, g);
            }
            else
            {
                if(f>sw)
                {
                    t = log(f+g);
                }
                else
                {
                    fg1 = f/(g + 1.0);
                    term = 2.0*fg1;
                    fg1sq = fg1*fg1;
                    t = term;
                    twoi1 = 1.0;
                    //Duda linea 300
                    do
                    {
                        twoi1 = twoi1 + 2.0;
                        term = term*fg1sq;
                        told = t;
                        t = t + term/twoi1;
                    }
                    while(fabs(t-told) > 1e-12);
                }
            }
            t = 2.0*(t/y + b)/u;
            if(l1 && z!=0.0)
            {
                qz = q/z;
                qz2 = qz*qz;
                qz = qz*qz2;
                dt = (3.0*x*t - 4.0*(a + qx*qsqfm1)/z)/u;
                if(l2)
                {
                    d2t = (3.0*t + 5.0*x*dt + 4.0*qz*qsqfm1)/u;
                }
                if(l3)
                {
                    d3t = (8.0*dt + 7.0*x*d2t - 12.0*qz*qz2*x*qsqfm1)/u;
                }
            }
        }
        else
        {
            dt = b;
            d2t = bb;
            d3t = aa;
        }
    }
    else
    {
        u0i = 1.0;
        if(l1)
        {
            u1i = 1.0;
        }
        if(l2)
        {
            u2i = 1.0;
        }
        if(l3)
        {
            u3i = 1.0;
        }
        term = 4.0;
        tq = q*qsqfm1;
        int i = 0;
        if(q<0.5)
        {
            tqsum = 1.0 - q*qsq;
        }
        if(q>=0.5)
        {
            tqsum = (1.0/(1.0 + q) + q)*qsqfm1;
        }
        ttmold = term/3.0;
        t = ttmold*tqsum;
        do
        {
            i=i+1;
            p=i;
            u0i = u0i*u;
            if(l1 && i>1)
            {
                u1i = u1i*u;
            }
            if(l2 && i>2)
            {
                u2i = u2i*u;
            }
            if(l3 && i>3)
            {
                u3i = u3i*u;
            }
            term = term*(p - 0.5)/p;
            tq = tq*qsq;
            tqsum = tqsum + tq;
            told = t;
            tterm = term/(2.0*p + 3.0);
            tqterm = tterm*tqsum;
            t = t - u0i*((1.5*p + 0.25)*tqterm/(p*p - 0.25)-ttmold*tq);
            ttmold = tterm;
            tqterm = tqterm*p;
            if(l1)
            {
                dt = dt + tqterm*u1i;
            }
            if(l2)
            {
                d2t = d2t + tqterm*u2i*(p - 1.0);
            }
            if(l3)
            {
                d3t = d3t + tqterm*u3i*(p - 1.0)*(p - 2.0);
            }
        }
        while(i<n || (fabs(t -told) > 1e-12));
        if(l3)
        {
            d3t = 8.0*x*(1.5*d2t - xsq*d3t);
        }
        if(l2)
        {
            d2t = 2.0*(2.0*xsq*d2t - dt);
        }
        if(l1)
        {
            dt = -2.0*x*dt;
        }
        t = t/xsq;
    }
    double *dev = (double*) malloc(4*sizeof(double));
    dev[0] = t;
    dev[1] = dt;
    dev[2] = d2t;
    dev[3] = d3t;
    return dev;
}

double *xlamb(double m, double q, double qsqfm1, double tin)
{

    double tol = 3e-7;
    double c0  = 1.7;
    double c1  = 0.5;
    double c2  = 0.03;
    double c3  = 0.15;
    double c41 = 1.0;
    double c42 = 0.24;

    double thr2 = atan2(qsqfm1, 2.0*q)/pi;

    double xpl = 0.0;
    double x = 0.0;

    double n, t0, dt, d2t, d3t, tdiff, w, xm = 0., tmin=0., xmold, xtest=0., tdiffm=0., d2t2=0., tdiff0, ij, t;
    double *aux, *dev;
    dev = (double*) malloc(3*sizeof(double));
    int i=0;
    if(m == 0.0)
    {
        n = 1.0;
        aux = tlamb(m,q,qsqfm1,0.0,0.0);
        t0 = aux[0];
        dt = aux[1];
        d2t = aux[2];
        d3t = aux[3];

        tdiff = tin -t0;
        if(tdiff<=0.0)
        {
            x = t0*tdiff/(-4.0*tin);
        }
        else
        {
            x = -tdiff/(tdiff + 4.0);
            w = x + c0*sqrt(2.0*(1.0 - thr2));
            if(w<0.0)
            {
                x = x - sqrt(d8rt(-w))*(x + sqrt(tdiff/(tdiff + 1.5*t0)));
            }
            w = 4.0/(4.0 + tdiff);
            x = x*(1.0 + x*(c1*w - c2*x*sqrt(w)));
        }
        free(aux);
    }
    else
    {
        xm = 1.0/(1.5*(m + 0.5)*pi);
        if(thr2<0.5)
        {
            xm = d8rt(2.0*thr2)*xm;
        }
        if(thr2>0.5)
        {
            xm = (2.0 - d8rt(2.0 - 2.0*thr2))*xm;
        }

        do
        {
            aux = tlamb(m,q,qsqfm1,xm,3);
            tmin = aux[0];
            dt = aux[1];
            d2t = aux[2];
            d3t = aux[3];
            if(d2t != 0.0)
            {
                xmold = xm;
                xm = xm - dt*d2t/(d2t*d2t - dt*d3t/2.0);
                xtest = fabs(xmold/xm - 1.0);
            }
            i++;
            free(aux);
        }
        while(i<12 && d2t!=0.0 && xtest>tol);

        if(i>12)
        {
            dev[0] = -1.0;
            dev[1] = x;
            dev[2] = xpl;
            return dev;
        }
        tdiffm = tin - tmin;
        if(tdiffm<0.0)
        {
            dev[0] = 0.0;
            dev[1] = x;
            dev[2] = xpl;
            return dev;
        }
        else if(tdiffm==0.0)
        {
            x = xm;
            dev[0] = 1.0;
            dev[1] = x;
            dev[2] = xpl;
            return dev;
        }
        else
        {
            n = 3.0;
            if(d2t==0.0)
            {
                d2t = 6.0*m*pi;
            }
            x = sqrt(tdiffm/(d2t/2.0 + tdiffm/((1.0 - xm)*(1.0 - xm))));
            w = xm + x;
            w = w*4.0/(4.0 + tdiffm) + (1.0 - w)*(1.0 - w);
            x = x*(1.0 - (1.0 + m + c41*(thr2 - 0.5))/(1.0 + c3*m)*x*(c1*w + c2*x*sqrt(w))) + xm;
            d2t2 = d2t/2.0;
            if(x>=1.0)
            {
                n = 1.0;
                aux = tlamb(m,q,qsqfm1,0.0,0.0);
                t0 = aux[0];
                dt = aux[1];
                d2t = aux[2];
                d3t = aux[3];
                tdiff0 = t0 - tmin;
                tdiff = tin - t0;
                if(tdiff<=0.0)
                {
                    x = xm - sqrt(tdiffm/(d2t2 - tdiffm*(d2t2/tdiff0 - 1.0/(xm*xm))));
                }
                else
                {
                    x = -tdiff/(tdiff + 4.0);
                    ij = 200.0;
                    w = x + c0*sqrt(2.0*(1.0 - thr2));
                    if(w<0.0)
                    {
                        x = x - sqrt(d8rt(-w))*(x + sqrt(tdiff/(tdiff+1.5*t0)));
                    }
                    w = 4.0/(4.0 + tdiff);
                    x = x*(1.0 + (1.0 + m + c42*(thr2 - 0.5))/(1.0 + c3*m)*x*(c1*w - c2*x*sqrt(w)));
                    if(x<=-1.0)
                    {
                        n = n-1.0;
                        if(n==1.0)
                        {
                            x = xpl;
                        }
                    }
                }
                free(aux);
            }
        }
    }
    for(i=0; i<3; i++)
    {
        aux = tlamb(m,q,qsqfm1,x,2.0);
        t = aux[0];
        dt = aux[1];
        d2t = aux[2];
        d3t = aux[3];
        t = tin-t;
        if(dt!=0.0)
        {
            x = x + t*dt/(dt*dt + t*d2t/2.0);
        }
        free(aux);
    }
    if((fabs(n - 3)>1e-12))
    {
        dev[0] = n;
        dev[1] = x;
        dev[2] = xpl;
        return dev;
    }
    n = 2.0;
    xpl = x;
    aux = tlamb(m,q,qsqfm1,0.0,0);
    t0 = aux[0];
    dt = aux[1];
    d2t = aux[2];
    d3t = aux[3];
    tdiff0 = t0 - tmin;
    tdiff = tin - t0;
    if(tdiff<=0.0)
    {
        x = xm - sqrt(tdiffm/(d2t2 - tdiffm*(d2t2/tdiff0 - 1.0/(xm*xm))));
    }
    else
    {
        x = -tdiff/(tdiff + 4.0);
        ij = 200.0;
        w = x + c0*sqrt(2.0*(1.0 - thr2));
        if(w<0.0)
        {
            x = x - sqrt(d8rt(-w))*(x + sqrt(tdiff/(tdiff+1.5*t0)));
        }
        w = 4.0/(4.0 + tdiff);
        x = x*(1.0 + (1.0 + m + c42*(thr2 - 0.5))/(1.0 + c3*m)*x*(c1*w - c2*x*sqrt(w)));
        if(x<=-1.0)
        {
            n = n-1;
            if(fabs(n-1) < 1e-12)
            {
                x = xpl;
            }
        }
    }
    dev[0] = n;
    dev[1] = x;
    dev[2] = xpl;
    free(aux);
    return dev;
}

double d8rt(double x)
{
    double d8rt = sqrt(sqrt(sqrt(x)));
    return d8rt;
}
