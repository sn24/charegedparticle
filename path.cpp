#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include<iostream>
using namespace std;
int func (double t, const double y[], double f[], void *params)           //initialize the system of the differential equations
{
    double mu0= *(double *)params;
    double mu1 = *(double *)params;
    double mu2= *(double *)params;
    f[0] = mu2*y[1]-y[2]*mu1;
    f[1] = mu0*y[2]-y[0]*mu2;
    f[2] = mu1*y[0]-y[1]*mu0;
    return GSL_SUCCESS;
}

int main ()
{
    double q,m,bx,by,bz,vx0,vy0,vz0;                                      // take the user inputs
    cin>>q>>m>>bx>>by>>bz>>vx0>>vy0>>vz0;
    double mu0 = q/m*bx;
    double mu1 = q/m*by;
    double mu2 = q/m*bz;
    gsl_odeiv2_system sys = {func, NULL, 3, &mu};
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    int i;
    double t = 0.0, t1 = 100.0;
    double y[3] = { vx0, vy0,vz0 };

    for (i = 1; i <= 100; i++)                                           // apply runge-kutta method
    {
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }

    }

    printf ( y[0], y[1],y[2]);
    gsl_odeiv2_driver_free (d);
    return 0;
}
