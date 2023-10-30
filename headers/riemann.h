#include <string>

#ifndef RIEMANN_INCLUDED
#define RIEMANN_INCLUDED

void riemann_solve (double dens_L, double vel_L, double pres_L, double s_vel_L, \
            double dens_R, double vel_R, double pres_R, double s_vel_R, \
            double pres_star, double adi_exp, double time, \
            double *velocity, double *density, double *pressure, double *x, \
            int N, double x_L, double x_R, std::string configuration);

#endif