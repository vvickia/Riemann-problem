#include <string>

#ifndef NEWTON_INCLUDED
#define NEWTON_INCLUDED

double fx (double x, double dens_L, double vel_L, double pres_L, double s_vel_L, \
                     double dens_R, double vel_R, double pres_R, double s_vel_R, \
           const double adi_exp, std::string configuration);

double dfx (double x, double dens_L, double pres_L, double s_vel_L, \
                      double dens_R, double pres_R, double s_vel_R, \
            const double adi_exp, std::string configuration);

double Newtons_iterative_method (double dens_L, double vel_L, double pres_L, double s_vel_L, \
                                 double dens_R, double vel_R, double pres_R, double s_vel_R, \
                                 const double adi_exp, std::string configuration);

#endif
