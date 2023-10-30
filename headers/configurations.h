#include <string>

#ifndef CONFIG_INCLUDED
#define CONFIG_INCLUDED

std::string configuration_selection (const double adi_exp, double dens_L, double vel_L, double pres_L, double s_vel_L,\
                                                           double dens_R, double vel_R, double pres_R, double s_vel_R);
void configuration_implementation (double dens_L, double vel_L, double pres_L, double s_vel_L, \
                                   double dens_R, double vel_R, double pres_R, double s_vel_R, \
                                   double pres_star, double dens_L_star, double dens_R_star, double vel_star, double dist_star, \
                                   double adi_exp, double time, std::string configuration, \
                                   double dist_RW_HL, double dist_RW_TL, double dist_RW_HR, double dist_RW_TR, \
                                   double dist_D_on_the_left, double dist_D_on_the_right, \
                                   double x, double& VEL, double& DENS, double& PRES);

#endif