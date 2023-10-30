#include <cmath>
#include <string>

#include "headers/configurations.h"
#include "headers/equations.h"

using namespace std;

// Introduce the notation:
//    dens_L, vel_L, pres_L = gas density, gas velocity, gas pressure on the LEFT,
//                          there are left border conditions.
//    dens_R, vel_R, pres_R = gas density, gas velocity, gas pressure on the RIGHT,
//                          there are right border conditions.
//    dens_L_star, dens_R_star = gas density on the left, right with *,
//    vel_star, pres_star = gas velocity, gas pressure with *,
//    s_vel_L, s_vel_R = sound velocity on the left, right,
//    adi_exp = adiabatic exponent (ratio of specific heats),
//    time is a time point under study.

void riemann_solve (double dens_L, double vel_L, double pres_L, double s_vel_L, \
                      double dens_R, double vel_R, double pres_R, double s_vel_R, \
                      double pres_star, double adi_exp, double time, \
                      double* velocity, double* density, double* pressure, double* x, \
                      int N, double x_L, double x_R, string configuration)
{
    double dens_L_star, dens_R_star, gas_velocity_star;
    double gas_velocity_RW_H_on_the_left, gas_velocity_RW_H_on_the_right, gas_velocity_RW_T_on_the_left, gas_velocity_RW_T_on_the_right, D_on_the_left, D_on_the_right;
    double dist_RW_HL, dist_RW_HR, dist_RW_TL, dist_RW_TR, dist_D_on_the_left, dist_D_on_the_right, dist_star;

    if (configuration == "A")
    {
        dens_L_star = A_rho_L_star (dens_L, pres_L, pres_star, adi_exp);
        dens_R_star = A_rho_R_star (dens_R, pres_R, pres_star, adi_exp);
        gas_velocity_star = A_v_xL_star (vel_L, s_vel_L, pres_L, pres_star, adi_exp);
        gas_velocity_RW_H_on_the_left = vel_L - s_vel_L; // leading edge speed
        gas_velocity_RW_T_on_the_left = gas_velocity_star - s_vel_L * pow(pres_star / pres_L, \
                                        ((adi_exp - 1) / (2 * adi_exp))); // falling edge speed
        D_on_the_right = A_D (vel_R, pres_R, s_vel_R, pres_star, adi_exp); // SW propagation speed
        // distances
        dist_RW_HL = gas_velocity_RW_H_on_the_left * time;
        dist_RW_TL = gas_velocity_RW_T_on_the_left * time;
        dist_D_on_the_right = D_on_the_right * time;
        dist_star = gas_velocity_star * time;
    }
    else if (configuration == "A*")
    {
        dens_L_star = A1_rho_L_star (dens_L, pres_L, pres_star, adi_exp);
        dens_R_star = A1_rho_R_star (dens_R, pres_R, pres_star, adi_exp);
        gas_velocity_star = vel_L - (1 / (dens_L * s_vel_L)) * (pres_star - pres_L) / \
                            sqrt((((adi_exp + 1) * pres_star) / (2 * adi_exp * pres_L)) + ((adi_exp - 1) / \
                            (2 * adi_exp)));
        gas_velocity_RW_H_on_the_right = vel_R + s_vel_R;
        gas_velocity_RW_T_on_the_right = gas_velocity_star + s_vel_R * pow(pres_star / pres_R,\
                                         ((adi_exp - 1) / (2 * adi_exp)));
        D_on_the_left = A1_D (vel_L, pres_L, s_vel_L, pres_star, adi_exp);
        // distances
        dist_RW_HR = gas_velocity_RW_H_on_the_right * time;
        dist_RW_TR = gas_velocity_RW_T_on_the_right * time;
        dist_D_on_the_left = D_on_the_left * time;
        dist_star = gas_velocity_star * time;
    }
    else if (configuration == "B")
    {
        dens_L_star = A1_rho_L_star (dens_L, pres_L, pres_star, adi_exp);
        dens_R_star = A_rho_R_star (dens_R, pres_R, pres_star, adi_exp);
        gas_velocity_star = B_v_xR_star (dens_R, vel_R, s_vel_R, pres_R, pres_star, adi_exp);
        // SW front velocities
        D_on_the_left = A1_D (vel_L, pres_L, s_vel_L, pres_star, adi_exp);
        D_on_the_right = A_D (vel_R, pres_R, s_vel_R, pres_star, adi_exp);
        // distances
        dist_D_on_the_left = D_on_the_left * time;
        dist_D_on_the_right = D_on_the_right * time;
        dist_star = gas_velocity_star * time;
    }
    else if (configuration == "C")
    {
        dens_L_star = A_rho_L_star (dens_L, pres_L, pres_star, adi_exp);
        dens_R_star = A1_rho_R_star (dens_R, pres_R, pres_star, adi_exp);
        gas_velocity_star = A_v_xL_star (vel_L, s_vel_L, pres_L, pres_star, adi_exp);
        gas_velocity_RW_H_on_the_left = vel_L - s_vel_L;
        gas_velocity_RW_T_on_the_left = gas_velocity_star - s_vel_L * pow(pres_star / pres_L,\
                                        ((adi_exp - 1) / (2 * adi_exp)));
        gas_velocity_RW_H_on_the_right = vel_R + s_vel_R;
        gas_velocity_RW_T_on_the_right = gas_velocity_star + s_vel_R * pow(pres_star / pres_R,\
                                         ((adi_exp - 1) / (2 * adi_exp)));
        // distances
        dist_RW_HL = gas_velocity_RW_H_on_the_left * time;
        dist_RW_TL = gas_velocity_RW_T_on_the_left * time;
        dist_RW_HR = gas_velocity_RW_H_on_the_right * time;
        dist_RW_TR = gas_velocity_RW_T_on_the_right * time;
        dist_star = gas_velocity_star * time;
    }

    double step = (x_R - x_L) / N;

    for (size_t i = 0; i != N; ++i)
    {
        x[i] = x_L + i * step;
        configuration_implementation (dens_L, vel_L, pres_L, s_vel_L, \
                                      dens_R, vel_R, pres_R, s_vel_R, \
                                      pres_star, dens_L_star, dens_R_star, gas_velocity_star, dist_star, \
                                      adi_exp, time, configuration, \
                                      dist_RW_HL, dist_RW_TL, dist_RW_HR, dist_RW_TR, \
                                      dist_D_on_the_left, dist_D_on_the_right, \
                                      x[i], velocity[i], density[i], pressure[i]);
    }
}
