#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

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

// Choice between configurations A, B, C and A* depending on the initial conditions

string configuration_selection (const double adi_exp, double dens_L, double vel_L, double pres_L, double s_vel_L,\
                                                      double dens_R, double vel_R, double pres_R, double s_vel_R)
{
    string configuration;

    double v_A = 2 * s_vel_L / (adi_exp - 1) * (pow((pres_R / pres_L), ((adi_exp - 1) / (2 * adi_exp)))\
                                                - 1);       // < 0
    double v_B = (pres_L - pres_R) / (dens_R * s_vel_R) / sqrt((adi_exp + 1) * pres_L / (2 * adi_exp * \
                  pres_R) + (adi_exp - 1) / (2 * adi_exp)); // > 0
    double v_C = -2 * (s_vel_L + s_vel_R) / (adi_exp - 1);  // < 0

    if ((vel_L - vel_R) < v_C)               // Configuration C particular
    {
        cerr << "A vacuum region (zero gas density) has formed between the rarefaction waves; \
                 the solution cannot be continued.\n";
        
        return "Error";
    }
    else
    {   
        if ((v_A < (vel_L - vel_R)) && ((vel_L - vel_R) < v_B))
        {
            if (pres_L > pres_R)
            {
                configuration = "A";         // Configuration A
            }
            else
            {
                configuration = "A*";        // Configuration A* (inverted)
            }
        }
        else if ((vel_L - vel_R) > v_B)      // Configuration B
        {
            configuration = "B";
        }
        else if ((vel_L - vel_R) < v_A)      // Configuration C
        {
            configuration = "C";
        }

        return configuration;
    }
}

void configuration_implementation (double dens_L, double vel_L, double pres_L, double s_vel_L, \
                                   double dens_R, double vel_R, double pres_R, double s_vel_R, \
                                   double pres_star, double dens_L_star, double dens_R_star, double vel_star, double dist_star, \
                                   double adi_exp, double time, string configuration, \
                                   double dist_RW_HL, double dist_RW_TL, double dist_RW_HR, double dist_RW_TR, \
                                   double dist_D_on_the_left, double dist_D_on_the_right, \
                                   double x, double& VEL, double& DENS, double& PRES)
{
    if (configuration == "A")
    {
        if (dist_star <= x)
        {
            if (x < dist_D_on_the_right)
            {
                VEL = vel_star;
                DENS = dens_R_star;
                PRES = pres_star;
            }
            else
            {
                VEL = vel_R;
                DENS = dens_R;
                PRES = pres_R;
            }
        }
        else
        {
            if (dist_RW_TL < x)
            {
                VEL = vel_star;
                DENS = dens_L_star;
                PRES = pres_star;
            }
            else if (dist_RW_HL <= x)
            {
                VEL = A_v_RW (x, vel_L, s_vel_L, adi_exp, time);
                DENS = A_rho_RW (x, dens_L, vel_L, s_vel_L, adi_exp, time);
                PRES = A_p_RW (x, vel_L, pres_L, s_vel_L, adi_exp, time);
            }
            else
            {
                VEL = vel_L;
                DENS = dens_L;
                PRES = pres_L;
            }   
        }
    }
    else if (configuration == "A*")
    {
        if (dist_star < x)
        {
            if (x < dist_RW_TR)
            {
                VEL = vel_star;
                DENS = dens_R_star;
                PRES = pres_star;
            }
            else if (x <= dist_RW_HR)
            {
                VEL = A1_v_RW (x, vel_R, s_vel_R, adi_exp, time);
                DENS = A1_rho_RW (x, dens_R, vel_R, s_vel_R, adi_exp, time);
                PRES = A1_p_RW (x, vel_R, pres_R, s_vel_R, adi_exp, time);
            }
            else
            {
                VEL = vel_R;
                DENS = dens_R;
                PRES = pres_R;
            }
        }
        else
        {
            if (x <= dist_D_on_the_left)
            {
                VEL = vel_L;
                DENS = dens_L;
                PRES = pres_L;
            }
            else
            {
                VEL = vel_star;
                DENS = dens_L_star;
                PRES = pres_star;
            }
        }
    }
    else if (configuration == "B")
    {
        if (dist_star < x)
        {
            if (x <= dist_D_on_the_right)
            {
                VEL = vel_star;
                DENS = dens_R_star;
                PRES = pres_star;
            }
            else
            {
                VEL = vel_R;
                DENS = dens_R;
                PRES = pres_R;
            }
        }
        else
        {
            if (x <= dist_D_on_the_left)
            {
                VEL = vel_L;
                DENS = dens_L;
                PRES = pres_L;
            }
            else
            {
                VEL = vel_star;
                DENS = dens_L_star;
                PRES = pres_star;
            }
        }
    }
    else if (configuration == "C")
    {
        if (dist_star < x)
        {
            if (x < dist_RW_TR)
            {
                VEL = vel_star;
                DENS = dens_R_star;
                PRES = pres_star;
            }
            else if (x <= dist_RW_HR)
            {
                VEL = C_v_RW (x, vel_R, s_vel_R, adi_exp, time);
                DENS = C_rho_RW (x, dens_R, vel_R, s_vel_R, adi_exp, time);
                PRES = C_p_RW (x, vel_R, pres_R, s_vel_R, adi_exp, time);
            }
            else
            {
                VEL = vel_R;
                DENS = dens_R;
                PRES = pres_R;
            }
        }
        else
        {
            if (dist_RW_TL < x)
            {
                VEL = vel_star;
                DENS = dens_L_star;
                PRES = pres_star;
            }
            else if (dist_RW_HL <= x)
            {
                VEL = A_v_RW (x, vel_L, s_vel_L, adi_exp, time);
                DENS = A_rho_RW (x, dens_L, vel_L, s_vel_L, adi_exp, time);
                PRES = A_p_RW (x, vel_L, pres_L, s_vel_L, adi_exp, time);
            }
            else
            {
                VEL = vel_L;
                DENS = dens_L;
                PRES = pres_L;
            }
        }
    }
    else
    {
        VEL = vel_star;
        DENS = dens_R_star;
        PRES = pres_star;
    }
}
