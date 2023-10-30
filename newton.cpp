#include <cmath>
#include <string>

#include "headers/equations.h"

using namespace std;

// Introduce the notation:
//    dens_L, vel_L, pres_L = gas density, gas velocity, gas pressure on the LEFT,
//                          there are left border conditions.
//    dens_R, vel_R, pres_R = gas density, gas velocity, gas pressure on the RIGHT,
//                          there are right border conditions.
//    s_vel_L, s_vel_R = sound velocity on the left, right,
//    adi_exp = adiabatic exponent (ratio of specific heats).

// Let's implement Newton's iterative method for this problem

#define eps 0.000001

// Computed function:
double fx (double x, double dens_L, double vel_L, double pres_L, double s_vel_L, \
                     double dens_R, double vel_R, double pres_R, double s_vel_R, \
           const double adi_exp, string configuration)
{
    double equation = 0;

    if (configuration == "A")
    {
        equation = A_v_xL_star (vel_L, s_vel_L, pres_L, x, adi_exp) - \
                   A_v_xR_star (dens_R, vel_R, s_vel_R, pres_R, x, adi_exp);
    }
    else if (configuration == "A*")
    {
        equation = A1_v_xL_star (dens_L, vel_L, s_vel_L, pres_L, x, adi_exp) - \
                   A1_v_xR_star (vel_R, s_vel_R, pres_R, x, adi_exp);
    }
    else if (configuration == "B")
    {
        equation = B_v_xL_star (dens_L, vel_L, s_vel_L, pres_L, x, adi_exp) - \
                   B_v_xR_star (dens_R, vel_R, s_vel_R, pres_R, x, adi_exp);
    }
    else if (configuration == "C")
    {
        equation = C_v_xL_star (vel_L, s_vel_L, pres_L, x, adi_exp) - \
                   C_v_xR_star (vel_R, s_vel_R, pres_R, x, adi_exp);
    }

    return equation;
}

// Derivative of a function:
double dfx (double x, double dens_L, double pres_L, double s_vel_L, double dens_R,\
            double pres_R, double s_vel_R, const double adi_exp, string configuration)
{
    double equation = 0;

    if (configuration == "A")
    {
        equation = A_dv (pres_L, x, s_vel_L, dens_R, pres_R, s_vel_R, adi_exp);
    }
    else if (configuration == "A*")
    {
        equation = A1_dv (dens_L, s_vel_L, pres_L, x, pres_R, s_vel_R, adi_exp);
    }
    else if (configuration == "B")
    {
        equation = B_dv (dens_L, s_vel_L, pres_L, x, dens_R, s_vel_R, pres_R, adi_exp);
    }
    else if (configuration == "C")
    {
        equation = C_dv (s_vel_L, pres_L, x, s_vel_R, pres_R, adi_exp);
    }

    return equation;
}

double Newtons_iterative_method (double dens_L, double vel_L, double pres_L,\
                                 double s_vel_L, double dens_R, double vel_R, double pres_R,\
                                 double s_vel_R, const double adi_exp, string configuration)
{
    // double p_0 = (pres_L + pres_R) / 2;
    double p_0 = (dens_R * s_vel_R * pres_L + dens_L * s_vel_L * pres_R + (vel_L - vel_R) * \
                 dens_L * dens_R * s_vel_L * s_vel_R) / (dens_L * s_vel_L + dens_R * s_vel_R);

    double p_1 = p_0 - fx (p_0, dens_L, vel_L, pres_L, s_vel_L, dens_R, vel_R, pres_R, s_vel_R, \
                            adi_exp, configuration) / \
                       dfx(p_0, dens_L, pres_L, s_vel_L, dens_R, pres_R, s_vel_R, \
                            adi_exp, configuration); // 1st approach

    int max_iterations = 10000;
    int counter = 0;

    while (abs(p_1 - p_0) > eps) // eps accuracy has not yet been achieved
    { 
        p_0 = p_1;
        p_1 = p_0 - fx (p_0, dens_L, vel_L, pres_L, s_vel_L, dens_R, vel_R, pres_R, s_vel_R, \
                            adi_exp, configuration) / \
                    dfx(p_0, dens_L, pres_L, s_vel_L, dens_R, pres_R, s_vel_R, \
                            adi_exp, configuration); // subsequent approaches
        counter++;
        if (counter == max_iterations)
            break;
    }

    return p_1;
}
