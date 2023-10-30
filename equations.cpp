#include <cmath>

using namespace std;

// The area between the RW and the SW containing a CD (configuration "A"):

double A_v_xL_star (double vel_L, double s_vel_L, double pres_L, double pres_unk, double adi_exp)
{
    double eq = vel_L - 2 * s_vel_L / (adi_exp - 1) * (pow((pres_unk / pres_L), ((adi_exp - 1) / \
                2 / adi_exp)) - 1);

    return eq;
}

double A_rho_L_star (double dens_L, double pres_L, double pres_unk, double adi_exp)
{
    double eq = dens_L * pow(pres_unk / pres_L, 1 / adi_exp);

    return eq;
}

double A_v_xR_star (double dens_R, double vel_R, double s_vel_R, double pres_R, double pres_unk, double adi_exp)
{
    double eq = vel_R + 1 / (dens_R * s_vel_R) * (pres_unk - pres_R) / sqrt((adi_exp + 1) / 2 / \
                adi_exp * pres_unk / pres_R + (adi_exp - 1) / 2 / adi_exp);

    return eq;
}

double A_rho_R_star (double dens_R, double pres_R, double pres_unk, double adi_exp)
{
    double eq = dens_R * ((adi_exp + 1) * (pres_unk / pres_R) + adi_exp - 1) / (adi_exp + 1 + \
                (adi_exp - 1) * (pres_unk / pres_R));

    return eq;
}

double A_D (double vel_R, double pres_R, double s_vel_R, double pres_unk, double adi_exp)
{
    double eq = vel_R + s_vel_R * sqrt(((adi_exp + 1) * pres_unk / (2 * adi_exp * pres_R)) + \
                ((adi_exp - 1) / (2 * adi_exp)));

    return eq;
}

// Inverted onfiguration "A" ("A*"):

double A1_v_xL_star (double dens_L, double vel_L, double s_vel_L, double pres_L, double pres_unk, double adi_exp)
{
    double eq = vel_L - 1 / (dens_L * s_vel_L) * (pres_unk - pres_L) / sqrt(((adi_exp + 1) / 2 / \
                adi_exp * pres_unk / pres_L) + (adi_exp - 1) / 2 / adi_exp);

    return eq;
}

double A1_rho_L_star (double dens_L, double pres_L, double pres_unk, double adi_exp)
{
    double eq = dens_L * ((adi_exp + 1) * (pres_unk / pres_L) + adi_exp - 1) / (adi_exp + 1 + \
                (adi_exp - 1) * (pres_unk / pres_L));

    return eq;
}

double A1_v_xR_star (double vel_R, double s_vel_R, double pres_R, double pres_unk, double adi_exp)
{
    double eq = vel_R + 2 * s_vel_R / (adi_exp - 1) * (pow((pres_unk / pres_R), ((adi_exp - 1) / \
                2 / adi_exp)) - 1);

    return eq;
}

double A1_rho_R_star (double dens_R, double pres_R, double pres_unk, double adi_exp)
{
    double eq = dens_R * pow(pres_unk / pres_R, 1 / adi_exp);

    return eq;
}

double A1_D (double vel_L, double pres_L, double s_vel_L, double pres_unk, double adi_exp)
{
    double eq = vel_L - s_vel_L * sqrt(((adi_exp + 1) * pres_unk / (2 * adi_exp * pres_L)) + \
                ((adi_exp - 1) / (2 * adi_exp)));

    return eq;
}

// The area between the fronts of SWs (configuration "B"):

double B_v_xL_star (double dens_L, double vel_L, double s_vel_L, double pres_L, double pres_unk, double adi_exp)
{
    double eq = vel_L - 1 / (dens_L * s_vel_L) * (pres_unk - pres_L) / sqrt((adi_exp + 1) / 2 / \
                adi_exp * pres_unk / pres_L + (adi_exp - 1) / 2 / adi_exp);

    return eq;
}

double B_v_xR_star (double dens_R, double vel_R, double s_vel_R, double pres_R, double pres_unk, double adi_exp)
{
    double eq = A_v_xR_star (dens_R, vel_R, s_vel_R, pres_R, pres_unk, adi_exp);

    return eq;
}

// Configuration "C"
// The area between the RW and the RW to the left of the CD:

double C_v_xL_star (double vel_L, double s_vel_L, double pres_L, double pres_unk, double adi_exp)
{
    double eq = A_v_xL_star (vel_L, s_vel_L, pres_L, pres_unk, adi_exp);

    return eq;
}

// The area between the CD and the right RW:

double C_v_xR_star (double vel_R, double s_vel_R, double pres_R, double pres_unk, double adi_exp)
{
    double eq = vel_R + 2 * s_vel_R / (adi_exp - 1) * (pow((pres_unk / pres_R), ((adi_exp - 1) / \
                2 / adi_exp)) - 1);

    return eq;
}


// Let's take derivatives!
//     Suppose we have already subtracted the right from the left side.
//     We will write derivatives for ready-made expressions:

double A_dv (double pres_L, double pres_unk, double s_vel_L, \
             double dens_R, double pres_R, double s_vel_R, double adi_exp)
{
    double eq = -s_vel_L / (adi_exp * pres_unk) * pow((pres_unk / pres_L), ((adi_exp - 1) / (2 * adi_exp))) \
                + (adi_exp + 1) * (pres_unk - pres_R) / (4 * s_vel_R * adi_exp * pres_R * dens_R * \
                pow(((adi_exp + 1) * pres_unk / (2 * adi_exp * pres_R) + (adi_exp - 1) / (2 * adi_exp)), \
                3 / 2)) - 1 / (s_vel_R * dens_R * sqrt((adi_exp + 1) * pres_unk / (2 * adi_exp * pres_R) + \
                (adi_exp - 1) / (2 * adi_exp)));

    return eq;
}

double A1_dv (double dens_L, double s_vel_L, double pres_L, double pres_unk, \
              double pres_R, double s_vel_R, double adi_exp)
{
    double eq = (adi_exp + 1) * (pres_unk - pres_L) / (4 * s_vel_L * adi_exp * pres_L * dens_L * \
                pow(((adi_exp + 1) * pres_unk / (2 * adi_exp * pres_L) + (adi_exp - 1) / (2 * adi_exp)), \
                3 / 2)) - 1 / (s_vel_L * dens_L * sqrt((adi_exp + 1) * pres_unk / (2 * adi_exp * pres_L) + \
                (adi_exp - 1) / (2 * adi_exp))) - s_vel_R / (adi_exp * pres_unk) * \
                pow((pres_unk / pres_R), ((adi_exp - 1) / (2 * adi_exp)));

    return eq;
}

double B_dv (double dens_L, double s_vel_L, double pres_L, double pres_unk, \
             double dens_R, double s_vel_R, double pres_R, double adi_exp)
{
    double eq = (adi_exp + 1) / (4 * adi_exp) * ((pres_unk - pres_L) / (s_vel_L * pres_L * dens_L * \
                pow(((adi_exp + 1) * pres_unk / (2 * adi_exp * pres_L) + (adi_exp - 1) / (2 * adi_exp)), \
                3 / 2)) + (pres_unk - pres_R) / (s_vel_R * pres_R * dens_R * pow(((adi_exp + 1) * \
                pres_unk / (2 * adi_exp * pres_R) + (adi_exp - 1) / (2 * adi_exp)), 3 / 2))) - 1 / \
                (s_vel_L * dens_L * sqrt((adi_exp + 1) * pres_unk / (2 * adi_exp * pres_L) + \
                (adi_exp - 1) / (2 * adi_exp))) - 1 / (s_vel_R * dens_R * sqrt((adi_exp + 1) * pres_unk / \
                (2 * adi_exp * pres_R) + (adi_exp - 1) / (2 * adi_exp)));

    return eq;
}

double C_dv (double s_vel_L, double pres_L, double pres_unk, \
             double s_vel_R, double pres_R, double adi_exp)
{
    double eq = -1 / (adi_exp * pres_unk) * (s_vel_L * pow((pres_unk / pres_L), ((adi_exp - 1) / \
                (2 * adi_exp))) + s_vel_R * pow((pres_unk / pres_R), ((adi_exp - 1) / (2 * adi_exp))));

    return eq;
}


// Other formulas

// RW propagating to the left (configuration "A"):

double A_v_RW (double x, double vel_L, double s_vel_L, double adi_exp, double time)
{
    double eq = vel_L * (adi_exp - 1) / (adi_exp + 1) + 2 * (x / time + s_vel_L) / (adi_exp + 1);

    return eq;
}

double A_rho_RW (double x, double dens_L, double vel_L, double s_vel_L, double adi_exp, double time)
{
    double eq = dens_L * pow(((2 / (adi_exp + 1)) + (adi_exp - 1) * (vel_L - (x / time)) / \
                ((adi_exp + 1) * s_vel_L)), (2 / (adi_exp - 1)));

    return eq;
}

double A_p_RW (double x, double vel_L, double pres_L, double s_vel_L, double adi_exp, double time)
{
    double eq = pres_L * pow(((2 / (adi_exp + 1)) + (adi_exp - 1) * (vel_L - (x / time)) / \
                ((adi_exp + 1) * s_vel_L)), ((2 * adi_exp) / (adi_exp - 1)));

    return eq;
}

// Configuration "A*":

double A1_v_RW (double x, double vel_R, double s_vel_R, double adi_exp, double time)
{
    double eq = vel_R * (adi_exp - 1) / (adi_exp + 1) + 2 * (x / time - s_vel_R) / (adi_exp + 1);

    return eq;
}

double A1_rho_RW (double x, double dens_R, double vel_R, double s_vel_R, double adi_exp, double time)
{
    double eq = dens_R * pow(((2 / (adi_exp + 1)) - (adi_exp - 1) * (vel_R - (x / time)) / \
                ((adi_exp + 1) * s_vel_R)), (2 / (adi_exp - 1)));

    return eq;
}

double A1_p_RW (double x, double vel_R, double pres_R, double s_vel_R, double adi_exp, double time)
{
    double eq = pres_R * pow(((2 / (adi_exp + 1)) - (adi_exp - 1) * (vel_R - (x / time)) / \
                ((adi_exp + 1) * s_vel_R)), ((2 * adi_exp) / (adi_exp - 1)));

    return eq;
}

// RW propagating to the right (configuration "C"):

double C_v_RW (double x, double vel_R, double s_vel_R, double adi_exp, double time)
{
    double eq = vel_R * (adi_exp - 1) / (adi_exp + 1) + 2 * (x / time - s_vel_R) / (adi_exp + 1);

    return eq;
}

double C_rho_RW (double x, double dens_R, double vel_R, double s_vel_R, double adi_exp, double time)
{
    double eq = dens_R * pow(((2 / (adi_exp + 1)) - (adi_exp - 1) * (vel_R - (x / time)) / \
                ((adi_exp + 1) * s_vel_R)), (2 / (adi_exp - 1)));

    return eq;
}

double C_p_RW (double x, double vel_R, double pres_R, double s_vel_R, double adi_exp, double time)
{
    double eq = pres_R * pow(((2 / (adi_exp + 1)) - (adi_exp - 1) * (vel_R - (x / time)) / \
                ((adi_exp + 1) * s_vel_R)), ((2 * adi_exp) / (adi_exp - 1)));

    return eq;
}
