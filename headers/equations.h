#ifndef EQ_INCLUDED
#define EQ_INCLUDED

double A_v_xL_star (double vel_L, double s_vel_L, double pres_L, double pres_unk, double adi_exp);
double A_rho_L_star (double dens_L, double pres_L, double pres_unk, double adi_exp);
double A_v_xR_star (double dens_R, double vel_R, double s_vel_R, double pres_R, double pres_unk, double adi_exp);
double A_rho_R_star (double dens_R, double pres_R, double pres_unk, double adi_exp);
double A_D (double vel_R, double pres_R, double s_vel_R, double pres_unk, double adi_exp);
double A1_v_xL_star (double dens_L, double vel_L, double s_vel_L, double pres_L, double pres_unk, double adi_exp);
double A1_rho_L_star (double dens_L, double pres_L, double pres_unk, double adi_exp);
double A1_v_xR_star (double vel_R, double s_vel_R, double pres_R, double pres_unk, double adi_exp);
double A1_rho_R_star (double dens_R, double pres_R, double pres_unk, double adi_exp);
double A1_D (double vel_L, double pres_L, double s_vel_L, double pres_unk, double adi_exp);
double B_v_xL_star (double dens_L, double vel_L, double s_vel_L, double pres_L, double pres_unk, double adi_exp);
double B_v_xR_star (double dens_R, double vel_R, double s_vel_R, double pres_R, double pres_unk, double adi_exp);
double C_v_xL_star (double vel_L, double s_vel_L, double pres_L, double pres_unk, double adi_exp);
double C_v_xR_star (double vel_R, double s_vel_R, double pres_R, double pres_unk, double adi_exp);

double A_dv (double pres_L, double pres_unk, double s_vel_L, double dens_R, double pres_R, double s_vel_R, double adi_exp);
double A1_dv (double dens_L, double s_vel_L, double pres_L, double pres_unk, double pres_R, double s_vel_R, double adi_exp);
double B_dv (double dens_L, double s_vel_L, double pres_L, double pres_unk, double dens_R, double s_vel_R, double pres_R, double adi_exp);
double C_dv (double s_vel_L, double pres_L, double pres_unk, double s_vel_R, double pres_R, double adi_exp);

double A_v_RW (double x, double vel_L, double s_vel_L, double adi_exp, double time);
double A_rho_RW (double x, double dens_L, double vel_L, double s_vel_L, double adi_exp, double time);
double A_p_RW (double x, double vel_L, double pres_L, double s_vel_L, double adi_exp, double time);
double A1_v_RW (double x, double vel_R, double s_vel_R, double adi_exp, double time);
double A1_rho_RW (double x, double dens_R, double vel_R, double s_vel_R, double adi_exp, double time);
double A1_p_RW (double x, double vel_R, double pres_R, double s_vel_R, double adi_exp, double time);
double C_v_RW (double x, double vel_R, double s_vel_R, double adi_exp, double time);
double C_rho_RW (double x, double dens_R, double vel_R, double s_vel_R, double adi_exp, double time);
double C_p_RW (double x, double vel_R, double pres_R, double s_vel_R, double adi_exp, double time);

#endif