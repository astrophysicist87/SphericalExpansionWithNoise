#ifndef EMISSION_DENSITY_H
#define EMISSION_DENSITY_H

/////////////////////////////////////////////////////
//**********************************
// General outline and methodology:
//**********************************
// Emission density S is first expanded
// schematically into series of fluctuations:
//		S = S^0 + delta_k * S_k^1
//			+ 0.5 * delta_k * delta_l * S_{kl}^2 + ...
//
// Each term depends on rapidity and
// thermodynamic quantities at freeze-out
//
// S_k^1 and S_{kl}^2 built up out of
// derivatives of S(x)
//
// Defined S == S_A * S_B to simplify derivatives
//
// All derivatives taken explicitly beforehand,
// and then written out in terms of appropriate
// functions.  All derivative coefficients evaluated
// at zero fluctuations.
//
// delta components are indexed as follows:
// delta T - 0, omega - 1, delta mu - 2
/////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "lib.h"
#include "defs1.h"

////////////////////////////////////////////////////////////////////
// Psi function stuff
inline double Psi_A(double xi)
{
	double prefactor = -ds*tauf*tauf*tauf / (8.0*M_PI*M_PI*M_PI);
	double shxi = sinh(xi);
	double arg = ( muf - mT*cosh(xi) ) / Tf;	//<=0
	return (
			prefactor * shxi * shxi * exp( arg )
			);
}

inline double Psi_B0(double a1, double a2, double z)
{
	return (
			4.0*M_PI*( ( a1 * z + a2 ) * sinh(z) - a2 * z * cosh(z) ) / (z*z)
			);
}

inline double Psi_B1(double a1, double a2, double z)
{
	double I0 = bessel_I0(0.5*z);
	double I1 = bessel_I1(0.5*z);
	return (
			M_PI*M_PI*( 2.0*(a1*z+a2)I1*I0 - a2*z*(I0*I0+I1*I1) ) / z
			);
}

inline double Psi_B2(double a1, double a2, double z)
{
	return (
			4.0*M_PI*( (a1*z+2.0*a2)*(z*sinh(z)+1)
						- (a1*z + a2*(z*z+2.0)) * cosh(z) ) / (z*z*z)
			);
}

inline double Psi_0(double xi)
{
	double a1 = mT*cosh(xi);
	double a2 = pT*sinh(xi);
	return ( Psi_A(xi) * Psi_B0(a1, a2, a2/Tf) );
}

inline double Psi_1(double xi)
{
	double a1 = mT*cosh(xi);
	double a2 = pT*sinh(xi);
	return ( Psi_A(xi) * Psi_B1(a1, a2, a2/Tf) );
}

inline double Psi_2(double xi)
{
	double a1 = mT*cosh(xi);
	double a2 = pT*sinh(xi);
	return ( Psi_A(xi) * Psi_B2(a1, a2, a2/Tf) );
}

//Functions useful for computing derivatives of Psi functions
inline double Psi_a_B0(double a1, double a2, double z)
{
	return (
				-4.0*M_PI * ( (a1*z + a2*(z*z+2.0))*sinh(z) - z*cosh(z)*(a1*z+2.0*a2) ) / (z*z*z)
			);
}

inline double Psi_b_B0(double a1, double a2, double z)
{
	return (
				4.0*M_PI * ( (z*z+2.0)*(a1*z+3.0*a2)*sinh(z) - z*cosh(z)*(2.0*a1*z+a2*(z*z+6.0)) ) / (z*z*z*z)
			);
}

inline double Psi_a_B1(double a1, double a2, double z)
{
	double I0 = bessel_I0(0.5*z);
	double I1 = bessel_I1(0.5*z);
	double I2 = bessel_I2(0.5*z);
	return (
				M_PI*M_PI*( z*I0*I0*(a1*z+a2) + z*I1*I1*(a1*z+3.0*a2)
							-2.0*I1*I0*(a1*z+a2* ( z*z+2.0 ) ) ) / (z*z)
			);
}

inline double Psi_b_B1(double a1, double a2, double z)
{
	return (
				-M_PI*M_PI*( z*I0*I0*(a1*z+a2*(z*z+3.0)) + z*I1*I1*(3.0*a1*z+a2*(z*z+11.0))
							 - 2.0*(z*z+2.0) * I1*I0*(a1*z+3*a2) ) / (z*z*z)
			);
}

inline double Psi_a_B2(double a1, double a2, double z)
{
	return (
				4.0*M_PI * ( ( (z*z+2.0)*cosh(z) - 2.0 )*(a1*z+3.0*a2)
							- z*sinh(z)*(2.0*a1*z + a2*(z*z+6.0)) ) / (z*z*z*z)
			);
}

inline double Psi_b_B2(double a1, double a2, double z)
{
	return (
				4.0*M_PI * ( ( z*sinh(z)*(z*z+6.0) + 6.0 )*(a1*z+4.0*a2)
							- cosh(z) * ( 3.0*a1*z*(z*z+2.0) + a2*(z*z*z*z + 12.0*z*z + 24.0) ) ) / (z*z*z*z*z)
			);
}

////////////////////////////////////////////////////////////////////

void set_SA_first_derivatives_vector(int iu)
{
	double u = u_pts[iu];

	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double prefactor = sh_xi*sh_xi*exp((muf - mT*u)/Tf);
	double lambda = 1.0 + (muf - mT*u)/Tf;
	double s_tilde = sf/Tf;

	dSA_dX[iu][0] = -prefactor*s_tilde*((lambda - 1.0)*chi_tilde_mu_mu + chi_tilde_T_mu);
	dSA_dX[iu][1] = -prefactor*mT*sh_xi/Tf;
	dSA_dX[iu][2] = prefactor*s_tilde*((lambda - 1.0)*chi_tilde_T_mu + chi_tilde_T_T);

	return;
}

void set_SA_second_derivatives_array(int iu)
{
	double u = u_pts[iu];

	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double prefactor = sh_xi*sh_xi*exp((muf - mT*u)/Tf);
	double lambda = 1.0 + (muf - mT*u)/Tf;

	dSA_dX_dY[iu][0][0] = prefactor*s_tilde*s_tilde*( (lambda*lambda-1.0)*chi_tilde_mu_mu*chi_tilde_mu_mu
													+ 2.0*chi_tilde_mu_mu*chi_tilde_T_mu
													+ chi_tilde_T_mu*chi_tilde_T_mu );
	dSA_dX_dY[iu][0][1] = prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu);
	dSA_dX_dY[iu][0][2] = -prefactor*s_tilde*s_tilde*(
						chi_tilde_T_mu*((lambda*lambda - 1.0)*chi_tilde_mu_mu + lambda*chi_tilde_T_mu)
						 + chi_tilde_T_T*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu)
					);

	dSA_dX_dY[iu][1][0] = prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu);
	dSA_dX_dY[iu][1][1] = prefactor*(mT/Tf)*((mT/Tf)*sh_xi*sh_xi-ch_xi);
	dSA_dX_dY[iu][1][2] = -prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_T_mu+chi_tilde_T_T);

	dSA_dX_dY[iu][2][0] = -prefactor*s_tilde*s_tilde*(
						chi_tilde_T_mu*((lambda*lambda - 1.0)*chi_tilde_mu_mu + lambda*chi_tilde_T_mu)
						 + chi_tilde_T_T*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu)
					);
	dSA_dX_dY[iu][2][1] = -prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_T_mu+chi_tilde_T_T);
	dSA_dX_dY[iu][2][2] = prefactor*s_tilde*s_tilde*( (lambda*lambda-1.0)*chi_tilde_T_mu*chi_tilde_T_mu
													+ 2.0*chi_tilde_T_mu*chi_tilde_T_T
													+ chi_tilde_T_T*chi_tilde_T_T );

	return;
}

void set_SB_first_derivatives_vector(int iu)
{
	double u = u_pts[iu];

	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double Tf2 = Tf*Tf;
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;

	double z_loc = pT*sh_xi/Tf;
	double local_SBa = S_B_a(z_loc, mT*ch_xi, pT*sh_xi);
	double kappa_S = (pT/Tf)*sh_xi, kappa_C = (pT/Tf)*ch_xi;

	double dz_dv0 = -s_tilde*kappa_S*chi_tilde_mu_mu;
	double dz_dv1 = kappa_C;
	double dz_dv2 = s_tilde*kappa_S*chi_tilde_T_mu;

	dSB_dX[iu][0] = local_SBa*dz_dv0;
	dSB_dX[iu][1] = local_SBa*dz_dv1;
	dSB_dX[iu][2] = local_SBa*dz_dv2;

	return ;
}

void set_SB_second_derivatives_array(int iu)
{
	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	//double sh_xi = sinh(xi), ch_xi = cosh(xi);
	double u = u_pts[iu];
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;

	double z_loc = pT*sh_xi/Tf;

	double local_SBa = S_B_a(z_loc, mT*ch_xi, pT*sh_xi);
	double local_SBb = S_B_b(z_loc, mT*ch_xi, pT*sh_xi);

	double dz_dv0 = -s_tilde*kappa_S*chi_tilde_mu_mu;
	double dz_dv1 = kappa_C;
	double dz_dv2 = s_tilde*kappa_S*chi_tilde_T_mu;

	double dz_dv0_dv0 = 2.0*s_tilde*s_tilde*kappa_S*chi_tilde_mu_mu*chi_tilde_mu_mu;
	double dz_dv0_dv1 = -s_tilde*kappa_C*chi_tilde_mu_mu;
	double dz_dv0_dv2 = -2.0*s_tilde*s_tilde*kappa_S*chi_tilde_T_mu*chi_tilde_mu_mu;
	double dz_dv1_dv0 = dz_dv0_dv1;
	double dz_dv1_dv1 = kappa_S;
	double dz_dv1_dv2 = s_tilde*kappa_C*chi_tilde_T_mu;
	double dz_dv2_dv0 = dz_dv0_dv2;
	double dz_dv2_dv1 = dz_dv1_dv2;
	double dz_dv2_dv2 = 2.0*s_tilde*s_tilde*kappa_S*chi_tilde_T_mu*chi_tilde_T_mu;

	dSB_dX_dY[iu][0][0] = local_SBa*dz_dv0_dv0 + local_SBb*dz_dv0*dz_dv0;
	dSB_dX_dY[iu][0][1] = local_SBa*dz_dv0_dv1 + local_SBb*dz_dv0*dz_dv1;
	dSB_dX_dY[iu][0][2] = local_SBa*dz_dv0_dv2 + local_SBb*dz_dv0*dz_dv2;

	dSB_dX_dY[iu][1][0] = local_SBa*dz_dv1_dv0 + local_SBb*dz_dv1*dz_dv0;
	dSB_dX_dY[iu][1][1] = local_SBa*dz_dv1_dv1 + local_SBb*dz_dv1*dz_dv1;
	dSB_dX_dY[iu][1][2] = local_SBa*dz_dv1_dv2 + local_SBb*dz_dv1*dz_dv2;

	dSB_dX_dY[iu][2][0] = local_SBa*dz_dv2_dv0 + local_SBb*dz_dv2*dz_dv0;
	dSB_dX_dY[iu][2][1] = local_SBa*dz_dv2_dv1 + local_SBb*dz_dv2*dz_dv1;
	dSB_dX_dY[iu][2][2] = local_SBa*dz_dv2_dv2 + local_SBb*dz_dv2*dz_dv2;

	return;
}

inline void set_N0()
{
	N0 = 0.0;
	for (int iu = 0; iu < n_u_pts; ++iu)
		N0 += u_wts[iu] * S0x[iu];
}

inline void set_int_wij_S0x_S0xp()
{
	double result = 0.0;
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
		result += u_wts[iu] * bar_w_ij(iu, iup) * S0x[iu] * S0x[iup];
	int_wij_S0x_S0xp = result;
}

inline void set_S0(int iu)
{
	SA[iu] = S_A(iu);
	SB[iu] = S_B(iu);
	S0x[iu] = SA[iu]*SB[iu];
}

inline void set_S_1_X(int iu, int iX)
{
	S1Xx[iX][iu] = dSA_dX[iu][iX]*SB[iu] + dSB_dX[iu][iX]*SA[iu];
}

inline void set_S_2_XY(int iu, int iX, int iY)
{
	S2XYx[iX][iY][iu] = 
			dSA_dX_dY[iu][iX][iY] * SB[iu]
			+ dSA_dX[iu][iX] * dSB_dX[iu][iY]
			+ dSB_dX[iu][iX] * dSA_dX[iu][iY]
			+ dSB_dX_dY[iu][iX][iY] * SA[iu]
			;
}

//smearing functions here
inline void set_theta_0_XY(int iu, int iX, int iY)
{
	double result = 0.0;

	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		result += u_wts[iup] * bar_w_ij(iu, iup) * S0x[iup] * S2XYx[iX][iY][iu];
//cout << "theta0 " << u_wts[iup] << "   " << bar_w_ij(iu, iup) << "   " << S0x[iup] << "   " << S2XYx[iX][iY][iu] << endl;
	}

	result -= int_wij_S0x_S0xp * S2XYx[iX][iY][iu] / N0;
//cout << "theta0 " << int_wij_S0x_S0xp << "   " << S2XYx[iX][iY][iu] << "   " << N0 << endl;
	theta0XY[iX][iY][iu] = result;
}


inline void set_theta_1_XY(int iu, int iup, int iX, int iY)
{
	double tmp1 = bar_w_ij(iu, iup) + 3.0 * int_wij_S0x_S0xp / (N0*N0);
	
	double tmp2 = 0.0;
	for (int iupp = 0; iupp < n_u_pts; ++iupp)
		tmp2 += u_wts[iupp] * bar_w_ij(iup, iupp) * S0x[iupp];

	theta1XY[iX][iY][iu][iup] = S1Xx[iX][iu] * S1Xx[iY][iup] * ( tmp1 - 4.0 * tmp2 / N0 );
}


// End of file

#endif
