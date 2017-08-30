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
/// All derivatives taken explicitly beforehand,
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
// Phi functions which use above Psi functions
////////////////////////////////////////////////////////////////////

inline void set_Phi_ij_grids()
{
	for (int ixi = 0; ixi < n_xi_pts; ++ixi)
	for (int ixip = 0; ixip < n_xi_pts; ++ixip)
	{
		double xi = xi_pts[ixi], xip = xi_pts[ixip];
		double sh_xi = sinh(xi), ch_xi = cosh(xi);
		double sh_xip = sinh(xip), ch_xip = cosh(xip);

		double Psi_0_xi = Psi_0(xi), Psi_0_xip = Psi_0(xip);
		double Psi_1_xi = Psi_1(xi), Psi_1_xip = Psi_1(xip);
		double Psi_2_xi = Psi_2(xi), Psi_2_xip = Psi_2(xip);

		Phi_0_0[ixi][ixip] = Psi_0_xi * Psi_0_xip;
		Phi_s_s[ixi][ixip] = 0.5*tauf*tauf*( ( Psi_0_xi - Psi_2_xi) * Psi_0_xip * sh_xi * sh_xi
											 + ( Psi_0_xip - Psi_2_xip) * Psi_0_xi * sh_xip * sh_xip );
		Phi_o_o[ixi][ixip] = 0.5*tauf*tauf*( Psi_2_xi * Psi_0_xip * sh_xi * sh_xi
											 + Psi_2_xip * Psi_0_xi * sh_xip * sh_xip
											 - 2.0 * Psi_1_xi * Psi_1_xip * sh_xi * sh_xip );
		Phi_o_t[ixi][ixip] = 0.5*tauf*tauf*( ch_xi - ch_xip ) * ( Psi_1_xi * Psi_0_xip * sh_xi
																	- Psi_1_xip * Psi_0_xi * sh_xip );
		Phi_t_t[ixi][ixip] = 0.5*tauf*tauf*Psi_0_xi * Psi_0_xip * ( ch_xi - ch_xip ) * ( ch_xi - ch_xip );
	}

	return;
}


////////////////////////////////////////////////////////////////////

void set_PsiA_first_derivatives_vector(int iu)
{
	double u = u_pts[iu];

	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double prefactor = ds*tauf*tauf*tauf*sh_xi*sh_xi*exp((muf - mT*u)/Tf) / (8.0*M_PI*M_PI*M_PI);
	double lambda = 1.0 + (muf - mT*u)/Tf;
	double s_tilde = sf/Tf;

	dPsiA_dX[iu][0] = -prefactor*s_tilde*((lambda - 1.0)*chi_tilde_mu_mu + chi_tilde_T_mu);
	dPsiA_dX[iu][1] = -prefactor*mT*sh_xi/Tf;
	dPsiA_dX[iu][2] = prefactor*s_tilde*((lambda - 1.0)*chi_tilde_T_mu + chi_tilde_T_T);

	return;
}

void set_PsiA_second_derivatives_array(int iu)
{
	double u = u_pts[iu];

	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double prefactor = ds*tauf*tauf*tauf*sh_xi*sh_xi*exp((muf - mT*u)/Tf) / (8.0*M_PI*M_PI*M_PI);
	double lambda = 1.0 + (muf - mT*u)/Tf;

	dPsiA_dX_dY[iu][0][0] = prefactor*s_tilde*s_tilde*( (lambda*lambda-1.0)*chi_tilde_mu_mu*chi_tilde_mu_mu
													+ 2.0*chi_tilde_mu_mu*chi_tilde_T_mu
													+ chi_tilde_T_mu*chi_tilde_T_mu );
	dPsiA_dX_dY[iu][0][1] = prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu);
	dPsiA_dX_dY[iu][0][2] = -prefactor*s_tilde*s_tilde*(
						chi_tilde_T_mu*((lambda*lambda - 1.0)*chi_tilde_mu_mu + lambda*chi_tilde_T_mu)
						 + chi_tilde_T_T*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu)
					);

	dPsiA_dX_dY[iu][1][0] = prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu);
	dPsiA_dX_dY[iu][1][1] = prefactor*(mT/Tf)*((mT/Tf)*sh_xi*sh_xi-ch_xi);
	dPsiA_dX_dY[iu][1][2] = -prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_T_mu+chi_tilde_T_T);

	dPsiA_dX_dY[iu][2][0] = -prefactor*s_tilde*s_tilde*(
						chi_tilde_T_mu*((lambda*lambda - 1.0)*chi_tilde_mu_mu + lambda*chi_tilde_T_mu)
						 + chi_tilde_T_T*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu)
					);
	dPsiA_dX_dY[iu][2][1] = -prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_T_mu+chi_tilde_T_T);
	dPsiA_dX_dY[iu][2][2] = prefactor*s_tilde*s_tilde*( (lambda*lambda-1.0)*chi_tilde_T_mu*chi_tilde_T_mu
													+ 2.0*chi_tilde_T_mu*chi_tilde_T_T
													+ chi_tilde_T_T*chi_tilde_T_T );

	return;
}

void set_PsiBk_first_derivatives_vector(int iu)
{
	double u = u_pts[iu];

	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double Tf2 = Tf*Tf;
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;

	double z_loc = pT*sh_xi/Tf;
	double local_PsiB0a = Psi_a_B0(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB1a = Psi_a_B1(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB2a = Psi_a_B2(z_loc, mT*ch_xi, pT*sh_xi);
	double kappa_S = (pT/Tf)*sh_xi, kappa_C = (pT/Tf)*ch_xi;

	double dz_dv0 = -s_tilde*kappa_S*chi_tilde_mu_mu;
	double dz_dv1 = kappa_C;
	double dz_dv2 = s_tilde*kappa_S*chi_tilde_T_mu;

	dPsiB0_dX[iu][0] = local_PsiB0a*dz_dv0;
	dPsiB0_dX[iu][1] = local_PsiB0a*dz_dv1;
	dPsiB0_dX[iu][2] = local_PsiB0a*dz_dv2;

	dPsiB1_dX[iu][0] = local_PsiB1a*dz_dv0;
	dPsiB1_dX[iu][1] = local_PsiB1a*dz_dv1;
	dPsiB1_dX[iu][2] = local_PsiB1a*dz_dv2;

	dPsiB2_dX[iu][0] = local_PsiB2a*dz_dv0;
	dPsiB2_dX[iu][1] = local_PsiB2a*dz_dv1;
	dPsiB2_dX[iu][2] = local_PsiB2a*dz_dv2;

	return ;
}

void set_PsiBk_second_derivatives_array(int iu)
{
	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	//double sh_xi = sinh(xi), ch_xi = cosh(xi);
	double u = u_pts[iu];
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;

	double z_loc = pT*sh_xi/Tf;

	double local_PsiB0a = Psi_a_B0(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB0b = Psi_b_B0(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB1a = Psi_a_B1(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB1b = Psi_b_B1(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB2a = Psi_a_B2(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB2b = Psi_b_B2(z_loc, mT*ch_xi, pT*sh_xi);

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

	//second derivatives of Psi_B0
	dPsiB0_dX_dY[iu][0][0] = local_PsiB0a*dz_dv0_dv0 + local_PsiB0b*dz_dv0*dz_dv0;
	dPsiB0_dX_dY[iu][0][1] = local_PsiB0a*dz_dv0_dv1 + local_PsiB0b*dz_dv0*dz_dv1;
	dPsiB0_dX_dY[iu][0][2] = local_PsiB0a*dz_dv0_dv2 + local_PsiB0b*dz_dv0*dz_dv2;

	dPsiB0_dX_dY[iu][1][0] = local_PsiB0a*dz_dv1_dv0 + local_PsiB0b*dz_dv1*dz_dv0;
	dPsiB0_dX_dY[iu][1][1] = local_PsiB0a*dz_dv1_dv1 + local_PsiB0b*dz_dv1*dz_dv1;
	dPsiB0_dX_dY[iu][1][2] = local_PsiB0a*dz_dv1_dv2 + local_PsiB0b*dz_dv1*dz_dv2;

	dPsiB0_dX_dY[iu][2][0] = local_PsiB0a*dz_dv2_dv0 + local_PsiB0b*dz_dv2*dz_dv0;
	dPsiB0_dX_dY[iu][2][1] = local_PsiB0a*dz_dv2_dv1 + local_PsiB0b*dz_dv2*dz_dv1;
	dPsiB0_dX_dY[iu][2][2] = local_PsiB0a*dz_dv2_dv2 + local_PsiB0b*dz_dv2*dz_dv2;

	//second derivatives of Psi_B1
	dPsiB1_dX_dY[iu][0][0] = local_PsiB1a*dz_dv0_dv0 + local_PsiB1b*dz_dv0*dz_dv0;
	dPsiB1_dX_dY[iu][0][1] = local_PsiB1a*dz_dv0_dv1 + local_PsiB1b*dz_dv0*dz_dv1;
	dPsiB1_dX_dY[iu][0][2] = local_PsiB1a*dz_dv0_dv2 + local_PsiB1b*dz_dv0*dz_dv2;

	dPsiB1_dX_dY[iu][1][0] = local_PsiB1a*dz_dv1_dv0 + local_PsiB1b*dz_dv1*dz_dv0;
	dPsiB1_dX_dY[iu][1][1] = local_PsiB1a*dz_dv1_dv1 + local_PsiB1b*dz_dv1*dz_dv1;
	dPsiB1_dX_dY[iu][1][2] = local_PsiB1a*dz_dv1_dv2 + local_PsiB1b*dz_dv1*dz_dv2;

	dPsiB1_dX_dY[iu][2][0] = local_PsiB1a*dz_dv2_dv0 + local_PsiB1b*dz_dv2*dz_dv0;
	dPsiB1_dX_dY[iu][2][1] = local_PsiB1a*dz_dv2_dv1 + local_PsiB1b*dz_dv2*dz_dv1;
	dPsiB1_dX_dY[iu][2][2] = local_PsiB1a*dz_dv2_dv2 + local_PsiB1b*dz_dv2*dz_dv2;

	//second derivatives of Psi_B2
	dPsiB2_dX_dY[iu][0][0] = local_PsiB2a*dz_dv0_dv0 + local_PsiB2b*dz_dv0*dz_dv0;
	dPsiB2_dX_dY[iu][0][1] = local_PsiB2a*dz_dv0_dv1 + local_PsiB2b*dz_dv0*dz_dv1;
	dPsiB2_dX_dY[iu][0][2] = local_PsiB2a*dz_dv0_dv2 + local_PsiB2b*dz_dv0*dz_dv2;

	dPsiB2_dX_dY[iu][1][0] = local_PsiB2a*dz_dv1_dv0 + local_PsiB2b*dz_dv1*dz_dv0;
	dPsiB2_dX_dY[iu][1][1] = local_PsiB2a*dz_dv1_dv1 + local_PsiB2b*dz_dv1*dz_dv1;
	dPsiB2_dX_dY[iu][1][2] = local_PsiB2a*dz_dv1_dv2 + local_PsiB2b*dz_dv1*dz_dv2;

	dPsiB2_dX_dY[iu][2][0] = local_PsiB2a*dz_dv2_dv0 + local_PsiB2b*dz_dv2*dz_dv0;
	dPsiB2_dX_dY[iu][2][1] = local_PsiB2a*dz_dv2_dv1 + local_PsiB2b*dz_dv2*dz_dv1;
	dPsiB2_dX_dY[iu][2][2] = local_PsiB2a*dz_dv2_dv2 + local_PsiB2b*dz_dv2*dz_dv2;

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

inline void set_Psi_k(int iu)
{
	SA[iu] = Psi_A(iu);
	SB[iu] = S_B(iu);
	S0x[iu] = SA[iu]*SB[iu];
}

inline void set_dPsik_dX(int iu, int iX)
{
	dPsi0_dX[iX][iu] = dPsiA_dX[iu][iX]*PsiB0[iu] + dPsiB0_dX[iu][iX]*PsiA[iu];
	dPsi1_dX[iX][iu] = dPsiA_dX[iu][iX]*PsiB1[iu] + dPsiB1_dX[iu][iX]*PsiA[iu];
	dPsi2_dX[iX][iu] = dPsiA_dX[iu][iX]*PsiB2[iu] + dPsiB2_dX[iu][iX]*PsiA[iu];
}

inline void set_dPsik_dX_dY(int iu, int iX, int iY)
{
	dPsi0_dX_dY[iX][iY][iu] = 
			dPsiA_dX_dY[iu][iX][iY] * PsiB0[iu]
			+ dPsiA_dX[iu][iX] * dPsiB0_dX[iu][iY]
			+ dPsiB0_dX[iu][iX] * dPsiA_dX[iu][iY]
			+ dPsiB0_dX_dY[iu][iX][iY] * PsiA[iu]
			;

	dPsi1_dX_dY[iX][iY][iu] = 
			dPsiA_dX_dY[iu][iX][iY] * PsiB1[iu]
			+ dPsiA_dX[iu][iX] * dPsiB1_dX[iu][iY]
			+ dPsiB1_dX[iu][iX] * dPsiA_dX[iu][iY]
			+ dPsiB1_dX_dY[iu][iX][iY] * PsiA[iu]
			;

	dPsi2_dX_dY[iX][iY][iu] = 
			dPsiA_dX_dY[iu][iX][iY] * PsiB2[iu]
			+ dPsiA_dX[iu][iX] * dPsiB2_dX[iu][iY]
			+ dPsiB2_dX[iu][iX] * dPsiA_dX[iu][iY]
			+ dPsiB2_dX_dY[iu][iX][iY] * PsiA[iu]
			;
}




inline void set_d_Phi_ij_dX_grids()
{
	for (int iX = 0; iX < 3; ++iX)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		double u = u_pts[iu], up = u_pts[iup];
		double sh_xi = sqrt(u*u-1.0), ch_xi = u;
		double sh_xip = sqrt(up*up-1.0), ch_xip = up;

		double dPsi_0_dX_u = dPsi0_dX[iX][iu], Psi_0_xip = Psi_0(xip, up);
		double dPsi_1_dX_u = dPsi1_dX[iX][iu], Psi_1_xip = Psi_1(xip, up);
		double dPsi_2_dX_u = dPsi2_dX[iX][iu], Psi_2_xip = Psi_2(xip, up);

		d_Phi_0_0_dX[ixi][ixip] = dPsi_0_dX_u * Psi_0_xip;
		d_Phi_s_s_dX[ixi][ixip] = 0.5*tauf*tauf*( ( dPsi_0_dX_u - dPsi_2_dX_u) * Psi_0_xip * sh_xi * sh_xi
											 + ( Psi_0_xip - Psi_2_xip) * dPsi_0_dX_u * sh_xip * sh_xip );
		d_Phi_o_o_dX[ixi][ixip] = 0.5*tauf*tauf*( dPsi_2_dX_u * Psi_0_xip * sh_xi * sh_xi
											 + Psi_2_xip * dPsi_0_dX_u * sh_xip * sh_xip
											 - 2.0 * dPsi_1_dX_u * Psi_1_xip * sh_xi * sh_xip );
		d_Phi_o_t_dX[ixi][ixip] = 0.5*tauf*tauf*( ch_xi - ch_xip ) * ( dPsi_1_dX_u * Psi_0_xip * sh_xi
																	- Psi_1_xip * dPsi_0_dX_u * sh_xip );
		d_Phi_t_t_dX[ixi][ixip] = 0.5*tauf*tauf*dPsi_0_dX_u * Psi_0_xip * ( ch_xi - ch_xip ) * ( ch_xi - ch_xip );
	}

	return;
}

inline void set_d_Phi_ij_dX_dY_grids()
{
	for (int iX = 0; iX < 3; ++iX)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		double u = u_pts[iu], up = u_pts[iup];
		double sh_xi = sqrt(u*u-1.0), ch_xi = u;
		double sh_xip = sqrt(up*up-1.0), ch_xip = up;

		double dPsi_0_dX_dY_u = dPsi0_dX_dY[iX][iu], Psi_0_up = Psi_0(xip, up);
		double dPsi_1_dX_dY_u = dPsi1_dX_dY[iX][iu], Psi_1_up = Psi_1(xip, up);
		double dPsi_2_dX_dY_u = dPsi2_dX_dY[iX][iu], Psi_2_up = Psi_2(xip, up);

		d_Phi_0_0_dX_dY[ixi][ixip] = dPsi_0_dX_dY_u * Psi_0_up;
		d_Phi_s_s_dX_dY[ixi][ixip] = 0.5*tauf*tauf*( ( dPsi_0_dX_dY_u - dPsi_2_dX_dY_u) * Psi_0_up * sh_xi * sh_xi
											 + ( Psi_0_up - Psi_2_up) * dPsi_0_dX_dY_u * sh_xip * sh_xip );
		d_Phi_o_o_dX_dY[ixi][ixip] = 0.5*tauf*tauf*( dPsi_2_dX_dY_u * Psi_0_up * sh_xi * sh_xi
											 + Psi_2_up * dPsi_0_dX_dY_u * sh_xip * sh_xip
											 - 2.0 * dPsi_1_dX_dY_u * Psi_1_up * sh_xi * sh_xip );
		d_Phi_o_t_dX_dY[ixi][ixip] = 0.5*tauf*tauf*( ch_xi - ch_xip ) * ( dPsi_1_dX_dY_u * Psi_0_up * sh_xi
																	- Psi_1_up * dPsi_0_dX_dY_u * sh_xip );
		d_Phi_t_t_dX_dY[ixi][ixip] = 0.5*tauf*tauf*dPsi_0_dX_dY_u * Psi_0_up * ( ch_xi - ch_xip ) * ( ch_xi - ch_xip );
	}

	return;
}

//notice that this function is different from the preceding one ( dY <--> dYp )
inline void set_d_Phi_ij_dX_dYp_grids()
{
	for (int iX = 0; iX < 3; ++iX)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		double u = u_pts[iu], up = u_pts[iup];
		double sh_xi = sqrt(u*u-1.0), ch_xi = u;
		double sh_xip = sqrt(up*up-1.0), ch_xip = up;

		double dPsi_0_dX_u = dPsi0_dX[iX][iu], dPsi_0_dYp_up = dPsi0_dX[iX][iup];
		double dPsi_1_dX_u = dPsi1_dX[iX][iu], dPsi_1_dYp_up = dPsi1_dX[iX][iup];
		double dPsi_2_dX_u = dPsi2_dX[iX][iu], dPsi_2_dYp_up = dPsi2_dX[iX][iup];

		d_Phi_0_0_dX_dYp[ixi][ixip] = dPsi_0_dX_u * dPsi_0_xip;
		d_Phi_s_s_dX_dYp[ixi][ixip] = 0.5*tauf*tauf*( ( dPsi_0_dX_u - dPsi_2_dX_u) * dPsi_0_dYp_up * sh_xi * sh_xi
											 + ( dPsi_0_dYp_up - dPsi_2_dYp_up) * dPsi_0_dX_u * sh_xip * sh_xip );
		d_Phi_o_o_dX_dYp[ixi][ixip] = 0.5*tauf*tauf*( dPsi_2_dX_u * dPsi_0_dYp_up * sh_xi * sh_xi
											 + dPsi_2_xip * dPsi_0_dX_u * sh_xip * sh_xip
											 - 2.0 * dPsi_1_dX_u * dPsi_1_dYp_up * sh_xi * sh_xip );
		d_Phi_o_t_dX_dYp[ixi][ixip] = 0.5*tauf*tauf*( ch_xi - ch_xip ) * ( dPsi_1_dX_u * dPsi_0_dYp_up * sh_xi
																	- dPsi_1_dYp_up * dPsi_0_dX_u * sh_xip );
		d_Phi_t_t_dX_dYp[ixi][ixip] = 0.5*tauf*tauf*dPsi_0_dX_u * dPsi_0_dYp_up * ( ch_xi - ch_xip ) * ( ch_xi - ch_xip );
	}

	return;
}



////////////////////////////////////////////////////////////////

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
