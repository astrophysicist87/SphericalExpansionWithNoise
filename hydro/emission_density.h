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

/*USAGE: debugger(__LINE__, __FILE__);*/
void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

inline double S_A(double xi)
{
	double sh_xi = sinh(xi);
	return ( sh_xi*sh_xi*exp( mT * cosh(xi) - muf) / Tf ) );
}

inline double S_B(double xi)
{
	double a1 = mT*ch_xi;
	double sh_xi = sinh(xi);
	double a2 = pT*sh_xi;

	double z = pT*sh_xi/Tf;
	double sh_z = sinh(z);
	double ch_z = cosh(z);

	return ( a1*sh_z/z - a2*( sh_z - z*ch_z ) / (z*z) );
}

// need these next two functions for
// taking derivatives of S_B
inline double S_B_a(double z, double a1, double a2)
{
	double z2 = z*z;
	double z3 = z2*z;

	return (
			( z*( a1*z-2.0*a2 )*cosh(z)
				+ ( a2*( z2+2.0 )-a1*z )*sinh(z) )
			/ z3 );
}

inline double S_B_b(double z, double a1, double a2)
{
	double z2 = z*z;
	double z3 = z2*z;
	double z4 = z3*z;

	return (
				z*(a2*(z2+6.0)-2.0*a1*z)*cosh(z)
					+ (z2+2.0)*(a1*z-3.0*a2)*sinh(z)
			);
}

vector<double> void set_SA_first_derivatives_vector(double u)
{
	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double prefactor = sh_xi*sh_xi*exp((muf - mT*u)/Tf);
	double lambda = 1.0 + (muf - mT*u)/Tf;
	double s_tilde = sf/Tf;

	vector<double> result;
	result.push_back( -prefactor*s_tilde*((lambda - 1.0)*chi_tilde_mu_mu + chi_tilde_T_mu) );
	result.push_back( -prefactor*mT*sh_xi/Tf );
	result.push_back( prefactor*s_tilde*((lambda - 1.0)*chi_tilde_T_mu + chi_tilde_T_T) );

	return (result);
}

vector<vector<double> > void set_SA_second_derivatives_array(double xi)
{
	vector<vector<double> > result;

	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double prefactor = sh_xi*sh_xi*exp((muf - mT*u)/Tf);
	double lambda = 1.0 + (muf - mT*u)/Tf;

	vector<double> row1;
	row1.push_back( prefactor*s_tilde*s_tilde*( (lambda*lambda-1.0)*chi_tilde_mu_mu*chi_tilde_mu_mu
													+ 2.0*chi_tilde_mu_mu*chi_tilde_T_mu
													+ chi_tilde_T_mu*chi_tilde_T_mu ) );
	row1.push_back( prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu) );
	row1.push_back( -prefactor*s_tilde*s_tilde*(
						chi_tilde_T_mu*((lambda*lambda - 1.0)*chi_tilde_mu_mu + lambda*chi_tilde_T_mu)
						 + chi_tidle_T_T*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu)
					) );
	result.push_back( row1 );

	vector<double> row2;
	row2.push_back( prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu) );
	row2.push_back( prefactor*(mT/Tf)*((mT/Tf)*sh_xi*sh_xi-ch_xi) );
	row2.push_back( -prefactor*(mT/Tf)*s_tilde*sh_xi(lambda*chi_tilde_T_mu+chi_tilde_T_T) );
	result.push_back( row2 );

	vector<double> row3;
	row3.push_back( -prefactor*s_tilde*s_tilde*(
						chi_tilde_T_mu*((lambda*lambda - 1.0)*chi_tilde_mu_mu + lambda*chi_tilde_T_mu)
						 + chi_tidle_T_T*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu)
					) );
	row3.push_back( -prefactor*(mT/Tf)*s_tilde*sh_xi(lambda*chi_tilde_T_mu+chi_tilde_T_T) );
	row3.push_back( prefactor*s_tilde*s_tilde*( (lambda*lambda-1.0)*chi_tilde_T_mu*chi_tilde_T_mu
													+ 2.0*chi_tilde_T_mu*chi_tilde_T_T
													+ chi_tilde_T_T*chi_tilde_T_T ) );
	result.push_back( row3 );

	return (result);
}

vector<double> void set_SB_first_derivatives_vector(double u)
{
	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double Tf2 = Tf*Tf;
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;

	double z_loc = pT*sh_xi/Tf;
	double local_SBa = S_B_a(z_loc, mT*ch_xi, pT*sh_xi);
	double kappa_S = (pT/Tf)*sh_xi, kappa_C = (pT/Tf)*ch_xi;

	double dz_dv0 = -s_tilde*kappa_S*chi_tilde_mu_mu;
	double dz_dv1 = kappa_C;
	double dz_dv2 = s_tilde*kappa_S*chi_tilde_T_mu;

	vector<double> result;
	result.push_back( SBa*dz_dv0 );
	result.push_back( SBa*dz_dv1 );
	result.push_back( SBa*dz_dv2 );

	return (result);
}

vector<vector<double> > void set_SB_second_derivatives_array(double u)
{
	vector<vector<double> > result;

	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	//double sh_xi = sinh(xi), ch_xi = cosh(xi);
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

	vector<double> row1;
	row1.push_back( SBa*dz_dv0_dv0 + SBb*dz_dv0*dz_dv0 );
	row1.push_back( SBa*dz_dv0_dv1 + SBb*dz_dv0*dz_dv1 );
	row1.push_back( SBa*dz_dv0_dv2 + SBb*dz_dv0*dz_dv2 );
	result.push_back( row1 );

	vector<double> row2;
	row2.push_back( SBa*dz_dv1_dv0 + SBb*dz_dv1*dz_dv0 );
	row2.push_back( SBa*dz_dv1_dv1 + SBb*dz_dv1*dz_dv1 );
	row2.push_back( SBa*dz_dv1_dv2 + SBb*dz_dv1*dz_dv2 );
	result.push_back( row2 );

	vector<double> row3;
	row3.push_back( SBa*dz_dv2_dv0 + SBb*dz_dv2*dz_dv0 );
	row3.push_back( SBa*dz_dv2_dv1 + SBb*dz_dv2*dz_dv1 );
	row3.push_back( SBa*dz_dv2_dv2 + SBb*dz_dv2*dz_dv2 );
	result.push_back( row3 );

	return (result);
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
	S0x[iu] = S_A[iu]*S_B[iu] );
}

inline void set_S_1_X(int iu, int iX)
{
	S1Xx[iX][iu] = dSA_dvi[iu][iX]*SB[iu] + dSB_dvi[iu][iX]*SA[iu];
}

inline void set_S_2_XY(int iu, int iX, int iY)
{
	S2XYx[iX][iY][iu] = 
			dSA_dX_dY[iu][iX][iY] * SB[iu]
			+ dSA_dX[iu][iX] * dSB_dvi[iu][iY]
			+ dSB_dX[iu][iX] * dSA_dvi[iu][iY]
			+ dSB_dX_dY[iu][iX][iY] * SA[iu]
			;
}

//smearing functions here
inline void set_theta_0_XY(int iu, int iX, int iY)
{
	double result = 0.0;
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		result += u_wts[iup] * bar_w_ij(iu, iup) * S0x[iup] * S2XYx[iX][iY][iu]
	}
	result -= int_wij_S0x_S0xp * S2XYx[iX][iY][iu] / N0;
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
