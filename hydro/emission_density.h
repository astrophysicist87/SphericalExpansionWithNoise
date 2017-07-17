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
	return ( exp( mT * cosh(xi) - muf) / Tf ) );
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
	vector<double> result;
	result.push_back( (-muf + mT*ch_xi)/Tf2 );
	result.push_back( -mT*sh_xi/Tf );
	result.push_back( 1.0/Tf );

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

	vector<double> row1;
	row1.push_back( ((muf - mT*ch_xi)*(muf + 2.0*Tf - mT*ch_xi))/Tf4 );
	row1.push_back( (mT*(muf + Tf - mT*ch_xi)*sh_xi)/Tf3 );
	row1.push_back( -((muf + Tf - mT*ch_xi)/Tf3) );
	result.push_back( row1 );

	vector<double> row2;
	row2.push_back( (mT*(muf + Tf - mT*ch_xi)*sh_xi)/Tf3 );
	row2.push_back( (mT*(-(Tf*ch_xi) + mT*sh_xi*sh_xi))/Tf2 );
	row2.push_back( -((mT*sh_xi)/Tf2) );
	result.push_back( row2 );

	vector<double> row3;
	row3.push_back( -((muf + Tf - mT*ch_xi)/Tf3) );
	row3.push_back( -((mT*sh_xi)/Tf2) );
	row3.push_back( 1.0/Tf2 );
	result.push_back( row3 );

	return (result);
}

vector<double> void set_SB_first_derivatives_vector(double u)
{
	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double Tf2 = Tf*Tf;
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double local_SBa = S_B_a(pT*sh_xi/Tf, mT*ch_xi, pT*sh_xi);

	double dz_dv0 = -((pT*sh_xi)/Tf2);
	double dz_dv1 = (pT*ch_xi)/Tf;

	vector<double> result;
	result.push_back( SBa*dz_dv0 );
	result.push_back( SBa*dz_dv1 );
	result.push_back( 0.0 );	//S_B is independent of delta mu

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

	double local_SBa = S_B_a(pT*sh_xi/Tf, mT*ch_xi, pT*sh_xi);
	double local_SBb = S_B_b(pT*sh_xi/Tf, mT*ch_xi, pT*sh_xi);

	double dz_dv0 = -((pT*sh_xi)/Tf2);
	double dz_dv1 = (pT*ch_xi)/Tf;

	double dz_dv0_dv0 = (2*pT*sh_xi)/Tf3;
	double dz_dv0_dv1 = -((pT*ch_xi)/Tf2);
	double dz_dv1_dv1 = (pT*sh_xi)/Tf`;

	vector<double> row1;
	row1.push_back( SBa*dz_dv0_dv0 + SBb*dz_dv0*dz_dv0 );
	row1.push_back( SBa*dz_dv0_dv1 + SBb*dz_dv0*dz_dv1 );
	row1.push_back( 0.0 );
	result.push_back( row1 );

	vector<double> row2;
	row2.push_back( SBa*dz_dv0_dv1 + SBb*dz_dv0*dz_dv1 );
	row2.push_back( SBa*dz_dv1_dv1 + SBb*dz_dv1*dz_dv1 );
	row2.push_back( 0.0 );
	result.push_back( row2 );

	vector<double> row3;
	row3.push_back( 0.0 );
	row3.push_back( 0.0 );
	row3.push_back( 0.0 );
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

inline void set_Q_X_k_x(int iX, int ik, int iu)
{
	QXkx[iX][ik][iu] = 
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
