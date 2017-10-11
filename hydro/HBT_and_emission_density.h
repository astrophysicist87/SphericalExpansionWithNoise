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
#include "Phi_functions_and_derivatives.h"

inline void set_N_00_ij()
{
	N_00_00 = 0.0;
	N_00_ss = 0.0;
	N_00_oo = 0.0;
	N_00_ot = 0.0;
	N_00_tt = 0.0;

	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		double u_loc = u_pts[iu];
		double up_loc = u_pts[iup];
		double measure = 1.0 / sqrt( (u_loc*u_loc-1.0) * (up_loc*up_loc-1.0) );	//for proper u-integration
		N_00_00 += u_wts[iu] * u_wts[iup] * measure * Phi_0_0[iu][iup];
		N_00_ss += u_wts[iu] * u_wts[iup] * measure * Phi_s_s[iu][iup];
		N_00_oo += u_wts[iu] * u_wts[iup] * measure * Phi_o_o[iu][iup];
		N_00_ot += u_wts[iu] * u_wts[iup] * measure * Phi_o_t[iu][iup];
		N_00_tt += u_wts[iu] * u_wts[iup] * measure * Phi_t_t[iu][iup];
	}

	return;
}



////////////////////////////////////////////////////////////////

//smearing functions here
inline void set_theta_0_ij_XY()
{
	/*if (1)
	{
		cerr << "This function not finished!  E.g., u-integral needs proper integration measure, etc.!!" << endl;
		exit(1);
	}*/

	for (int iX = 0; iX < 3; ++iX)
	for (int iY = 0; iY < 3; ++iY)
	for (int iu = 0; iu < n_u_pts; ++iu)
	{
		theta_0_ss_XY[iX][iY][iu] = int_dup_d_Phi_s_s_dX_dY[iX][iY][iu] / N_00_00
									- (N_00_ss / ( N_00_00*N_00_00 )) * int_dup_d_Phi_0_0_dX_dY[iX][iY][iu];
		theta_0_oo_XY[iX][iY][iu] = int_dup_d_Phi_o_o_dX_dY[iX][iY][iu] / N_00_00
									- (N_00_oo / ( N_00_00*N_00_00 )) * int_dup_d_Phi_0_0_dX_dY[iX][iY][iu];
		theta_0_ot_XY[iX][iY][iu] = int_dup_d_Phi_o_t_dX_dY[iX][iY][iu] / N_00_00
									- (N_00_ot / ( N_00_00*N_00_00 )) * int_dup_d_Phi_0_0_dX_dY[iX][iY][iu];
		theta_0_tt_XY[iX][iY][iu] = int_dup_d_Phi_t_t_dX_dY[iX][iY][iu] / N_00_00
									- (N_00_tt / ( N_00_00*N_00_00 )) * int_dup_d_Phi_0_0_dX_dY[iX][iY][iu];
	}

	return;
}


inline void set_theta_1_ij_XY()
{
	/*if (1)
	{
		cerr << "This function not finished!  E.g., u-integral needs proper integration measure, etc.!!" << endl;
		exit(1);
	}*/

	for (int iX = 0; iX < 3; ++iX)
	for (int iYp = 0; iYp < 3; ++iYp)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		theta_1_ss_XY[iX][iYp][iu][iup] = d_Phi_s_s_dX_dY[iX][iYp][iu][iup] / N_00_00
										- (N_00_ss / ( N_00_00*N_00_00 )) * d_Phi_0_0_dX_dY[iX][iYp][iu][iup]
										- 4.0 * ( int_dup_d_Phi_0_0_dX[iX][iu] / ( N_00_00*N_00_00 ) )
										* (
											int_dup_d_Phi_s_s_dX[iYp][iup]
												 - ( N_00_ss / N_00_00 ) * int_dup_d_Phi_0_0_dX[iYp][iup]
											);
		theta_1_oo_XY[iX][iYp][iu][iup] = d_Phi_o_o_dX_dY[iX][iYp][iu][iup] / N_00_00
										- (N_00_oo / ( N_00_00*N_00_00 )) * d_Phi_0_0_dX_dY[iX][iYp][iu][iup]
										- 4.0 * ( int_dup_d_Phi_0_0_dX[iX][iu] / ( N_00_00*N_00_00 ) )
										* (
											int_dup_d_Phi_o_o_dX[iYp][iup]
												 - ( N_00_oo / N_00_00 ) * int_dup_d_Phi_0_0_dX[iYp][iup]
											);
		theta_1_ot_XY[iX][iYp][iu][iup] = d_Phi_o_t_dX_dY[iX][iYp][iu][iup] / N_00_00
										- (N_00_ot / ( N_00_00*N_00_00 )) * d_Phi_0_0_dX_dY[iX][iYp][iu][iup]
										- 4.0 * ( int_dup_d_Phi_0_0_dX[iX][iu] / ( N_00_00*N_00_00 ) )
										* (
											int_dup_d_Phi_o_t_dX[iYp][iup]
												 - ( N_00_ot / N_00_00 ) * int_dup_d_Phi_0_0_dX[iYp][iup]
											);
		theta_1_tt_XY[iX][iYp][iu][iup] = d_Phi_t_t_dX_dY[iX][iYp][iu][iup] / N_00_00
										- (N_00_tt / ( N_00_00*N_00_00 )) * d_Phi_0_0_dX_dY[iX][iYp][iu][iup]
										- 4.0 * ( int_dup_d_Phi_0_0_dX[iX][iu] / ( N_00_00*N_00_00 ) )
										* (
											int_dup_d_Phi_t_t_dX[iYp][iup]
												 - ( N_00_tt / N_00_00 ) * int_dup_d_Phi_0_0_dX[iYp][iup]
											);
	}

	return;
}


// End of file

#endif
