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
		int index2D = indexer2D(iu, iup, n_u_pts, n_u_pts);
		double u_loc = u_pts[iu];
		double up_loc = u_pts[iup];
		double measure = 1.0 / sqrt( (u_loc*u_loc-1.0) * (up_loc*up_loc-1.0) );	//for proper u-integration
		N_00_00 += u_wts[iu] * u_wts[iup] * measure * Phi_0_0[index2D];
		N_00_ss += u_wts[iu] * u_wts[iup] * measure * Phi_s_s[index2D];
		N_00_oo += u_wts[iu] * u_wts[iup] * measure * Phi_o_o[index2D];
		N_00_ot += u_wts[iu] * u_wts[iup] * measure * Phi_o_t[index2D];
		N_00_tt += u_wts[iu] * u_wts[iup] * measure * Phi_t_t[index2D];
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
		int index3D = indexer3D(iX, iY, iu, 3, 3, n_u_pts);
		theta_0_ss_XY[index3D] = int_dup_d_Phi_s_s_dX_dY[index3D] / N_00_00
									- (N_00_ss / ( N_00_00*N_00_00 )) * int_dup_d_Phi_0_0_dX_dY[index3D];
		theta_0_oo_XY[index3D] = int_dup_d_Phi_o_o_dX_dY[index3D] / N_00_00
									- (N_00_oo / ( N_00_00*N_00_00 )) * int_dup_d_Phi_0_0_dX_dY[index3D];
		theta_0_ot_XY[index3D] = int_dup_d_Phi_o_t_dX_dY[index3D] / N_00_00
									- (N_00_ot / ( N_00_00*N_00_00 )) * int_dup_d_Phi_0_0_dX_dY[index3D];
		theta_0_tt_XY[index3D] = int_dup_d_Phi_t_t_dX_dY[index3D] / N_00_00
									- (N_00_tt / ( N_00_00*N_00_00 )) * int_dup_d_Phi_0_0_dX_dY[index3D];
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
		int index2Da = indexer2D(iX, iu, 3, n_u_pts);
		int index2Db = indexer2D(iYp, iup, 3, n_u_pts);
		int index4D = indexer4D(iX, iYp, iu, iup, 3, 3, n_u_pts, n_u_pts);
		theta_1_ss_XY[index4D] = d_Phi_s_s_dX_dY[index4D] / N_00_00
										- (N_00_ss / ( N_00_00*N_00_00 )) * d_Phi_0_0_dX_dY[index4D]
										- 4.0 * ( int_dup_d_Phi_0_0_dX[index2Da] / ( N_00_00*N_00_00 ) )
										* (
											int_dup_d_Phi_s_s_dX[index2Db]
												 - ( N_00_ss / N_00_00 ) * int_dup_d_Phi_0_0_dX[index2Db]
											);
		theta_1_oo_XY[index4D] = d_Phi_o_o_dX_dY[index4D] / N_00_00
										- (N_00_oo / ( N_00_00*N_00_00 )) * d_Phi_0_0_dX_dY[index4D]
										- 4.0 * ( int_dup_d_Phi_0_0_dX[index2Da] / ( N_00_00*N_00_00 ) )
										* (
											int_dup_d_Phi_o_o_dX[index2Db]
												 - ( N_00_oo / N_00_00 ) * int_dup_d_Phi_0_0_dX[index2Db]
											);
		theta_1_ot_XY[index4D] = d_Phi_o_t_dX_dY[index4D] / N_00_00
										- (N_00_ot / ( N_00_00*N_00_00 )) * d_Phi_0_0_dX_dY[index4D]
										- 4.0 * ( int_dup_d_Phi_0_0_dX[index2Da] / ( N_00_00*N_00_00 ) )
										* (
											int_dup_d_Phi_o_t_dX[index2Db]
												 - ( N_00_ot / N_00_00 ) * int_dup_d_Phi_0_0_dX[index2Db]
											);
		theta_1_tt_XY[index4D] = d_Phi_t_t_dX_dY[index4D] / N_00_00
										- (N_00_tt / ( N_00_00*N_00_00 )) * d_Phi_0_0_dX_dY[index4D]
										- 4.0 * ( int_dup_d_Phi_0_0_dX[index2Da] / ( N_00_00*N_00_00 ) )
										* (
											int_dup_d_Phi_t_t_dX[index2Db]
												 - ( N_00_tt / N_00_00 ) * int_dup_d_Phi_0_0_dX[index2Db]
											);
	}

	return;
}


// End of file

#endif
