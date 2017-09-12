#ifndef PHI_FUNCTIONS_AND_DERIVATIVES_H
#define PHI_FUNCTIONS_AND_DERIVATIVES_H

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
#include "Psi_functions_and_derivatives.h"

////////////////////////////////////////////////////////////////////
// Phi functions which use above Psi functions
////////////////////////////////////////////////////////////////////

inline void set_Phi_ij_grids()
{
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		double u = u_pts[iu], up = u_pts[iup];
		double sh_xi = sqrt(u*u-1.0), ch_xi = u;
		double sh_xip = sqrt(up*up-1.0), ch_xip = up;

		double Psi_0_u = Psi_0[iu], Psi_0_up = Psi0[iup];
		double Psi_1_u = Psi_1[iu], Psi_1_up = Psi1[iup];
		double Psi_2_u = Psi_2[iu], Psi_2_up = Psi2[iup];

		Phi_0_0[iu][iup] = Psi_0_u * Psi_0_up;
		Phi_s_s[iu][iup] = 0.5*tauf*tauf*( ( Psi_0_u - Psi_2_u) * Psi_0_up * sh_xi * sh_xi
											 + ( Psi_0_up - Psi_2_up) * Psi_0_u * sh_xip * sh_xip );
		Phi_o_o[iu][iup] = 0.5*tauf*tauf*( Psi_2_u * Psi_0_up * sh_xi * sh_xi
											 + Psi_2_up * Psi_0_u * sh_xip * sh_xip
											 - 2.0 * Psi_1_u * Psi_1_up * sh_xi * sh_xip );
		Phi_o_t[iu][iup] = 0.5*tauf*tauf*( ch_xi - ch_xip ) * ( Psi_1_u * Psi_0_up * sh_xi
																	- Psi_1_up * Psi_0_u * sh_xip );
		Phi_t_t[iu][iup] = 0.5*tauf*tauf*Psi_0_u * Psi_0_up * ( ch_xi - ch_xip ) * ( ch_xi - ch_xip );
	}

	return;
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

		double dPsi_0_dX_u = dPsi0_dX[iX][iu], Psi_0_up = Psi0[iup];
		double dPsi_1_dX_u = dPsi1_dX[iX][iu], Psi_1_up = Psi1[iup];
		double dPsi_2_dX_u = dPsi2_dX[iX][iu], Psi_2_up = Psi2[iup];

		d_Phi_0_0_dX[iX][iu][iup] = dPsi_0_dX_u * Psi_0_up;
		d_Phi_s_s_dX[iX][iu][iup] = 0.5*tauf*tauf*( ( dPsi_0_dX_u - dPsi_2_dX_u) * Psi_0_up * sh_xi * sh_xi
											 + ( Psi_0_up - Psi_2_up) * dPsi_0_dX_u * sh_xip * sh_xip );
		d_Phi_o_o_dX[iX][iu][iup] = 0.5*tauf*tauf*( dPsi_2_dX_u * Psi_0_up * sh_xi * sh_xi
											 + Psi_2_up * dPsi_0_dX_u * sh_xip * sh_xip
											 - 2.0 * dPsi_1_dX_u * Psi_1_up * sh_xi * sh_xip );
		d_Phi_o_t_dX[iX][iu][iup] = 0.5*tauf*tauf*( ch_xi - ch_xip ) * ( dPsi_1_dX_u * Psi_0_up * sh_xi
																	- Psi_1_up * dPsi_0_dX_u * sh_xip );
		d_Phi_t_t_dX[iX][iu][iup] = 0.5*tauf*tauf*dPsi_0_dX_u * Psi_0_up * ( ch_xi - ch_xip ) * ( ch_xi - ch_xip );
	}

	return;
}

inline void set_d_Phi_ij_dX_dY_grids()
{
	for (int iX = 0; iX < 3; ++iX)
	for (int iY = 0; iY < 3; ++iY)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		double u = u_pts[iu], up = u_pts[iup];
		double sh_xi = sqrt(u*u-1.0), ch_xi = u;
		double sh_xip = sqrt(up*up-1.0), ch_xip = up;

		double dPsi_0_dX_dY_u = dPsi0_dX_dY[iX][iY][iu], Psi_0_up = Psi_0[iup];
		double dPsi_1_dX_dY_u = dPsi1_dX_dY[iX][iY][iu], Psi_1_up = Psi_1[iup];
		double dPsi_2_dX_dY_u = dPsi2_dX_dY[iX][iY][iu], Psi_2_up = Psi_2[iup];

		d_Phi_0_0_dX_dY[iX][iY][iu][iup] = dPsi_0_dX_dY_u * Psi_0_up;
		d_Phi_s_s_dX_dY[iX][iY][iu][iup] = 0.5*tauf*tauf*( ( dPsi_0_dX_dY_u - dPsi_2_dX_dY_u) * Psi_0_up * sh_xi * sh_xi
											 + ( Psi_0_up - Psi_2_up) * dPsi_0_dX_dY_u * sh_xip * sh_xip );
		d_Phi_o_o_dX_dY[iX][iY][iu][iup] = 0.5*tauf*tauf*( dPsi_2_dX_dY_u * Psi_0_up * sh_xi * sh_xi
											 + Psi_2_up * dPsi_0_dX_dY_u * sh_xip * sh_xip
											 - 2.0 * dPsi_1_dX_dY_u * Psi_1_up * sh_xi * sh_xip );
		d_Phi_o_t_dX_dY[iX][iY][iu][iup] = 0.5*tauf*tauf*( ch_xi - ch_xip ) * ( dPsi_1_dX_dY_u * Psi_0_up * sh_xi
																	- Psi_1_up * dPsi_0_dX_dY_u * sh_xip );
		d_Phi_t_t_dX_dY[iX][iY][iu][iup] = 0.5*tauf*tauf*dPsi_0_dX_dY_u * Psi_0_up * ( ch_xi - ch_xip ) * ( ch_xi - ch_xip );
	}

	return;
}

//notice that this function is different from the preceding one ( dY <--> dYp )
inline void set_d_Phi_ij_dX_dYp_grids()
{
	for (int iX = 0; iX < 3; ++iX)
	for (int iYp = 0; iYp < 3; ++iYp)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
	{
		double u = u_pts[iu], up = u_pts[iup];
		double sh_xi = sqrt(u*u-1.0), ch_xi = u;
		double sh_xip = sqrt(up*up-1.0), ch_xip = up;

		double dPsi_0_dX_u = dPsi0_dX[iX][iu], dPsi_0_dYp_up = dPsi0_dX[iYp][iup];
		double dPsi_1_dX_u = dPsi1_dX[iX][iu], dPsi_1_dYp_up = dPsi1_dX[iYp][iup];
		double dPsi_2_dX_u = dPsi2_dX[iX][iu], dPsi_2_dYp_up = dPsi2_dX[iYp][iup];

		d_Phi_0_0_dX_dYp[iX][iYp][iu][iup] = dPsi_0_dX_u * dPsi_0_dYp_up;
		d_Phi_s_s_dX_dYp[iX][iYp][iu][iup] = 0.5*tauf*tauf*( ( dPsi_0_dX_u - dPsi_2_dX_u) * dPsi_0_dYp_up * sh_xi * sh_xi
											 + ( dPsi_0_dYp_up - dPsi_2_dYp_up) * dPsi_0_dX_u * sh_xip * sh_xip );
		d_Phi_o_o_dX_dYp[iX][iYp][iu][iup] = 0.5*tauf*tauf*( dPsi_2_dX_u * dPsi_0_dYp_up * sh_xi * sh_xi
											 + dPsi_2_dYp_up * dPsi_0_dX_u * sh_xip * sh_xip
											 - 2.0 * dPsi_1_dX_u * dPsi_1_dYp_up * sh_xi * sh_xip );
		d_Phi_o_t_dX_dYp[iX][iYp][iu][iup] = 0.5*tauf*tauf*( ch_xi - ch_xip ) * ( dPsi_1_dX_u * dPsi_0_dYp_up * sh_xi
																	- dPsi_1_dYp_up * dPsi_0_dX_u * sh_xip );
		d_Phi_t_t_dX_dYp[iX][iYp][iu][iup] = 0.5*tauf*tauf*dPsi_0_dX_u * dPsi_0_dYp_up * ( ch_xi - ch_xip ) * ( ch_xi - ch_xip );
	}

	return;
}

inline void set_integrated_d_Phi_ij_dX_grids()
{
	for (int iX = 0; iX < 3; ++iX)
	for (int iu = 0; iu < n_u_pts; ++iu)
	{
		for (int iup = 0; iup < n_u_pts; ++iup)
		{
			double up_loc = u_pts[iup];
			int_dup_d_Phi_0_0_dX[iX][iu] += u_wts[iup] * d_Phi_0_0_dX[iX][iu][iup] / sqrt(up_loc*up_loc-1.0);
			int_dup_d_Phi_s_s_dX[iX][iu] += u_wts[iup] * d_Phi_s_s_dX[iX][iu][iup] / sqrt(up_loc*up_loc-1.0);
			int_dup_d_Phi_o_o_dX[iX][iu] += u_wts[iup] * d_Phi_o_o_dX[iX][iu][iup] / sqrt(up_loc*up_loc-1.0);
			int_dup_d_Phi_o_t_dX[iX][iu] += u_wts[iup] * d_Phi_o_t_dX[iX][iu][iup] / sqrt(up_loc*up_loc-1.0);
			int_dup_d_Phi_t_t_dX[iX][iu] += u_wts[iup] * d_Phi_t_t_dX[iX][iu][iup] / sqrt(up_loc*up_loc-1.0);
		}
	}

	return;
}

inline void set_integrated_d_Phi_ij_dX_dY_grids()
{
	for (int iX = 0; iX < 3; ++iX)
	for (int iY = 0; iY < 3; ++iY)
	for (int iu = 0; iu < n_u_pts; ++iu)
	{
		for (int iup = 0; iup < n_u_pts; ++iup)
		{
			double up_loc = u_pts[iup];
			int_dup_d_Phi_0_0_dX_dY[iX][iY][iu] += u_wts[iup] * d_Phi_0_0_dX_dY[iX][iY][iu][iup] / sqrt(up_loc*up_loc-1.0);
			int_dup_d_Phi_s_s_dX_dY[iX][iY][iu] += u_wts[iup] * d_Phi_s_s_dX_dY[iX][iY][iu][iup] / sqrt(up_loc*up_loc-1.0);
			int_dup_d_Phi_o_o_dX_dY[iX][iY][iu] += u_wts[iup] * d_Phi_o_o_dX_dY[iX][iY][iu][iup] / sqrt(up_loc*up_loc-1.0);
			int_dup_d_Phi_o_t_dX_dY[iX][iY][iu] += u_wts[iup] * d_Phi_o_t_dX_dY[iX][iY][iu][iup] / sqrt(up_loc*up_loc-1.0);
			int_dup_d_Phi_t_t_dX_dY[iX][iY][iu] += u_wts[iup] * d_Phi_t_t_dX_dY[iX][iY][iu][iup] / sqrt(up_loc*up_loc-1.0);
		}
	}

	return;
}

// End of file

#endif
