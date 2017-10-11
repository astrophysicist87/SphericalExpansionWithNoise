#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <algorithm>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "lib.h"
#include "defs1.h"
#include "defs2.h"
#include "HBT_and_emission_density.h"
#include "gauss_quadrature.h"

bool do_1p_calc;
bool do_HBT_calc;
bool scale_out_y_dependence = false;
const int particle_to_study = 2;	//1 is pion, 2 is proton

const double hbarC = 197.33;
const double k_infinity = 20.0;
const double u_infinity = 100.0;
const int n_Dy = 1;

const double mu_pion = 1.e-3 / hbarC;
const int n_u_pts = 100;
const int n_k_pts = 100;
const int n_x_pts = 51;
const int n_tau_pts = 201;

double muis[3] = {420.0, 620.0, 820.0};

double mT = 1000.0 * 0.26 / hbarC;
double pT = 0.0, betaT = 0.0;

int main(int argc, char *argv[])
{
	//set parameters corresponding to all different possible trajectories
	Ti = 250.0 / hbarC;		//initial trajectory temperature
	current_itau = 0;

	int chosen_trajectory = atoi(argv[1]);
	chosen_trajectory--;		//use for array indexing throughout

	pT = sqrt(mT*mT - (139.57*139.57/hbarC/hbarC));
	betaT = pT / mT;

	/////////////////////////////////////////////////
	//just for now...
	//initialize_all(chosen_trajectory, particle_to_study);
	//set_everything_else();

	//cout << setw(20) << setprecision(16) << tauf << "   " << pT << "   " << mT << "   " << Tf << "   " << muf << endl
	//		<< sf << "   " << chi_tilde_mu_mu << "   " << chi_tilde_T_mu << "   " << chi_tilde_T_T << endl;
	//if (1) exit(0);
	///////////////////

	/////////////////////////////////////////////////////////////////
	//Check Psi functions and their derivatives in this section    //
	/////////////////////////////////////////////////////////////////

	//success
	/*for (int iu = 0; iu < n_u_pts; ++iu)
		cout << setw(20) << setprecision(16) << u_pts[iu] << "   "
				<< PsiA[iu] << "   " << PsiB0[iu] << "   " << PsiB1[iu] << "   " << PsiB2[iu] 
				<< "   " << Psi0[iu] << "   " << Psi1[iu] << "   " << Psi2[iu] << endl;
	if (1) exit (0);*/

	//success
	/*for (int iX = 0; iX < 3; ++iX)
	for (int iu = 0; iu < n_u_pts; ++iu)
		cout << setw(20) << setprecision(16) << iX << "   " << u_pts[iu] << "   "
				<< dPsiA_dX[iX][iu] << "   " << dPsiB0_dX[iX][iu] << "   " << dPsiB1_dX[iX][iu] << "   " << dPsiB2_dX[iX][iu] << "   "
				<< dPsi0_dX[iX][iu] << "   " << dPsi1_dX[iX][iu] << "   " << dPsi2_dX[iX][iu] << endl;
	if (1) exit (0);*/

	//success
	/*for (int iX = 0; iX < 3; ++iX)
	for (int iY = 0; iY < 3; ++iY)
	for (int iu = 0; iu < n_u_pts; ++iu)
		cout << setw(20) << setprecision(16) << iX << "   " << iY << "   " << u_pts[iu] << "   "
				<< dPsiA_dX_dY[iX][iY][iu] << "   " << dPsiB0_dX_dY[iX][iY][iu] << "   " << dPsiB1_dX_dY[iX][iY][iu] << "   " << dPsiB2_dX_dY[iX][iY][iu] << endl;
	if (1) exit (0);*/

	/////////////////////////////////////////////////////////////////
	//Check Phi functions and their derivatives in this section    //
	/////////////////////////////////////////////////////////////////

	//success
	/*for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
		cout << setw(20) << setprecision(16) << u_pts[iu] << "   " << u_pts[iup] << "   "
			<< Phi_0_0[iu][iup] << "   " << Phi_s_s[iu][iup] << "   " << Phi_o_o[iu][iup] << "   "
			<< Phi_o_t[iu][iup] << "   " << Phi_t_t[iu][iup] << endl;
	if (1) exit (0);*/

	//success
	/*for (int iX = 0; iX < 3; ++iX)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
		cout << setw(20) << setprecision(16) << iX << "   " << u_pts[iu] << "   " << u_pts[iup] << "   "
			<< d_Phi_0_0_dX[iX][iu][iup] << "   " << d_Phi_s_s_dX[iX][iu][iup] << "   " << d_Phi_o_o_dX[iX][iu][iup] << "   "
			<< d_Phi_o_t_dX[iX][iu][iup] << "   " << d_Phi_t_t_dX[iX][iu][iup] << endl;
	if (1) exit (0);*/

	//success
	/*for (int iX = 0; iX < 3; ++iX)
	for (int iY = 0; iY < 3; ++iY)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
		cout << setw(20) << setprecision(16) << iX << "   " << iY << "   " << u_pts[iu] << "   " << u_pts[iup] << "   "
			<< d_Phi_0_0_dX_dY[iX][iY][iu][iup] << "   " << d_Phi_s_s_dX_dY[iX][iY][iu][iup] << "   " << d_Phi_o_o_dX_dY[iX][iY][iu][iup] << "   "
			<< d_Phi_o_t_dX_dY[iX][iY][iu][iup] << "   " << d_Phi_t_t_dX_dY[iX][iY][iu][iup] << endl;
	if (1) exit (0);*/

	//success
	/*for (int iX = 0; iX < 3; ++iX)
	for (int iYp = 0; iYp < 3; ++iYp)
	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int iup = 0; iup < n_u_pts; ++iup)
		cout << setw(20) << setprecision(16) << iX << "   " << iYp << "   " << u_pts[iu] << "   " << u_pts[iup] << "   "
			<< d_Phi_0_0_dX_dYp[iX][iYp][iu][iup] << "   " << d_Phi_s_s_dX_dYp[iX][iYp][iu][iup] << "   " << d_Phi_o_o_dX_dYp[iX][iYp][iu][iup] << "   "
			<< d_Phi_o_t_dX_dYp[iX][iYp][iu][iup] << "   " << d_Phi_t_t_dX_dYp[iX][iYp][iu][iup] << endl;
	if (1) exit (0);*/

	//success
	/*for (int iX = 0; iX < 3; ++iX)
	for (int iu = 0; iu < n_u_pts; ++iu)
		cout << setw(20) << setprecision(16) << iX << "   " << u_pts[iu] << "   "
			<< int_dup_d_Phi_0_0_dX[iX][iu] << "   " << int_dup_d_Phi_s_s_dX[iX][iu] << "   " << int_dup_d_Phi_o_o_dX[iX][iu] << "   "
			<< int_dup_d_Phi_o_t_dX[iX][iu] << "   " << int_dup_d_Phi_t_t_dX[iX][iu] << endl;
	if (1) exit (0);*/

	//success
	/*for (int iX = 0; iX < 3; ++iX)
	for (int iY = 0; iY < 3; ++iY)
	for (int iu = 0; iu < n_u_pts; ++iu)
		cout << setw(20) << setprecision(16) << iX << "   " << iY << "   " << u_pts[iu] << "   "
			<< int_dup_d_Phi_0_0_dX_dY[iX][iY][iu] << "   " << int_dup_d_Phi_s_s_dX_dY[iX][iY][iu] << "   " << int_dup_d_Phi_o_o_dX_dY[iX][iY][iu] << "   "
			<< int_dup_d_Phi_o_t_dX_dY[iX][iY][iu] << "   " << int_dup_d_Phi_t_t_dX_dY[iX][iY][iu] << endl;
	if (1) exit (0);*/

	//do full calculation
	//complex<double> mean_delta_R2ij = get_mean_delta_R2ij(chosen_trajectory, particle_to_study);
	set_mean_delta_R2ij(chosen_trajectory, particle_to_study);

	//cout << "mean_delta_R2ij = " << mean_delta_R2ij << endl;

	return 0;
}
