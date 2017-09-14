#ifndef DEFS2_H
#define DEFS2_H

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
#include "legendre.h"
#include "HBT_and_emission_density.h"

inline int indexer(int ik, int it1, int it2)
{
	return ( (2*n_tau_pts * ik+it1)*2*n_tau_pts+it2 );
}

inline void get_roots(double tau, double k, vector<complex<double> > & roots, complex<double> &zeta0)
{
	double T = interpolate1D(all_tau_pts, all_T_pts, tau, 0, false);
	double mu = interpolate1D(all_tau_pts, all_mu_pts, tau, 0, false);
	double n0 = n(T, mu);
	double w0 = w(T, mu);
	double vsig2 = vsigma2(T, mu);
	double vs2_0 = vs2(T, mu);
	double vn2_0 = vn2(T, mu);

	double b = -2.0*vsig2;
	double c = (k*k+0.25)*vsig2-(1.0-2.0*vsig2);
	double d = (k*k+0.25)*(vn2_0-vs2_0)*mu*n0 / w0;	
	zeta0 = d;

	complex<double> r1, r2, r3;

	solve_cubic_equation(b, c, d, r1, r2, r3);

	roots.push_back( r1 );
	roots.push_back( r2 );
	roots.push_back( r3 );

//cout << "nip it in the butt: " << b << "   " << c << "   " << d << "   " << r1 << "   " << r2 << "   " << r3 << endl;

	return;
}

inline void get_roots_at_tauf(double k, vector<complex<double> > & roots, complex<double> &zeta0)
{
	double n0 = n(Tf, muf);
	double w0 = w(Tf, muf);
	double vsig2 = vsigma2(Tf, muf);
	double vs2_0 = vs2(Tf, muf);
	double vn2_0 = vn2(Tf, muf);

	double b = -2.0*vsig2;
	double c = (k*k+0.25)*vsig2-(1.0-2.0*vsig2);
	double d = (k*k+0.25)*(vn2_0-vs2_0)*muf*n0 / w0;	
	zeta0 = d;

	complex<double> r1, r2, r3;

	solve_cubic_equation(b, c, d, r1, r2, r3);

	roots.push_back( r1 );
	roots.push_back( r2 );
	roots.push_back( r3 );

	return;
}

inline void set_G3_and_tauDtau_G3_matrices()
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < 2*n_tau_pts; ++it)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		vector<complex<double> > roots;
		complex<double> zeta0;
		double k = k_pts[ik];
		double tau = all_tau_pts[it];
		double taup = all_tau_pts[itp];

		get_roots(tau, k, roots, zeta0);	//assumes coefficients are evaluated at tau, not tau'

		complex<double> r1 = roots[0];
		complex<double> r2 = roots[1];
		complex<double> r3 = roots[2];
		complex<double> a1 = (r2*r3 + zeta0) * pow(tau/taup, r1) / ((r1-r2)*(r1-r3));
		complex<double> a2 = (r3*r1 + zeta0) * pow(tau/taup, r2) / ((r2-r3)*(r2-r1));
		complex<double> a3 = (r1*r2 + zeta0) * pow(tau/taup, r3) / ((r3-r1)*(r3-r2));

		G3_tau_taup[ik][it][itp] = a1 + a2 + a3;
		tauDtau_G3_tau_taup[ik][it][itp] = r1*a1 + r2*a2 + r3*a3;
	}
	return;
}

inline void set_G3_and_tauDtau_G3_matrices_at_tauf()
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		vector<complex<double> > roots;
		complex<double> zeta0;
		double k = k_pts[ik];
		double taup = all_tau_pts[itp];

		get_roots_at_tauf(k, roots, zeta0);	//assumes coefficients are evaluated at tau, not tau'

		complex<double> r1 = roots[0];
		complex<double> r2 = roots[1];
		complex<double> r3 = roots[2];
		complex<double> a1 = (r2*r3 + zeta0) * pow(tauf/taup, r1) / ((r1-r2)*(r1-r3));
		complex<double> a2 = (r3*r1 + zeta0) * pow(tauf/taup, r2) / ((r2-r3)*(r2-r1));
		complex<double> a3 = (r1*r2 + zeta0) * pow(tauf/taup, r3) / ((r3-r1)*(r3-r2));

		G3_tauf_taup[ik][itp] = a1 + a2 + a3;

		tauDtau_G3_tauf_taup[ik][itp] = r1*a1 + r2*a2 + r3*a3;
	}
	return;
}

inline void set_A1_pts()
{
	for (int it = 0; it < 2*n_tau_pts; ++it)
	{
		double T_loc = all_T_pts[it];
		double mu_loc = all_mu_pts[it];
		A1_pts[it] = (mu_loc/T_loc) + sPERn;
	}
	return;
}

inline void set_A2_pts()
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < 2*n_tau_pts; ++it)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double T_loc = all_T_pts[it];
		double mu_loc = all_mu_pts[it];

		double mu_BY_T = mu_loc/T_loc;

		A2_pts[indexer(ik,it,itp)] = mu_BY_T * G3_tau_taup[ik][it][itp] + sPERn * tauDtau_G3_tau_taup[ik][it][itp];
//cout << "set_A2_pts(): " << ik << "   " << it << "   " << itp << "   " << mu_BY_T << "   " << G3_tau_taup[ik][it][itp] << "   " << sPERn << "   " << tauDtau_G3_tau_taup[ik][it][itp] << endl;
	}
	return;
}

inline void set_B_pts()
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
		B_pts[ik][itp] = sPERn * tauDtau_G3_tauf_taup[ik][itp];
	return;
}

inline void set_C_pts()
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
		C_pts[ik][itp] = G3_tauf_taup[ik][itp];
	return;
}

inline double F_11_11()
{
	//tau1==tau2==tauf
	double sum = 0.0;
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double taup_loc = all_tau_pts[itp];
		sum += all_tau_wts[itp] * transport_pts[itp] * A1_pts[itp] * A1_pts[itp] / pow(taup_loc, 7.0);
	}
	return (sum);
}

//NOT CHECKED YET
inline complex<double> F_12_11(int ik2)
{
	//tau1==tau2==tauf
	complex<double> sum = 0.0;
	for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
	for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		//sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it1p] * A1_pts[it1p] * A2_pts[ik2][it2p][it1p] / (tau2p*tau2p*pow(tau1p, 6.0) );
		sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it1p] * A1_pts[it1p] * A2_pts[indexer(ik2,it2p,it1p)] / (tau2p*tau2p*pow(tau1p, 6.0) );
	}
	return (sum);
}

//NOT CHECKED YET
inline complex<double> F_21_11(int ik1)
{
	//tau1==tau2==tauf
	complex<double> sum = 0.0;
	for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
	for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		//sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A1_pts[it2p] * A2_pts[ik1][it1p][it2p] / (tau1p*tau1p*pow(tau2p, 6.0) );
		sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A1_pts[it2p] * A2_pts[indexer(ik1,it1p,it2p)] / (tau1p*tau1p*pow(tau2p, 6.0) );
	}
	return (sum);
}

//NOT CHECKED YET
//inline complex<double> F_22_11(int ik1, int ik2)
inline void set_F_22_11_pts()
{
	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		//tau1==tau2==tauf
		complex<double> sum = 0.0;
		for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
		for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
		{
			double tau1p = all_tau_pts[it1p];
			double tau2p = all_tau_pts[it2p];
			double min_t1_t2 = min(tau1p, tau2p);
			double hw = 0.5 * (min_t1_t2 - tau0);
			double cen = 0.5 * (min_t1_t2 + tau0);
			//complex<double> tmpA2_k1 = A2_pts[ik1][it1p][it2p];
			//complex<double> tmpA2_k2 = A2_pts[ik2][it1p][it2p];
			complex<double> tmpA2_k1 = A2_pts[indexer(ik1,it1p,it2p)];
			complex<double> tmpA2_k2 = A2_pts[indexer(ik1,it1p,it2p)];
			for (int itpp = 0; itpp < n_x_pts; ++itpp)
			{
				double taupp = cen + hw * x_pts[itpp];
				sum += hw * x_wts[itpp] * all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[itpp]
						* tmpA2_k1 * tmpA2_k2 / (tau1p*tau1p*tau2p*tau2p*pow(taupp, 5.0));
			}
		}
		F_22_11_pts[ik1][ik2] = sum;
	}
	return;
}

//NOT CHECKED YET
inline complex<double> F_1_12(int ik2)
{
	//tau1==tau2==tauf
	complex<double> sum = 0.0;
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double taup = all_tau_pts[itp];
		sum += all_tau_wts[itp] * transport_pts[itp] * A1_pts[itp] * B_pts[ik2][itp] / pow(taup, 6.0);
	}
	return (sum);
}

//NOT CHECKED YET
inline complex<double> F_2_12(int ik1, int ik2)
{
	//tau1==tau2==tauf
	complex<double> sum = 0.0;
	for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
	for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		//sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A2_pts[ik1][it1p][it2p] * B_pts[ik2][it2p] / (tau1p*tau1p*pow(tau2p, 5.0) );
		sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A2_pts[indexer(ik1,it1p,it2p)] * B_pts[ik2][it2p] / (tau1p*tau1p*pow(tau2p, 5.0) );
	}
	return (sum);
}

//NOT CHECKED YET
inline complex<double> F_1_13(int ik2)
{
	//tau1==tau2==tauf
	complex<double> sum = 0.0;
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double taup = all_tau_pts[itp];
		sum += all_tau_wts[itp] * transport_pts[itp] * A1_pts[itp] * C_pts[ik2][itp] / pow(taup, 6.0);
	}
	return (sum);
}

//NOT CHECKED YET
inline complex<double> F_2_13(int ik1, int ik2)
{
	//tau1==tau2==tauf
	complex<double> sum = 0.0;
	for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
	for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		//sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A2_pts[ik1][it1p][it2p] * C_pts[ik2][it2p] / (tau1p*tau1p*pow(tau2p, 5.0) );
		sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A2_pts[indexer(ik1,it1p,it2p)] * C_pts[ik2][it2p] / (tau1p*tau1p*pow(tau2p, 5.0) );
	}
	return (sum);
}

//////////////////////////////////////////////////////////////////
// all two-point functions <del_i del_j> require 2D MF-inversion,
// so compute these in functions below (a few 2pt functions may
// require further integrations)
//////////////////////////////////////////////////////////////////

inline void set_T_11()
{
	double local_F_11_11 = F_11_11();
	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		F_12_11_pts[ik] = F_12_11(ik);
		F_21_11_pts[ik] = F_21_11(ik);
	}

	set_F_22_11_pts();

	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		Tarray[0][0][ik1][ik2] = tauf*tauf*legendre_integral_array[ik1][ik2]
							* (local_F_11_11 + F_12_11_pts[ik2] + F_21_11_pts[ik1] + F_22_11_pts[ik1][ik2]) / (2.0*M_PI);
//cout << ik1 << "   " << ik2 << "   " << legendre_integral_array[ik1][ik2] << "   " << local_F_11_11 << "   " << F_12_11_pts[ik2] << "   " << F_21_11_pts[ik1] << "   " << F_22_11_pts[ik1][ik2] << endl;
	}
	return;
}

inline void set_T_22()
{
	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		complex<double> tmp = 0.0;
		for (int itp = 0; itp < 2*n_tau_pts; ++itp)
		{
			double taup = all_tau_pts[itp];
			tmp += tauf*tauf*legendre_integral_array[ik1][ik2]
							* transport_pts[itp] * B_pts[ik1][itp] * B_pts[ik2][itp] / (2.0*M_PI*pow(taup,5.0));	//only one time-index since first one evaluated at tau_f
		}
		Tarray[1][1][ik1][ik2] = tmp;
	}
	return;
}

inline void set_T_33()
{
	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		complex<double> tmp = 0.0;
		for (int itp = 0; itp < 2*n_tau_pts; ++itp)
		{
			double taup = all_tau_pts[itp];
			tmp += tauf*tauf*legendre_integral_array[ik1][ik2]
							* transport_pts[itp] * C_pts[ik1][itp] * C_pts[ik2][itp] / (2.0*M_PI*taup*taup);	//only one time-index since first one evaluated at tau_f
		}
		Tarray[2][2][ik1][ik2] = tmp;
	}
	return;
}

inline void set_T_12()
{
	for (int ik = 0; ik < n_k_pts; ++ik)
		F_1_12_pts[ik] = F_1_12(ik);

	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		F_2_12_pts[ik1][ik2] = F_2_12(ik1,ik2);
		Tarray[0][1][ik1][ik2] = tauf*tauf*legendre_integral_array[ik1][ik2]
							* (F_1_12_pts[ik2] + F_2_12_pts[ik1][ik2]) / (2.0*M_PI);
	}
	return;
}

inline void set_T_13()
{
	for (int ik = 0; ik < n_k_pts; ++ik)
		F_1_13_pts[ik] = F_1_13(ik);

	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		F_2_13_pts[ik1][ik2] = F_2_13(ik1,ik2);
		Tarray[0][2][ik1][ik2] = tauf*tauf*legendre_integral_array[ik1][ik2]
							* (F_1_13_pts[ik2] + F_2_13_pts[ik1][ik2]) / (2.0*M_PI);
	}
	return;
}

inline void set_T_23()
{
	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		complex<double> tmp = 0.0;
		for (int itp = 0; itp < 2*n_tau_pts; ++itp)
		{
			double taup = all_tau_pts[itp];
			tmp += tauf*tauf*legendre_integral_array[ik1][ik2]
							* transport_pts[itp] * B_pts[ik1][itp] * C_pts[ik2][itp] / (2.0*M_PI*pow(taup,5.0));	//only one time-index since first one evaluated at tau_f
		}
		Tarray[1][2][ik1][ik2] = tmp;
	}
	return;
}


inline void set_T_XY()
{
	//set diagonals
	set_T_11();
	set_T_22();
	set_T_33();

	//set off-diagonals
	set_T_12();
	set_T_13();
	set_T_23();

	//use symmetry to set the rest
	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		Tarray[1][0][ik1][ik2] = Tarray[0][1][ik1][ik2];
		Tarray[2][0][ik1][ik2] = Tarray[0][2][ik1][ik2];
		Tarray[2][1][ik1][ik2] = Tarray[1][2][ik1][ik2];
	}

	return;
}

inline void set_everything_else()
{
	//transport related functions here...
	set_G3_and_tauDtau_G3_matrices();
	set_G3_and_tauDtau_G3_matrices_at_tauf();

	set_A1_pts();
	set_A2_pts();
	set_B_pts();
	set_C_pts();

	//HBT related functions here...
	//set earliest stuff first
	set_PsiA();
	set_PsiA_first_derivatives_vector();
	set_PsiA_second_derivatives_array();
	set_PsiBk();
	set_PsiBk_first_derivatives_vector();
	set_PsiBk_second_derivatives_array();

	//next layer of dependence
	set_Psi_k();
	set_dPsik_dX();
	set_dPsik_dX_dY();

	//next layer
	set_Phi_ij_grids();
	set_d_Phi_ij_dX_grids();
	set_d_Phi_ij_dX_dY_grids();
	set_d_Phi_ij_dX_dYp_grids();
	set_integrated_d_Phi_ij_dX_grids();
	set_integrated_d_Phi_ij_dX_dY_grids();
	set_N_00_ij();

	//final layer
	//set_theta_0_ij_XY();
	//set_theta_1_ij_XY();

	return;
}

inline void set_mean_delta_R2ij(int chosen_trajectory, int particle_to_study)
{
	initialize_all(chosen_trajectory, particle_to_study);

	complex<double> result_ss = 0.0, result_oo = 0.0, result_ot = 0.0, result_tt = 0.0;

	//legendre stuff here
	set_Q_X_k(QXk, k_pts, u_pts);
	compute_legendre_integral(legendre_integral_array, k_pts);

	//set other stuff here
	set_everything_else();

	set_T_XY();

	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	for (int iX = 0; iX < 3; ++iX)
	for (int iY = 0; iY < 3; ++iY)
	{
		double k1 = k_pts[ik1];
		double k2 = k_pts[ik2];
		for (int iu = 0; iu < n_u_pts; ++iu)
		{
			result_ss += k_wts[ik1]*k_wts[ik2]*k1*k2*tanh(M_PI*k1)*tanh(M_PI*k2)
						*Tarray[iX][iY][ik1][ik2]*QXk[iX][ik1][iu]*QXk[iY][ik2][iu]
						*theta_0_ss_XY[iX][iY][iu];
			result_oo += k_wts[ik1]*k_wts[ik2]*k1*k2*tanh(M_PI*k1)*tanh(M_PI*k2)
						*Tarray[iX][iY][ik1][ik2]*QXk[iX][ik1][iu]*QXk[iY][ik2][iu]
						*theta_0_oo_XY[iX][iY][iu];
			result_ot += k_wts[ik1]*k_wts[ik2]*k1*k2*tanh(M_PI*k1)*tanh(M_PI*k2)
						*Tarray[iX][iY][ik1][ik2]*QXk[iX][ik1][iu]*QXk[iY][ik2][iu]
						*theta_0_ot_XY[iX][iY][iu];
			result_tt += k_wts[ik1]*k_wts[ik2]*k1*k2*tanh(M_PI*k1)*tanh(M_PI*k2)
						*Tarray[iX][iY][ik1][ik2]*QXk[iX][ik1][iu]*QXk[iY][ik2][iu]
						*theta_0_tt_XY[iX][iY][iu];
//if (abs(Tarray[iX][iY][ik1][ik2]) > 1.e-10)
//	cout << "line 1: " << iX << "   " << iY<< "   " << ik1 << "   " << ik2 << "   " << Tarray[iX][iY][ik1][ik2] << "   " << iu << "   " << theta0XY[iX][iY][iu] << endl;
		}
		for (int iu1 = 0; iu1 < n_u_pts; ++iu1)
		for (int iu2 = 0; iu2 < n_u_pts; ++iu2)
		{
			result_ss += k_wts[ik1]*k_wts[ik2]*k1*k2*tanh(M_PI*k1)*tanh(M_PI*k2)
						*Tarray[iX][iY][ik1][ik2]*QXk[iX][ik1][iu1]*QXk[iY][ik2][iu2]
						*theta_1_ss_XY[iX][iY][iu1][iu2];
			result_oo += k_wts[ik1]*k_wts[ik2]*k1*k2*tanh(M_PI*k1)*tanh(M_PI*k2)
						*Tarray[iX][iY][ik1][ik2]*QXk[iX][ik1][iu1]*QXk[iY][ik2][iu2]
						*theta_1_oo_XY[iX][iY][iu1][iu2];
			result_ot += k_wts[ik1]*k_wts[ik2]*k1*k2*tanh(M_PI*k1)*tanh(M_PI*k2)
						*Tarray[iX][iY][ik1][ik2]*QXk[iX][ik1][iu1]*QXk[iY][ik2][iu2]
						*theta_1_ot_XY[iX][iY][iu1][iu2];
			result_tt += k_wts[ik1]*k_wts[ik2]*k1*k2*tanh(M_PI*k1)*tanh(M_PI*k2)
						*Tarray[iX][iY][ik1][ik2]*QXk[iX][ik1][iu1]*QXk[iY][ik2][iu2]
						*theta_1_tt_XY[iX][iY][iu1][iu2];
//cout << "line 2: " << iu1 << "   " << iu2 << "   " << theta1XY[iX][iY][iu1][iu2] << endl;
		}
	}

	return;
}


//End of file

#endif
