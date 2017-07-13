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

inline void set_G3_and_tauDtau_G3_matrices()
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < 2*n_tau_pts; ++it)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		vector<complex<double> > roots;
		vector<complex<double> > zeta0;
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

		tauDtau_G3_tau_taup[ik][it][itp] = r1*a1 + r2*a2 + f3*a3;
	}
	return;
}

inline double set_A1()
{
	for (int it = 0; it < 2*n_tau_pts; ++it)
	{
		double T_loc = all_T_pts[it];
		double mu_loc = all_mu_pts[it];
		A1_pts[it] = (mu_loc/T_loc) + sPERn];
	}
	return;
}

inline double A2(int ik, int it, int itp)
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < 2*n_tau_pts; ++it)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double T_loc = all_T_pts[it];
		double mu_loc = all_mu_pts[it];
		//double T_loc = interpolate1D(all_tau_pts, all_T_pts, tau_loc, 2*n_tau_pts, 0, false);
		//double mu_loc = interpolate1D(all_tau_pts, all_mu_pts, tau_loc, 2*n_tau_pts, 0, false);

		double mu_BY_T = mu_loc/T_loc;

		A2_pts[ik][it][itp] = mu_BY_T * G3_tau_taup[ik][it][itp] + sPERn * tauDtau_G3_tau_taup[ik][it][itp] );
	}
	return;
}

inline double B(int ik, int it, int itp)
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < 2*n_tau_pts; ++it)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
		B_pts[ik][it][itp] = sPERn * tauDtau_G3_tau_taup[ik][it][itp];
	return;
}

inline double C(int ik, int it, int itp)
{
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < 2*n_tau_pts; ++it)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
		C_pts[ik][it][itp] = G3_tau_taup[ik][it][itp];
	return;
}}

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
inline double F_12_11(int ik2)
{
	//tau1==tau2==tauf
	double sum = 0.0;
	for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
	for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it1p] * A1_pts[it1p] * A2_pts[ik2][it2p][it1p] / (tau2p*tau2p*pow(tau1p, 6.0) );
	}
	return (sum);
}

//NOT CHECKED YET
inline double F_21_11(int ik1)
{
	//tau1==tau2==tauf
	double sum = 0.0;
	for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
	for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A1_pts[it2p] * A2_pts[ik2][it1p][it2p] / (tau1p*tau1p*pow(tau2p, 6.0) );
	}
	return (sum);
}

//NOT CHECKED YET
inline double F_22_11(int ik1, int ik2)
{
	//tau1==tau2==tauf
	double sum = 0.0;
	for (int it = 0; it < 2*n_tau_pts; ++it)
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		double min_t1_t2 = min(tau1p, tau2p);
		double hw = 0.5 * (min_t1_t2 - tau0);
		double cen = 0.5 * (min_t1_t2 + tau0);
		for (int itpp = 0; itpp < n_x_pts; ++itpp)
		{
			double taupp = cen + hw * x_pts[itpp];
			sum += hw * x_wts[itpp] * all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[itpp]
					* A2_pts[ik2][it1p][it2p] * A2_pts[ik2][it1p][it2p] / (tau1p*tau1p*tau2p*tau2p*pow(taupp, 5.0));
		}
	}
	return (sum);
}

//NOT CHECKED YET
inline double F_1_12(int ik2)
{
	//tau1==tau2==tauf
	double sum = 0.0;
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double taup = all_tau_pts[itp];
		sum += all_tau_wts[itp] * transport_pts[itp] * A1_pts[itp] * B_pts[ik2][itf][itp] / pow(taup, 6.0);
	}
	return (sum);
}

//NOT CHECKED YET
inline double F_2_12(int ik1, int ik2)
{
	//tau1==tau2==tauf
	double sum = 0.0;
	for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
	for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A2_pts[ik1][it1p][it2p] * B_pts[ik2][itf][it2p] / (tau1p*tau1p*pow(tau2p, 5.0) );
	}
	return (sum);
}

//NOT CHECKED YET
inline double F_1_13(int ik2)
{
	//tau1==tau2==tauf
	double sum = 0.0;
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double taup = all_tau_pts[itp];
		sum += all_tau_wts[itp] * transport_pts[itp] * A1_pts[itp] * C_pts[ik2][itf][itp] / pow(taup, 6.0);
	}
	return (sum);
}

//NOT CHECKED YET
inline double F_2_13(int ik1, int ik2)
{
	//tau1==tau2==tauf
	double sum = 0.0;
	for (int it1p = 0; it1p < 2*n_tau_pts; ++it1p)
	for (int it2p = 0; it2p < 2*n_tau_pts; ++it2p)
	{
		double tau1p = all_tau_pts[it1p];
		double tau2p = all_tau_pts[it2p];
		sum += all_tau_wts[it1p] * all_tau_wts[it2p] * transport_pts[it2p] * A2_pts[ik1][it1p][it2p] * C_pts[ik2][itf][it2p] / (tau1p*tau1p*pow(tau2p, 5.0) );
	}
	return (sum);
}

//////////////////////////////////////////////////////////////////
// all two-point functions <del_i del_j> require 2D MF-inversion,
// so compute these in functions below (a few 2pt functions may
// require further integrations)
//////////////////////////////////////////////////////////////////

inline void set_MFintegrand_delta_1_delta_1(double ** array)
{
	double local_F_11_11 = F_11_11();
	double local_F_12_11[n_k_pts], local_F_21_11[n_k_pts];
	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		local_F_12_11[ik] = F_12_11(ik);
		local_F_21_11[ik] = F_21_11(ik);
	}

	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		array[ik1][ik2] = tauf*tauf*legendre_integral_array[ik1][ik2]
							* (local_F_11_11 + local_F_12_11[ik2]+ local_F_12_11[ik1] + F_22_11(ik1, ik2));
	}
	return;
}

inline void set_MFintegrand_delta_1_delta_2(double ** array)
{

}


inline void set_MFintegrand_delta_1_delta_3(double ** array)
{

}


inline void set_MFintegrand_delta_1_delta_1(double ** array)
{

}


inline void set_MFintegrand_delta_2_delta_2(double ** array)
{

}


inline void set_MFintegrand_delta_2_delta_3(double ** array)
{

}


inline void set_MFintegrand_delta_3_delta_1(double ** array)
{

}


inline void set_MFintegrand_delta_3_delta_2(double ** array)
{

}


inline void set_MFintegrand_delta_3_delta_3(double ** array)
{

}





//End of file

#endif
