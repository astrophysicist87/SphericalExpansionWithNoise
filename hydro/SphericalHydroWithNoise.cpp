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

#include "defs.h"
#include "gauss_quadrature.h"

bool do_1p_calc;
bool do_HBT_calc;
bool scale_out_y_dependence = false;
const int particle_to_study = 2;	//1 is pion, 2 is proton

const double hbarC = 197.33;
const double k_infinity = 10.0;
const double xi_infinity = 4.0;
const int n_Dy = 1;

const double mu_pion = 1.e-3 / hbarC;
const int n_xi_pts = 200;
const int n_k_pts = 200;
const int n_x_pts = 101;
const int n_tau_pts = 200;

int main(int argc, char *argv[])
{
	//set parameters corresponding to all different possible trajectories
	Ti = 250.0 / hbarC;		//initial trajectory temperature
	double muis[3] = {420.0, 620.0, 820.0};
	current_itau = 0;

	int chosen_trajectory = atoi(argv[1]);
	chosen_trajectory--;		//use for array indexing throughout
	do_1p_calc = (bool)atoi(argv[2]);
	do_HBT_calc = (bool)atoi(argv[3]);

	//performs all necessary initializations
	initialize_all(chosen_trajectory, particle_to_study);

	//
	compute_thermodynamic_twopoint_functions();

	if (do_1p_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		//cout << "norm = " << setprecision(15) << norm << endl;

		vector<complex<double> > Fts_vec, Fto_vec, Ftn_vec;
		vector<complex<double> > Ctss_vec, Ctso_vec, Ctsn_vec, Ctos_vec, Ctoo_vec, Cton_vec, Ctns_vec, Ctno_vec, Ctnn_vec;

		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];
			current_kwt = k_wts[ik];
			Fts_vec.push_back(Ftilde_s(k));
			Fto_vec.push_back(Ftilde_omega(k));
			Ftn_vec.push_back(Ftilde_n(k));
			Ctss_vec.push_back(Ctilde_s_s(k));
			Ctso_vec.push_back(Ctilde_s_omega(k));
			Ctsn_vec.push_back(Ctilde_s_n(k));
			Ctos_vec.push_back(Ctilde_omega_s(k));
			Ctoo_vec.push_back(Ctilde_omega_omega(k));
			Cton_vec.push_back(Ctilde_omega_n(k));
			Ctns_vec.push_back(Ctilde_n_s(k));
			Ctno_vec.push_back(Ctilde_n_omega(k));
			Ctnn_vec.push_back(Ctilde_n_n(k));
		}

		for (int iDy = 0; iDy < n_Dy; iDy++)
		{
			double Delta_y = (double)iDy * 0.1;
			current_DY = Delta_y;
			complex<double> sum(0,0);
			complex<double> sum00(0,0), sum01(0,0), sum02(0,0), sum10(0,0), sum11(0,0), sum12(0,0), sum20(0,0), sum21(0,0), sum22(0,0);
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				double k = k_pts[ik];
				current_kwt = k_wts[ik];

				complex<double> Fts = Fts_vec[ik];
				complex<double> Fto = Fto_vec[ik];
				complex<double> Ftn = Ftn_vec[ik];
				complex<double> Ctss = Ctss_vec[ik];
				complex<double> Ctso = Ctso_vec[ik];
				complex<double> Ctsn = Ctsn_vec[ik];
				complex<double> Ctos = Ctos_vec[ik];
				complex<double> Ctoo = Ctoo_vec[ik];
				complex<double> Cton = Cton_vec[ik];
				complex<double> Ctns = Ctns_vec[ik];
				complex<double> Ctno = Ctno_vec[ik];
				complex<double> Ctnn = Ctnn_vec[ik];

				sum += k_wts[ik] * exp(i * k * Delta_y)
						* ( Fts * conj(Fts) * Ctss
							+ Fts * conj(Fto) * Ctso
							+ Fts * conj(Ftn) * Ctsn
							+ Fto * conj(Fts) * Ctos
							+ Fto * conj(Fto) * Ctoo
							+ Fto * conj(Ftn) * Cton
							+ Ftn * conj(Fts) * Ctns
							+ Ftn * conj(Fto) * Ctno
							+ Ftn * conj(Ftn) * Ctnn );
				sum00 += k_wts[ik] * exp(i * k * Delta_y) * Fts * conj(Fts) * Ctss;
				sum01 += k_wts[ik] * exp(i * k * Delta_y) * Fts * conj(Fto) * Ctso;
				sum02 += k_wts[ik] * exp(i * k * Delta_y) * Fts * conj(Ftn) * Ctsn;
				sum10 += k_wts[ik] * exp(i * k * Delta_y) * Fto * conj(Fts) * Ctos;
				sum11 += k_wts[ik] * exp(i * k * Delta_y) * Fto * conj(Fto) * Ctoo;
				sum12 += k_wts[ik] * exp(i * k * Delta_y) * Fto * conj(Ftn) * Cton;
				sum20 += k_wts[ik] * exp(i * k * Delta_y) * Ftn * conj(Fts) * Ctns;
				sum21 += k_wts[ik] * exp(i * k * Delta_y) * Ftn * conj(Fto) * Ctno;
				sum22 += k_wts[ik] * exp(i * k * Delta_y) * Ftn * conj(Ftn) * Ctnn;
			}

			complex<double> result = (exp(mu_part/Tf)*ds*tauf*Tf / (2.0*M_PI*M_PI * norm)) * sum;
			cout << setprecision(15) << Delta_y << "   " << 1.0*result.real() << endl;
		}
	}

	if (do_HBT_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		alpha0 = get_alpha0();
		psi0 = get_psi0();

		vector<complex<double> > Fbts_vec, Fbto_vec, Fbtn_vec;
		vector<complex<double> > Ctss_vec, Ctso_vec, Ctsn_vec, Ctos_vec, Ctoo_vec, Cton_vec, Ctns_vec, Ctno_vec, Ctnn_vec;

		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];
			Fbts_vec.push_back(Fbtilde_s(k));
			Fbto_vec.push_back(Fbtilde_omega(k));
			Fbtn_vec.push_back(Fbtilde_n(k));
			Ctss_vec.push_back(Ctilde_s_s(k));
			Ctso_vec.push_back(Ctilde_s_omega(k));
			Ctsn_vec.push_back(Ctilde_s_n(k));
			Ctos_vec.push_back(Ctilde_omega_s(k));
			Ctoo_vec.push_back(Ctilde_omega_omega(k));
			Cton_vec.push_back(Ctilde_omega_n(k));
			Ctns_vec.push_back(Ctilde_n_s(k));
			Ctno_vec.push_back(Ctilde_n_omega(k));
			Ctnn_vec.push_back(Ctilde_n_n(k));
		}

		for (int iDy = 0; iDy < n_Dy; iDy++)
		{
			double Delta_y = (double)iDy * 0.05;
			//option #1
			double y1 = Delta_y;
			double y2 = 0.0;
			//option #2
			//double y1 = 0.5*Delta_y;
			//double y2 = -0.5*Delta_y;

			//scale out y-dependence, if desired
			double cy1 = cosh(y1);
			double cy2 = cosh(y2);
			double cDy = cosh(y1-y2);
			double scale_out_y_dep_factor = 1.0;
			if (scale_out_y_dependence)
				scale_out_y_dep_factor = cy1*cy1*cy2*cy2 / (cDy*cDy);

			complex<double> sum(0,0);
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				double k = k_pts[ik];

				complex<double> Fbts = Fbts_vec[ik];
				complex<double> Fbto = Fbto_vec[ik];
				complex<double> Fbtn = Fbtn_vec[ik];

				complex<double> Ctss = Ctss_vec[ik];
				complex<double> Ctso = Ctso_vec[ik];
				complex<double> Ctsn = Ctsn_vec[ik];
				complex<double> Ctos = Ctos_vec[ik];
				complex<double> Ctoo = Ctoo_vec[ik];
				complex<double> Cton = Cton_vec[ik];
				complex<double> Ctns = Ctns_vec[ik];
				complex<double> Ctno = Ctno_vec[ik];
				complex<double> Ctnn = Ctnn_vec[ik];

				sum += k_wts[ik] * exp(i * k * Delta_y)
						* ( Fbts * conj(Fbts) * Ctss
							+ Fbts * conj(Fbto) * Ctso
							+ Fbts * conj(Fbtn) * Ctsn
							+ Fbto * conj(Fbts) * Ctos
							+ Fbto * conj(Fbto) * Ctoo
							+ Fbto * conj(Fbtn) * Cton
							+ Fbtn * conj(Fbts) * Ctns
							+ Fbtn * conj(Fbto) * Ctno
							+ Fbtn * conj(Fbtn) * Ctnn );
			}

			sum *= pow(tauf, 4.0)/ (cy1*cy1*cy2*cy2);

			double A = 1.0; //fm^2
			double dNdy = exp(mu_part/Tf)*ds*tauf*A*Tf*Tf*Tf*norm / (4.0*M_PI*M_PI);
			double mean_R2l_vs_Dy = 0.5*tauf*tauf*psi0 / (cDy*cDy);
			double mean_R2l_vs_y1 = 0.5*tauf*tauf*psi0 / (cy1*cy1);
			double mean_R2l_vs_y2 = 0.5*tauf*tauf*psi0 / (cy2*cy2);
			complex<double> result = (exp(mu_part/Tf)*ds*tauf*Tf*Tf*Tf*norm / (2.0*M_PI*M_PI))*scale_out_y_dep_factor * sum / (mean_R2l_vs_Dy);
			complex<double> result2 = (exp(mu_part/Tf)*ds*tauf*Tf*Tf*Tf*norm / (2.0*M_PI*M_PI))*scale_out_y_dep_factor * sum / (mean_R2l_vs_y1*mean_R2l_vs_y2);
			cout << Delta_y << "   " << dNdy << "   " << mean_R2l_vs_Dy << "   " << mean_R2l_vs_y1 << "   " << mean_R2l_vs_y2 << "   " << result.real() << "   " << result2.real() << endl;
		}
	}

	return 0;
}
