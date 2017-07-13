#ifndef DEFS_H
#define DEFS_H

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


int integration_mode = 1;
bool include_baryon_chemical_potential_fluctations = true;
extern const double hbarC;
extern const double mu_pion;
extern const int n_xi_pts;
extern const int n_k_pts;
extern const int n_tau_pts;



inline double norm_int(double x)
{
	double cx = cosh(x);
	return (incompleteGamma3(mByT * cx) / (cx*cx));
}

inline double Fs(double x)
{
	double cx = cosh(x);
	
	double c1 = s_at_mu_part * chi_tilde_mu_mu;
	double c2 = -s_at_mu_part * double(include_baryon_chemical_potential_fluctations) * (chi_tilde_T_mu + chi_tilde_mu_mu * mu_part / Tf);

	return ( (c1 * incompleteGamma4(mByT * cx) + c2 * incompleteGamma3(mByT * cx) ) / (cx*cx) );
}

inline double Fomega(double x)
{
	double cx = cosh(x);
	
	return (Tf * tanh(x)*incompleteGamma4(mByT * cx) / (cx*cx));
}

inline double Fn(double x)
{
	double cx = cosh(x);
	
	double c1 = -s_at_mu_part * chi_tilde_T_mu;
	double c2 = s_at_mu_part * double(include_baryon_chemical_potential_fluctations) * (chi_tilde_T_T + chi_tilde_T_mu * mu_part / Tf);

	return ( (c1 * incompleteGamma4(mByT * cx) + c2 * incompleteGamma3(mByT * cx) ) / (cx*cx) );
}

inline complex<double> Ftilde_s(double k)
{
	return (integrate_1D_FT(Fs, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k));
}

inline complex<double> Ftilde_omega(double k)
{
	return (integrate_1D_FT(Fomega, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k));
}

inline complex<double> Ftilde_n(double k)
{
	return (integrate_1D_FT(Fn, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k));
}

inline complex<double> Gtilde_s(double k, double t_p)
{
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	complex<double> b = beta(k);
	double t_by_tp = tauf / t_p;
	double f0 = vs2_at_tauf;
	double pref = (vsigma2_at_tauf-vs2_at_tauf)*pow(t_by_tp, -a_at_tauf);
	complex<double> f1 = (a_at_tauf/b)*sinh(b*log(t_by_tp)) + cosh(b*log(t_by_tp));

	return (
		(-i * k *mu_loc / (T_loc * vsigma2_at_tauf) ) * (f0 + pref*f1)	//dimensionless
	);
}

inline complex<double> Gtilde_omega(double k, double t_p)
{
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	double t_by_tp = tauf / t_p;
	complex<double> b = beta(k);
	
	return (
		(k*k*sPERn/b) * ( vsigma2_at_tauf-vn2_at_tauf ) * pow(t_by_tp, -a_at_tauf) * sinh(b*log(t_by_tp))	//dimensionless
	);
}

inline complex<double> Gtilde_n(double k, double t_p)
{
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	complex<double> b = beta(k);
	double t_by_tp = tauf / t_p;
	double f0 = vn2_at_tauf;
	double pref = (vsigma2_at_tauf-vn2_at_tauf)*pow(t_by_tp, -a_at_tauf);
	complex<double> f1 = (a_at_tauf/b)*sinh(b*log(t_by_tp)) + cosh(b*log(t_by_tp));

	return (
		(i * k / vsigma2_at_tauf ) * (f0 + pref*f1)		//dimensionless
	);
}

inline complex<double> tau_integration(complex<double> (*Gtilde_X)(double, double), complex<double> (*Gtilde_Y)(double, double), double k)
{
	complex<double> locsum(0,0);
	if (integration_mode == 1)
	{
		//lower section
		for (int it = 0; it < n_tau_pts; ++it)
		{
			double t_loc = tau_pts_lower[it];
			double T_loc = T_pts_lower[it];
			double mu_loc = mu_pts_lower[it];
			double arg = n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) );
			locsum += tau_wts_lower[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(arg, 2.0)
					* (*Gtilde_X)(k, t_loc) * (*Gtilde_Y)(-k, t_loc);
		}
		//upper section
		for (int it = 0; it < n_tau_pts; ++it)
		{
			double t_loc = tau_pts_upper[it];
			double T_loc = T_pts_upper[it];
			double mu_loc = mu_pts_upper[it];
			double arg = n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) );
			locsum += tau_wts_upper[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(arg, 2.0)
					* (*Gtilde_X)(k, t_loc) * (*Gtilde_Y)(-k, t_loc);
		}
	}
	else
	{
		//both sections together
		for (int it = 0; it < n_tau_pts; ++it)
		{
			double t_loc = tau_pts[it];
			double T_loc = T_pts[it];
			double mu_loc = mu_pts[it];
			double arg = n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) );
			locsum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(arg, 2.0)
					* (*Gtilde_X)(k, t_loc) * (*Gtilde_Y)(-k, t_loc);
		}
	}

	return (locsum);
}

inline complex<double> Ctilde_s_s(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_s, Gtilde_s, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_s_omega(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_s, Gtilde_omega, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_s_n(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_s, Gtilde_n, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_omega_s(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_omega, Gtilde_s, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_omega_omega(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_omega, Gtilde_omega, k);
	return ( sum );	//fm^2
}

inline complex<double> Ctilde_omega_n(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_omega, Gtilde_n, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_n_s(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_n, Gtilde_s, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_n_omega(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_n, Gtilde_omega, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_n_n(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_n, Gtilde_n, k);

	return ( sum );	//fm^2
}

// End of file

#endif
