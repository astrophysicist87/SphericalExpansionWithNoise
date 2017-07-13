#ifndef DEFS1_H
#define DEFS1_H

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

double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, m, sf, s_at_mu_part;
double mByT, alpha0, psi0;
double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

double exp_delta, exp_gamma, exp_nu;
double T0, mu0, Tc, Pc, nc, sc, wc, muc;
double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si, ni;
double a_at_tauf, vs2_at_tauf, vn2_at_tauf, vsigma2_at_tauf;

double current_kwt, current_DY;
int current_itau;
double mu_proton, mu_part;

double * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;
double * tau_pts_lower, * tau_wts_lower;
double * tau_pts_upper, * tau_wts_upper;
double * T_pts, * mu_pts;
double * T_pts_lower, * mu_pts_lower;
double * T_pts_upper, * mu_pts_upper;
double * all_T_pts, * all_mu_pts;
double * Delta_lambda_pts, * vn2_pts, * vs2_pts, * vsigma2_pts, * n_Tmu_pts, * s_Tmu_pts, * w_Tmu_pts;

//double G3_tau_taup[][][], tauDtau_G3_tau_taup[][][];
//double A1pts[], A2pts[][][], Bpts[][][], Cpts[][][];
double *** G3_tau_taup, *** tauDtau_G3_tau_taup, * transport_pts;
double * A1pts, *** A2pts, *** Bpts, *** Cpts;

//general functions
inline double Omega(double x)
{
	return (
		0.48*tanh(0.23*x) + (1.04 / M_PI) * atan(0.65*x)
	);
}

//equation of state and other thermodynamic relations
inline double P(double T, double mu)
{
	return(
		A4*pow(T, 4.0) + A2*T*T*mu*mu + A0*pow(mu, 4.0) - C0*T*T - B
	);
}

inline void set_phase_diagram_and_EOS_parameters()
{
	Nf = 2.0;											//number of massless flavors
	//T0 = 170.0;											//phase transition curve, T scale
	//mu0 = 1218.48;										//phase transition curve, mu scale
	T0 = 170.0 / hbarC;											//phase transition curve, T scale
	mu0 = 1218.48 / hbarC;										//phase transition curve, mu scale
	A4 = M_PI*M_PI*(16.0 + 10.5*Nf) / 90.0;				//coeff in P(T,mu) (i.e., EOS)
	A2 = Nf / 18.0;										//same
	A0 = Nf / (324.0 * M_PI * M_PI);					//same
	C0 = mu0*mu0*( A2 - 2.0*A0*mu0*mu0 / (T0*T0) );		//same
	B = 0.8 * pow(T0, 4.0);								//same

	//some parameters for the thermal conductivity critical enhancement
	//xi0 = 0.37;
	//qD = 1.0/xi0;
	//xibar0 = 0.69;
	qD = M_PI * T0;
	xi0 = 1.0 / qD;
	xibar0 = 0.701;
	RD = 1.05;
	exp_delta = 4.815;
	exp_gamma = 1.24;
	exp_nu = 0.63;
	etaBYs = 1.0 / (4.0 * M_PI);

	return;
}

// functions to guess seed values for T and mu,
// followed by functions to iteratively solve equations
// for T and mu

inline double guess_T(double tau)
{
	return (Ti * pow(taui / tau, 1.0/3.0));
}

inline double guess_mu(double tau)
{
	return (mui * pow(taui / tau, 1.0/3.0));
}

///////////////////////////////////////////////////////////
//two separate definitions of s and n each, for convenience
///////////////////////////////////////////////////////////
inline double s(double tau)
{
	return (si * taui / tau);
}

inline double n(double tau)
{
	return (ni * taui / tau);
}

inline double s(double T, double mu)
{
	return (
		-2.0 * C0 * T + 4.0 * A4 * pow(T, 3.0) + 2.0 * A2 * T * mu * mu
	);
}

inline double n(double T, double mu)
{
	return (
		2.0 * A2 * T * T * mu + 4.0 * A0 * pow(mu, 3.0)
	);
}
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

inline double w(double T, double mu)
{
	return (T * s(T, mu) + mu * n(T, mu));
}

inline void set_critical_point_parameters()
{
	//set critical point quantities
	//Tc = 160.0;
	//muc = 411.74;
	Tc = 160.0 / hbarC;
	muc = 411.74 / hbarC;
	Pc = P(Tc, muc);
	sc = s(Tc, muc);
	nc = n(Tc, muc);
	wc = w(Tc, muc);
	
	return;
}

inline double chi_TT(double T, double mu)
{
	return (
		2.0 * ( -C0 + 6.0 * A4 * T * T + A2 * mu * mu )
	);
}

inline double chi_Tmu(double T, double mu)
{
	return (
		4.0 * A2 * T * mu
	);
}

inline double chi_mumu(double T, double mu)
{
	return (
		2.0 * ( A2 * T * T + 6.0 * A0 * mu * mu )
	);
}


inline double xi(double T, double mu)
{
	double t1 = (1.0/3.0)*((exp_delta - 1.0)/(2.0 - exp_gamma)) * pow(abs((T/Tc) - 1.0), exp_gamma);
	double t2 = 5.0 * exp_delta * pow(abs((n(T,mu)/nc) - 1.0), exp_delta - 1.0);
	return (
		xibar0 * pow( t1+t2, -exp_nu/exp_gamma )		//fm^1
	);
}

inline double chi_B(double T, double mu)
{
	double t1 = (1.0/3.0)*((exp_delta - 1.0)/(2.0 - exp_gamma)) * pow(abs((T/Tc) - 1.0), exp_gamma);
	double t2 = 5.0 * exp_delta * pow(abs((n(T,mu)/nc) - 1.0), exp_delta - 1.0);
	return (
		( nc*nc / ( (exp_delta + 1.0) * Pc ) ) * pow(t1+t2, -1.0)		//MeV^2
	);
}

inline double cP(double T, double mu)
{
	double tp = (T/mu) * (2.0*A4*T*T + A2*mu*mu - C0) / (A2*T*T + 2.0*A0*mu*mu);
	double tm = 2.0*A2*T*mu / (A2*T*T + 6.0*A0*mu*mu);
	return (
		T * chi_B(T,mu) * (tp - tm) * (tp - tm)			//MeV^3
	);
}

inline double Delta_DT(double T, double mu)
{
	double xi_loc = xi(T, mu);
	return (
		RD * T * Omega(qD * xi_loc) / (6.0 * M_PI * etaBYs * s(T, mu) * xi_loc)
	);
}

inline double Delta_lambda(double T, double mu)
{
	return (cP(T,mu) * Delta_DT(T,mu));	//MeV^2
}

//speeds of sound for FTd-Green functions
inline double d(double T, double mu)
{
	return (
		2.0*A4*A2*pow(T, 4.0) + (12.0*A4*A0 - A2*A2)*mu*mu*T*T + 2.0*A2*A0*pow(mu, 4.0)
			- C0*(A2*T*T + 6.0*A0*mu*mu) / 3.0
	);
}

inline double vn2(double T, double mu)
{
	return (
		(1.0/3.0) - ( 2.0*C0*( A2*T*T + 6.0*A0*mu*mu ) ) / ( 9.0*d(T, mu) )
	);
}

inline double vs2(double T, double mu)
{
	return (
		(1.0/3.0) + ( ( 4.0*A2*C0*T*T ) / ( 9.0*d(T, mu) ) )
	);
}

inline double vsigma2(double T, double mu)
{
	double s_loc = s(T, mu);
	double n_loc = n(T, mu);
	double w_loc = w(T, mu);
	double vn2_loc = vn2(T, mu);
	double vs2_loc = vs2(T, mu);
	return (
		(T*s_loc*vn2_loc + mu*n_loc*vs2_loc) / w_loc
	);
}

////////////////////////////////////////////////////////////////////////////////
// Functions/stuff to solve for time-dependence of T and mu and/or Tf and muf
////////////////////////////////////////////////////////////////////////////////

struct rparams
{
	double tau;
};

int print_state (size_t iter, gsl_multiroot_fsolver * s)
{
	printf ("iter = %3u x = % .3f % .3f "
		"f(x) = % .3e % .3e\n",
		iter,
		gsl_vector_get (s->x, 0), 
		gsl_vector_get (s->x, 1),
		gsl_vector_get (s->f, 0), 
		gsl_vector_get (s->f, 1));
}

//compute final time-step Tf and muf

int input_get_Tf_and_muf_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	double tau_local = ((struct rparams *) params)->tau;

	const double x0 = gsl_vector_get (x, 0);	//T
	const double x1 = gsl_vector_get (x, 1);	//mu

	const double y0 = s(x0, x1)/n(x0, x1) - sPERn;	//defines fixed s/n curve
	const double y1 = P(x0,x1);						//defines P==0 curve

	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);

	return GSL_SUCCESS;
}

void compute_Tf_and_muf()
{
	const gsl_multiroot_fsolver_type *gsl_T;
	gsl_multiroot_fsolver *gsl_s;

	int status;
	size_t i, iter = 0;

	const size_t n = 2;
	struct rparams p = {taui};
	gsl_multiroot_function f = {&input_get_Tf_and_muf_f, n, &p};

	double x_init[2] = {Ti, mui};
	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set (x, 0, x_init[0]);
	gsl_vector_set (x, 1, x_init[1]);

	gsl_T = gsl_multiroot_fsolver_hybrids;
	gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, 2);
	gsl_multiroot_fsolver_set (gsl_s, &f, x);

	//print_state (iter, gsl_s);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate (gsl_s);

		//print_state (iter, gsl_s);

		if (status)   /* check if solver is stuck */
			break;

		status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	//printf ("status = %s\n", gsl_strerror (status));

	//finally, store results
	Tf = gsl_vector_get (gsl_s->x, 0);
	muf = gsl_vector_get (gsl_s->x, 1);

	gsl_multiroot_fsolver_free (gsl_s);
	gsl_vector_free (x);
	
	return;
}

//compute T and mu at each time step

int input_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	double tau_local = ((struct rparams *) params)->tau;

	const double x0 = gsl_vector_get (x, 0);	//T
	const double x1 = gsl_vector_get (x, 1);	//mu

	const double y0 = s(x0, x1) - s(tau_local);
	const double y1 = n(x0, x1) - n(tau_local);

	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);

	return GSL_SUCCESS;
}

void populate_T_and_mu_vs_tau()
{
	//both sections together
	for (int it = 0; it < n_tau_pts; ++it)
	{
		const gsl_multiroot_fsolver_type *gsl_T;
		gsl_multiroot_fsolver *gsl_s;

		int status;
		size_t i, iter = 0;

		const size_t n = 2;
		struct rparams p = {tau_pts[it]};
		gsl_multiroot_function f = {&input_f, n, &p};

		double x_init[2] = {guess_T(tau_pts[it]), guess_mu(tau_pts[it])};
		gsl_vector *x = gsl_vector_alloc (n);

		gsl_vector_set (x, 0, x_init[0]);
		gsl_vector_set (x, 1, x_init[1]);

		gsl_T = gsl_multiroot_fsolver_hybrids;
		gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, 2);
		gsl_multiroot_fsolver_set (gsl_s, &f, x);

		//print_state (iter, gsl_s);

		do
		{
			iter++;
			status = gsl_multiroot_fsolver_iterate (gsl_s);

			//print_state (iter, gsl_s);

			if (status)   /* check if solver is stuck */
				break;

			status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
		}
		while (status == GSL_CONTINUE && iter < 1000);

		//printf ("status = %s\n", gsl_strerror (status));

		//finally, store results
		T_pts[it] = gsl_vector_get (gsl_s->x, 0);
		mu_pts[it] = gsl_vector_get (gsl_s->x, 1);
		//cout << "full: " << T_pts[it] << "   " << mu_pts[it] << endl;

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
	}
	
	return;
}

void populate_T_and_mu_vs_tau_part2()
{
	//lower section
	for (int it = 0; it < n_tau_pts; ++it)
	{
		const gsl_multiroot_fsolver_type *gsl_T;
		gsl_multiroot_fsolver *gsl_s;

		int status;
		size_t i, iter = 0;

		const size_t n = 2;
		struct rparams p = {tau_pts_lower[it]};
		gsl_multiroot_function f = {&input_f, n, &p};

		double x_init[2] = {guess_T(tau_pts_lower[it]), guess_mu(tau_pts_lower[it])};
		gsl_vector *x = gsl_vector_alloc (n);

		gsl_vector_set (x, 0, x_init[0]);
		gsl_vector_set (x, 1, x_init[1]);

		gsl_T = gsl_multiroot_fsolver_hybrids;
		gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, 2);
		gsl_multiroot_fsolver_set (gsl_s, &f, x);

		//print_state (iter, gsl_s);

		do
		{
			iter++;
			status = gsl_multiroot_fsolver_iterate (gsl_s);

			//print_state (iter, gsl_s);

			if (status)   /* check if solver is stuck */
				break;

			status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
		}
		while (status == GSL_CONTINUE && iter < 1000);

		//printf ("status = %s\n", gsl_strerror (status));

		//finally, store results
		T_pts_lower[it] = gsl_vector_get (gsl_s->x, 0);
		mu_pts_lower[it] = gsl_vector_get (gsl_s->x, 1);
		//cout << "lower: " << T_pts_lower[it] << "   " << mu_pts_lower[it] << endl;

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
		//if (1) exit (0);
	}

	//upper section
	for (int it = 0; it < n_tau_pts; ++it)
	{
		const gsl_multiroot_fsolver_type *gsl_T;
		gsl_multiroot_fsolver *gsl_s;

		int status;
		size_t i, iter = 0;

		const size_t n = 2;
		struct rparams p = {tau_pts_upper[it]};
		gsl_multiroot_function f = {&input_f, n, &p};

		double x_init[2] = {guess_T(tau_pts_upper[it]), guess_mu(tau_pts_upper[it])};
		gsl_vector *x = gsl_vector_alloc (n);

		gsl_vector_set (x, 0, x_init[0]);
		gsl_vector_set (x, 1, x_init[1]);

		gsl_T = gsl_multiroot_fsolver_hybrids;
		gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, 2);
		gsl_multiroot_fsolver_set (gsl_s, &f, x);

		//print_state (iter, gsl_s);

		do
		{
			iter++;
			status = gsl_multiroot_fsolver_iterate (gsl_s);

			//print_state (iter, gsl_s);

			if (status)   /* check if solver is stuck */
				break;

			status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
		}
		while (status == GSL_CONTINUE && iter < 1000);

		//printf ("status = %s\n", gsl_strerror (status));

		//finally, store results
		T_pts_upper[it] = gsl_vector_get (gsl_s->x, 0);
		mu_pts_upper[it] = gsl_vector_get (gsl_s->x, 1);
		//cout << "upper: " << T_pts_upper[it] << "   " << mu_pts_upper[it] << endl;

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
	}
	
	return;
}

//a miscellaneous function to help split up tau integral
void break_up_integral(int nt, double &max_DL, double &tauc, double &width)
{
	double Delta_t = (tauf - taui - 1.0) / (double)nt;
	vector<double> Delta_lambda_vec;
	vector<double> tpts;

	//check Delta_lambda
	for (int it = 1; it < nt; ++it)
	{
		double t_loc = taui + 0.5 + (double)it * Delta_t;
		double T_loc = interpolate1D(tau_pts, T_pts, t_loc, n_tau_pts, 1, false);
		double mu_loc = interpolate1D(tau_pts, mu_pts, t_loc, n_tau_pts, 1, false);
		double DL = Delta_lambda(T_loc, mu_loc);
		Delta_lambda_vec.push_back(DL);
		tpts.push_back(t_loc);
		//cout << t_loc << "   " << DL << endl;
	}
	std::vector<double>::iterator result = std::max_element(Delta_lambda_vec.begin(), Delta_lambda_vec.end());
    int idx = std::distance(Delta_lambda_vec.begin(), result);
	max_DL = Delta_lambda_vec[idx];
	tauc = tpts[idx];
	width = 0.0;
	double den = 0.0;
	for (int ii = 0; ii < nt; ++ii)
	{
		width += pow(tpts[ii] - tauc, 2.0) * Delta_lambda_vec[ii];
		den += Delta_lambda_vec[ii];
	}
	width = sqrt(width / den);
	
	//cout << max_DL << "   " << tauc << "   " << width << endl;
	return;
}

inline void set_all_thermodynamic_points()
{
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double local_T = T_pts_lower[it];
		double local_mu = mu_pts_lower[it];
		all_T_pts[it] = local_T;
		all_mu_pts[it] = local_mu;		

		Delta_lambda_pts[it] = Delta_lambda(local_T, local_mu);
		vn2_pts[it] = vn2(local_T, local_mu);
		vs2_pts[it] = vs2(local_T, local_mu);
		vsigma2_pts[it] = vsigma2(local_T, local_mu);
		n_Tmu_pts[it] = n(local_T, local_mu);
		s_Tmu_pts[it] = s(local_T, local_mu);
		w_Tmu_pts[it] = w(local_T, local_mu);
	}
	for (int it = n_tau_pts; it < 2*n_tau_pts; ++it)
	{
		double local_T = T_pts_upper[it];
		double local_mu = mu_pts_upper[it];
		all_T_pts[it] = local_T;
		all_mu_pts[it] = local_mu;		
		
		Delta_lambda_pts[it] = Delta_lambda(local_T, local_mu);
		vn2_pts[it] = vn2(local_T, local_mu);
		vs2_pts[it] = vs2(local_T, local_mu);
		vsigma2_pts[it] = vsigma2(local_T, local_mu);
		n_Tmu_pts[it] = n(local_T, local_mu);
		s_Tmu_pts[it] = s(local_T, local_mu);
		w_Tmu_pts[it] = w(local_T, local_mu);
	}
	return;
}

inline void initialize_all(int chosen_trajectory, int particle_to_study)
{
	mui = muis[chosen_trajectory];
	mui /= hbarC;

	set_phase_diagram_and_EOS_parameters();

	//other constants
	m = (particle_to_study == 1) ? 139.57 / hbarC : 939.0 / hbarC;
	taui = 0.5;		//fm/c

	si = s(Ti, mui);
	ni = n(Ti, mui);
	sPERn = si / ni;

	compute_Tf_and_muf();

	sf = s(Tf, muf);
	tauf = si * taui / sf;

	set_critical_point_parameters();

	mu_proton = muf;
	mu_part = (particle_to_study == 1) ? mu_pion : mu_proton;
	//mu_part = muf;
	s_at_mu_part = s(Tf, mu_part);

    // initialize other parameters
    nu = 1.0 / (3.0 * M_PI);
	nuVB = nu;
    ds = 2.0;
	mByT = m / Tf;

	//set the susceptibilities
	double chi_mu_mu = chi_mumu(Tf, mu_part);
	double chi_T_mu = chi_Tmu(Tf, mu_part);
	double chi_T_T = chi_TT(Tf, mu_part);
	Delta = chi_mu_mu * chi_T_T - chi_T_mu * chi_T_mu;
	chi_tilde_mu_mu = chi_mu_mu / Delta;
	chi_tilde_T_mu = chi_T_mu / Delta;
	chi_tilde_T_T = chi_T_T / Delta;

	//set parameters for FTd-Green's functions
	a_at_tauf = alpha(Tf, muf);
	vs2_at_tauf = vs2(Tf, muf);
	vn2_at_tauf = vn2(Tf, muf);
	vsigma2_at_tauf = vsigma2(Tf, muf);

    // set up grid points for integrations
    xi_pts_minf_inf = new double [n_xi_pts];
    tau_pts_lower = new double [n_tau_pts];
    tau_pts_upper = new double [n_tau_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts_lower = new double [n_tau_pts];
    tau_wts_upper = new double [n_tau_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];
	x_pts = new double [n_x_pts];
	x_wts = new double [n_x_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -xi_infinity, xi_infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, 0.0, k_infinity, k_pts, k_wts);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);
    tmp = gauss_quadrature(n_x_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

	T_pts_lower = new double [n_tau_pts];
	T_pts_upper = new double [n_tau_pts];
	T_pts = new double [n_tau_pts];
	mu_pts_lower = new double [n_tau_pts];
	mu_pts_upper = new double [n_tau_pts];
	mu_pts = new double [n_tau_pts];

	Delta_lambda_pts = new double [2*n_tau_pts];
	vn2_pts = new double [2*n_tau_pts];
	vs2_pts = new double [2*n_tau_pts];
	vsigma2_pts = new double [2*n_tau_pts];
	n_Tmu_pts = new double [2*n_tau_pts];
	s_Tmu_pts = new double [2*n_tau_pts];
	w_Tmu_pts = new double [2*n_tau_pts];

	//computes tau-dependence of T and mu for remainder of calculation
	populate_T_and_mu_vs_tau();

	int nt = 1000000;
	double max_DL, tauc, width;
	break_up_integral(nt, max_DL, tauc, width);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauc, tau_pts_lower, tau_wts_lower);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, tauc, tauf, tau_pts_upper, tau_wts_upper);
	populate_T_and_mu_vs_tau_part2();

	set_all_thermodynamic_points();

	//allocate and initialize some other needed arrays
	G3_tau_taup = new double ** [n_k_pts];
	tauDtau_G3_tau_taup = new double ** [n_k_pts];
	transport_pts = new double [2*n_tau_pts];
	A1pts = new double [2*n_tau_pts];
	A2pts = new double ** [n_k_pts];
	Bpts = new double ** [n_k_pts];
	Cpts = new double ** [n_k_pts];
	//
	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		G3_tau_taup[ik] = new double * [2*n_tau_pts];
		tauDtau_G3_tau_taup[ik] = new double * [2*n_tau_pts];
		A2pts[ik] = new double * [2*n_tau_pts];
		Bpts[ik] = new double * [2*n_tau_pts];
		Cpts[ik] = new double * [2*n_tau_pts];
		//
		for (int it = 0; it < 2*n_tau_pts; ++it)
		{
			G3_tau_taup[ik][it] = new double [2*n_tau_pts];
			tauDtau_G3_tau_taup[ik][it] = new double [2*n_tau_pts];
			A2pts[ik][it] = new double [2*n_tau_pts];
			Bpts[ik][it] = new double [2*n_tau_pts];
			Cpts[ik][it] = new double [2*n_tau_pts];
			//
			for (int itp = 0; itp < 2*n_tau_pts; ++itp)
			{
				G3_tau_taup[ik][it][itp] = 0.0;
				tauDtau_G3_tau_taup[ik][it][itp] = 0.0;
				A2pts[ik][it][itp] = 0.0;
				Bpts[ik][it][itp] = 0.0;
				Cpts[ik][it][itp] = 0.0;
			}
		}
	}

	return;
}

// End of file

#endif
