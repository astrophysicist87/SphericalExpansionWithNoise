#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
//#include <algorithm>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>

using namespace std;

struct params1D
{
	double conjugate;	//u <--> k
	double (*function)(double);
};

struct params2D
{
	double k1, k2;
};

inline double gt_k1_gt_k2 (double x, void * params_ptr)
{
	params2D local_parameters = *(params2D *) (params_ptr);
	double k1 = local_parameters.k1;
	double k2 = local_parameters.k2;
	return ( gsl_sf_conicalP_1(k1, x) * gsl_sf_conicalP_1(k2, x) / sqrt(x*x - 1.0) );
}

inline double MFkernel(double u, void * params_ptr)
{
	params1D local_parameters = *(params1D *) (params_ptr);
	double local_conjugate = local_parameters.conjugate;
	return (gsl_sf_conicalP_0(local_conjugate, u) * (local_parameters.function)(u));
}

inline double inverseMFkernel(double k, void * params_ptr)
{
	params1D local_parameters = *(params1D *) (params_ptr);
	double local_conjugate = local_parameters.conjugate;
	return (k * tanh(M_PI * k) * gsl_sf_conicalP_0(k, local_conjugate) * (local_parameters.function)(k));
}
}

inline double mehler_fock_transform(double k, double (*f)(double))
{
	gsl_function F;
	F.function = &MFkernel;

	struct params1D my_params;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	double result, error;
	my_params.conjugate = k;
	my_params.function = *f;

	F.params = &my_params;

	gsl_integration_qagiu (&F, 1.0, 0, 1e-7, 1000, w, &result, &error); 

	gsl_integration_workspace_free (w);
	return (result);
}

inline double mehler_fock_inverse_transform(double u, double (*f)(double))
{
	gsl_function F;
	F.function = &MFkernel;

	struct params1D my_params;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	double result, error;
	my_params.conjugate = u;
	my_params.function = *f;

	F.params = &my_params;

	gsl_integration_qagiu (&F, 1.0, 0, 1e-7, 1000, w, &result, &error); 

	gsl_integration_workspace_free (w);
	return (result);
}

inline void mehler_fock_transform(double * k_pts, double (*f)(double), double * results, int n_pts)
{
	gsl_function F;
	F.function = &MFkernel;

	struct params1D my_params;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	double result, error;
	my_params.function = *f;

	for (int ipt = 0; ipt < n_pts; ++ipt)
	{
		my_params.conjugate = k_pts[ipt];
		F.params = &my_params;
		gsl_integration_qagiu (&F, 1.0, 0, 1e-7, 1000, w, &result, &error);
		results[ipt] = result;
	}

	gsl_integration_workspace_free (w);

	return;
}

inline void mehler_fock_inverse_transform(double * u_pts, double (*f)(double), double * results, int n_pts)
{
	gsl_function F;
	F.function = &MFkernel;

	struct params1D my_params;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	double result, error;
	my_params.function = *f;

	for (int ipt = 0; ipt < n_pts; ++ipt)
	{
		my_params.conjugate = u_pts[ipt];
		F.params = &my_params;
		gsl_integration_qagiu (&F, 1.0, 0, 1e-7, 1000, w, &result, &error);
		results[ipt] = result;
	}

	gsl_integration_workspace_free (w);

	return (result);
}

void compute_legendre_integral(vector<vector<double> > array, vector<double> k_pts_arr)
{
	int n_k_pts = k_pts_arr.size();

	//assume array has been properly defined elsewhere
	double k1 = k_pts_arr[ik1];
	double k2 = k_pts_arr[ik2];
	gsl_function F;
	F.function = &f;	//should be gt_k1_gt_k2, I think

	params2D my_params;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		double result, error;
		my_params.k1 = k1;
		my_params.k2 = k2;

		F.params = &my_params;

		gsl_integration_qagiu (&F, 1.0, 0, 1e-7, 1000, w, &result, &error); 

		array[ik1][ik2] = result;
	}

	gsl_integration_workspace_free (w);

	return;
}

void invert_2D_MFspace(vector<vector<double> > array, vector<double> k_pts_arr)
{
	int n_k_pts = k_pts_arr.size();

	//assume array has been properly defined elsewhere
	double k1 = k_pts_arr[ik1];
	double k2 = k_pts_arr[ik2];
	gsl_function F;
	F.function = &f;	//should be gt_k1_gt_k2, I think

	params2D my_params;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		double result, error;
		my_params.k1 = k1;
		my_params.k2 = k2;

		F.params = &my_params;

		gsl_integration_qagiu (&F, 1.0, 0, 1e-7, 1000, w, &result, &error); 

		array[ik1][ik2] = result;
	}

	gsl_integration_workspace_free (w);

	return;
}

void set_Q_X_k(vector<vector<double> > result, vector<double> k_pts_arr, vector<double> u_pts_arr)
{
	int n_k_pts = k_pts_arr.size();
	int n_u_pts = u_pts_arr.size();

	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		double u_loc = u_pts_arr[iu];
		result[0][ik] = gsl_sf_conicalP_0(k_pts_arr[ik], u_pts_arr[iu]);
		result[1][ik] = sqrt(u_loc*u_loc-1.0)*gsl_sf_conicalP_cyl_reg(1, k_pts_arr[ik], u_pts_arr[iu]);
		result[2][ik] = result[0][ik];
	}

	return;
}

//End of file
