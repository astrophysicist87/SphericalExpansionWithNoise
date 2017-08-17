#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>

using namespace std;

#include "legendre.h"

inline double gt_k1_gt_k2 (double x, void * params_ptr)
{
	params2D local_parameters = *(params2D *) (params_ptr);
	double k1 = local_parameters.k1;
	double k2 = local_parameters.k2;
	return ( gsl_sf_conicalP_1(k1, x) * gsl_sf_conicalP_1(k2, x) / sqrt(x*x - 1.0) );
}


//CHECKED!!!
void compute_legendre_integral(vector<vector<double> > & array, vector<double> k_pts_arr)
{
	int n_k_pts = k_pts_arr.size();

	//assume array has been properly defined elsewhere
	gsl_function F;
	F.function = &gt_k1_gt_k2;

	params2D my_params;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	for (int ik1 = 0; ik1 < n_k_pts; ++ik1)
	for (int ik2 = 0; ik2 < n_k_pts; ++ik2)
	{
		double result, error;
		double k1 = k_pts_arr[ik1];
		double k2 = k_pts_arr[ik2];
		my_params.k1 = k1;
		my_params.k2 = k2;

		F.params = &my_params;

		gsl_integration_qagiu (&F, 1.0, 0, 1e-7, 1000, w, &result, &error); 

		array[ik1][ik2] = result;
//cout << "Check compute_legendre_integral(): " << k_pts_arr[ik1] << "   " << k_pts_arr[ik2] << "   " << result << endl;
	}

	gsl_integration_workspace_free (w);

	return;
}

//CHECKED!!!
void set_Q_X_k(vector<vector<vector<double> > > & result, vector<double> k_pts_arr, vector<double> u_pts_arr)
{
	int n_k_pts = k_pts_arr.size();
	int n_u_pts = u_pts_arr.size();

	for (int iu = 0; iu < n_u_pts; ++iu)
	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		double u_loc = u_pts_arr[iu];
		result[0][ik][iu] = gsl_sf_conicalP_0(k_pts_arr[ik], u_pts_arr[iu]);
		result[1][ik][iu] = sqrt(u_loc*u_loc-1.0)*gsl_sf_conicalP_cyl_reg(1, k_pts_arr[ik], u_pts_arr[iu]);
		result[2][ik][iu] = result[0][ik][iu];
//cout << "Check set_Q_X_k(): " << u_loc << "   " << k_pts_arr[ik] << "   " << result[0][ik][iu] << "   " << result[1][ik][iu] << endl;
	}

	return;
}

//End of file
