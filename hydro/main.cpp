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
#include "emission_density.h"
#include "gauss_quadrature.h"

bool do_1p_calc;
bool do_HBT_calc;
bool scale_out_y_dependence = false;
const int particle_to_study = 2;	//1 is pion, 2 is proton

const double hbarC = 197.33;
const double k_infinity = 10.0;
const double u_infinity = 10.0;
const int n_Dy = 1;

const double mu_pion = 1.e-3 / hbarC;
const int n_u_pts = 50;
const int n_k_pts = 50;
const int n_x_pts = 51;
const int n_tau_pts = 50;

double muis[3] = {420.0, 620.0, 820.0};

double bar_w_ij(int iu, int iup)
{
	return 1.0;
}

int main(int argc, char *argv[])
{
	//set parameters corresponding to all different possible trajectories
	Ti = 250.0 / hbarC;		//initial trajectory temperature
	current_itau = 0;

	int chosen_trajectory = atoi(argv[1]);
	chosen_trajectory--;		//use for array indexing throughout

	//do full calculation
	complex<double> mean_delta_R2ij = get_mean_delta_R2ij(chosen_trajectory, particle_to_study);

	return 0;
}
