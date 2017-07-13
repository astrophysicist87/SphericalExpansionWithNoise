#ifndef EMISSION_DENSITY_H
#define EMISSION_DENSITY_H

/////////////////////////////////////////////////////
//**********************************
// General outline and methodology:
//**********************************
// Emission density S is first expanded
// schematically into series of fluctuations:
//		S = S^0 + delta_k * S_k^1
//			+ 0.5 * delta_k * delta_l * S_{kl}^2 + ...
//
// Each term depends on rapidity and
// thermodynamic quantities at freeze-out
//
// S_k^1 and S_{kl}^2 built up out of
// derivatives of S(x)
//
// Defined S == S_A * S_B to simplify derivatives
//
// All derivatives taken explicitly beforehand,
// and then written out in terms of appropriate
// functions.  All derivative coefficients evaluated
// at zero fluctuations.
//
// delta components are indexed as follows:
// delta T - 0, omega - 1, delta mu - 2
/////////////////////////////////////////////////////


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

inline double S_A(double xi)
{
	return ( exp( mT * ch_xi - muf) / Tf ) );
}

inline double S_B(double xi)
{
	double a1 = mT*ch_xi;
	double sh_xi = sinh(xi);
	double a2 = pT*sh_xi;

	double z = pT*sh_xi/Tf;
	double sh_z = sinh(z);
	double ch_z = cosh(z);

	return ( a1*sh_z/z - a2*( sh_z - z*ch_z ) / (z*z) );
}

inline double S0(double xi)
{
	return ( S_A(xi)*S_B(xi) );
}


// need these next two functions for
// taking derivatives of S_B
inline double S_B_a(double z, double a1, double a2)
{
	double z2 = z*z;
	double z3 = z2*z;

	return (
			( z*( a1*z-2.0*a2 )*cosh(z)
				+ ( a2*( z2+2.0 )-a1*z )*sinh(z) )
			/ z3 );
}

inline double S_B_b(double z, double a1, double a2)
{
	double z2 = z*z;
	double z3 = z2*z;
	double z4 = z3*z;

	return (
				z*(a2*(z2+6.0)-2.0*a1*z)*cosh(z)
					+ (z2+2.0)*(a1*z-3.0*a2)*sinh(z)
			);
}

vector<double> void set_first_derivatives_vector(double xi)
{
	double sh_xi = sinh(x), ch_xi = cosh(xi);
	vector<double> result;
	result.push_back( (-muf + mT*ch_xi)/Tf2 );
	result.push_back( -mT*sh_xi/Tf );
	result.push_back( 1.0/Tf );

	return (result);
}

vector<vector<double> > void set_second_derivatives_array(double xi)
{
	vector<vector<double> > result;

	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	double sh_xi = sinh(x), ch_xi = cosh(xi);

	vector<double> row1;
	row1.push_back( ((muf - mT*ch_xi)*(muf + 2.0*Tf - mT*ch_xi))/Tf4 );
	row1.push_back( (mT*(muf + Tf - mT*ch_xi)*sh_xi)/Tf3 );
	row1.push_back( -((muf + Tf - mT*ch_xi)/Tf3) );
	result.push_back( row1 );

	vector<double> row2;
	row2.push_back( (mT*(muf + Tf - mT*ch_xi)*sh_xi)/Tf3 );
	row2.push_back( (mT*(-(Tf*ch_xi) + mT*sh_xi*sh_xi))/Tf2 );
	row2.push_back( -((mT*sh_xi)/Tf2) );
	result.push_back( row2 );

	vector<double> row3;
	row3.push_back( -((muf + Tf - mT*ch_xi)/Tf3) );
	row3.push_back( -((mT*sh_xi)/Tf2) );
	row3.push_back( 1.0/Tf2 );
	result.push_back( row3 );

	return (result);
}











// End of file

#endif
