#ifndef LIB_H
#define LIB_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <string>

#include <gsl/gsl_integration.h>

#define VNAME(x) #x
#define VDUMP(x) std::cout << #x << " " << x << std::endl

using namespace std;

const std::complex<double> i(0, 1);
inline double cot(double x){return (cos(x)/sin(x));}

double integrate_1D(double (*f)(double), double * xpts, double * xwts, int nx);
complex<double> integrate_1D_FT(double (*f)(double), double * xpts, double * xwts, int nx, double k);
double integrate_2D(double (*f)(double, double), double * xpts, double * ypts, double * xwts, double * ywts, int nx, int ny);

//main interpolation routine
double interpolate1D(vector<double> x, vector<double> y, double x0, int kind, bool uniform_spacing,
			bool returnflag = false, double default_return_value = 0.0);

//subsidiary interpolation routines
double interpLinearDirect(vector<double> x, vector<double> y, double x0,
				bool returnflag = false, double default_return_value = 0.0);
double interpLinearNondirect(vector<double> x, vector<double> y, double x0,
				bool returnflag = false, double default_return_value = 0.0);
double interpCubicDirect(vector<double> x, vector<double> y, double x0,
				bool returnflag = false, double default_return_value = 0.0);
double interpCubicNonDirect(vector<double> x, vector<double> y, double x0,
				bool returnflag = false, double default_return_value = 0.0);
long binarySearch(vector<double> A, double value,
				bool skip_out_of_range = true, bool verbose = false);
void solve_cubic_equation(double b, double c, double d, complex<double> &r1, complex<double> &r2, complex<double> &r3);

template <typename T>
void create_matrix_2D(
	vector<vector<T> > * matrix_to_create,
	int dim1, int dim2, bool verbose = false)
{
	/*for (int i = 0; i < dim1; ++i)
	{
		vector<T> tmp(dim2);
		(*matrix_to_create).push_back( tmp );
	}*/
	(*matrix_to_create) = vector<vector<T> > (dim1);
	for (int i = 0; i < dim1; ++i)
	{
		(*matrix_to_create)[i] = vector<T> (dim2);
		for (int j = 0; j < dim2; ++j)
			(*matrix_to_create)[i][j] = T(0.0);
	}
	if (verbose)
	{
		/*cout << "create_matrix_2D():" << endl;
		cout << dim1 << "   " << (*matrix_to_create).size() << endl;
		cout << dim2 << "   " << (*matrix_to_create)[0].size() << endl;*/
		cout << "create_matrix_2D(): (" << dim1 << " x " << dim2
				<< ") <--> (" << (*matrix_to_create).size() << " x " << (*matrix_to_create)[0].size() << ")" << endl;
	}
	return;
}

template <typename T>
void create_matrix_3D(
	vector<vector<vector<T> > > * matrix_to_create,
	int dim1, int dim2, int dim3, bool verbose = false)
{
	/*for (int i = 0; i < dim1; ++i)
	{
		vector<vector<T> > tmp2(dim2);
		for (int j = 0; j < dim2; ++j)
		{
			vector<T> tmp3(dim3);
			tmp2.push_back( tmp3 );
		}
		(*matrix_to_create).push_back( tmp2 );
	}*/
	(*matrix_to_create) = vector<vector<vector<T> > > (dim1);
	for (int i = 0; i < dim1; ++i)
	{
		(*matrix_to_create)[i] = vector<vector<T> > (dim2);
		for (int j = 0; j < dim2; ++j)
		{
			(*matrix_to_create)[i][j] = vector<T> (dim3);
			for (int k = 0; k < dim3; ++k)
				(*matrix_to_create)[i][j][k] = T(0.0);
		}
	}
	if (verbose)
	{
		/*cout << "create_matrix_3D():" << endl;
		cout << dim1 << "   " << (*matrix_to_create).size() << endl;
		cout << dim2 << "   " << (*matrix_to_create)[0].size() << endl;
		cout << dim3 << "   " << (*matrix_to_create)[0][0].size() << endl;*/
		cout << "create_matrix_3D(): (" << dim1 << " x " << dim2 << " x "  << dim3
				<< ") <--> (" << (*matrix_to_create).size() << " x " << (*matrix_to_create)[0].size() << " x " << (*matrix_to_create)[0][0].size() << ")" << endl;
	}
	return;
}

template <typename T>
void create_matrix_4D(
	vector<vector<vector<vector<T> > > > * matrix_to_create,
	int dim1, int dim2, int dim3, int dim4, bool verbose = false)
{
	/*for (int i = 0; i < dim1; ++i)
	{
		vector<vector<vector<T> > > tmp2(dim2);
		for (int j = 0; j < dim2; ++j)
		{
			vector<vector<T> > tmp3(dim3);
			for (int k = 0; k < dim3; ++k)
			{
				vector<T> tmp4(dim4);
				tmp3.push_back( tmp4 );
			}
			tmp2.push_back( tmp3 );
		}
		(*matrix_to_create).push_back( tmp2 );
	}*/
	(*matrix_to_create) = vector<vector<vector<vector<T> > > > (dim1);
	for (int i = 0; i < dim1; ++i)
	{
		(*matrix_to_create)[i] = vector<vector<vector<T> > > (dim2);
		for (int j = 0; j < dim2; ++j)
		{
			(*matrix_to_create)[i][j] = vector<vector<T> > (dim3);
			for (int k = 0; k < dim3; ++k)
			{
				(*matrix_to_create)[i][j][k] = vector<T> (dim4);
				for (int l = 0; l < dim4; ++l)
					(*matrix_to_create)[i][j][k][l] = T(0.0);
			}
		}
	}
	if (verbose)
	{
		/*cout << "create_matrix_4D():" << endl;
		cout << dim1 << " <--> " << (*matrix_to_create).size() << endl;
		cout << dim2 << " <--> " << (*matrix_to_create)[0].size() << endl;
		cout << dim3 << " <--> " << (*matrix_to_create)[0][0].size() << endl;
		cout << dim3 << " <--> " << (*matrix_to_create)[0][0][0].size() << endl;*/
		cout << "create_matrix_4D(): (" << dim1 << " x " << dim2 << " x "  << dim3 << " x "  << dim4
				<< ") <--> (" << (*matrix_to_create).size() << " x " << (*matrix_to_create)[0].size() << " x " << (*matrix_to_create)[0][0].size() << " x " << (*matrix_to_create)[0][0][0].size() << ")" << endl;
	}
	return;
}
		( ( ( ( ( ipt * n_pphi_pts + ipphi ) * qtnpts + iqt ) * qxnpts + iqx ) * qynpts + iqy ) * qznpts + iqz ) * 2 + itrig
inline long int indexer2D(int i, int j, int ni, int nj)
{
	return (
		i * nj + nj
		);
}

inline long int indexer3D(int i, int j, int k, int ni, int nj, int nk)
{
	return (
		( i * nj + j ) * nk + k
		);
}

inline long int indexer4D(int i, int j, int k, int l, int ni, int nj, int nk, int nl)
{
	return (
		( ( i * nj + j ) * nk + k ) * nl + l
		);
}

#endif
