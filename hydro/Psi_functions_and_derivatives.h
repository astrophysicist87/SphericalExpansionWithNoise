#ifndef PSI_FUNCTIONS_AND_DERIVATIVES_H
#define PSI_FUNCTIONS_AND_DERIVATIVES_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <cmath>

#include <gsl/gsl_sf_bessel.h>

using namespace std;

#include "lib.h"
#include "defs1.h"

////////////////////////////////////////////////////////////////////
// Psi function stuff
////////////////////////////////////////////////////////////////////
inline double Psi_A(double u)
{
	double prefactor = -ds*tauf*tauf*tauf / (8.0*M_PI*M_PI*M_PI);
	double shxi = sqrt(u*u-1.0);
	double arg = ( muf - mT*u ) / Tf;	//<=0
	return (
			prefactor * shxi * shxi * exp( arg )
			);
}

inline double Psi_B0(double a1, double a2, double z)
{
	return (
			4.0*M_PI*( ( a1 * z + a2 ) * sinh(z) - a2 * z * cosh(z) ) / (z*z)
			);
}

inline double Psi_B1(double a1, double a2, double z)
{
	double I0 = gsl_sf_bessel_I0(0.5*z);
	double I1 = gsl_sf_bessel_I1(0.5*z);
	return (
			M_PI*M_PI*( 2.0*(a1*z+a2)I1*I0 - a2*z*(I0*I0+I1*I1) ) / z
			);
}

inline double Psi_B2(double a1, double a2, double z)
{
	return (
			4.0*M_PI*( (a1*z+2.0*a2)*(z*sinh(z)+1)
						- (a1*z + a2*(z*z+2.0)) * cosh(z) ) / (z*z*z)
			);
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
inline double Psi_0(double u)
{
	double a1 = mT*u;
	double a2 = pT*sqrt(u*u-1.0);
	return ( Psi_A(xi) * Psi_B0(a1, a2, a2/Tf) );
}

inline double Psi_1(double u)
{
	double a1 = mT*u;
	double a2 = pT*sqrt(u*u-1.0);
	return ( Psi_A(xi) * Psi_B1(a1, a2, a2/Tf) );
}

inline double Psi_2(double u)
{
	double a1 = mT*u;
	double a2 = pT*sqrt(u*u-1.0);
	return ( Psi_A(xi) * Psi_B2(a1, a2, a2/Tf) );
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

//Functions useful for computing derivatives of Psi functions
inline double Psi_a_B0(double a1, double a2, double z)
{
	return (
				-4.0*M_PI * ( (a1*z + a2*(z*z+2.0))*sinh(z) - z*cosh(z)*(a1*z+2.0*a2) ) / (z*z*z)
			);
}

inline double Psi_b_B0(double a1, double a2, double z)
{
	return (
				4.0*M_PI * ( (z*z+2.0)*(a1*z+3.0*a2)*sinh(z) - z*cosh(z)*(2.0*a1*z+a2*(z*z+6.0)) ) / (z*z*z*z)
			);
}

inline double Psi_a_B1(double a1, double a2, double z)
{
	double I0 = gsl_sf_bessel_I0(0.5*z);
	double I1 = gsl_sf_bessel_I1(0.5*z);
	double I2 = gsl_sf_bessel_In(2, 0.5*z);
	return (
				M_PI*M_PI*( z*I0*I0*(a1*z+a2) + z*I1*I1*(a1*z+3.0*a2)
							-2.0*I1*I0*(a1*z+a2* ( z*z+2.0 ) ) ) / (z*z)
			);
}

inline double Psi_b_B1(double a1, double a2, double z)
{
	return (
				-M_PI*M_PI*( z*I0*I0*(a1*z+a2*(z*z+3.0)) + z*I1*I1*(3.0*a1*z+a2*(z*z+11.0))
							 - 2.0*(z*z+2.0) * I1*I0*(a1*z+3*a2) ) / (z*z*z)
			);
}

inline double Psi_a_B2(double a1, double a2, double z)
{
	return (
				4.0*M_PI * ( ( (z*z+2.0)*cosh(z) - 2.0 )*(a1*z+3.0*a2)
							- z*sinh(z)*(2.0*a1*z + a2*(z*z+6.0)) ) / (z*z*z*z)
			);
}

inline double Psi_b_B2(double a1, double a2, double z)
{
	return (
				4.0*M_PI * ( ( z*sinh(z)*(z*z+6.0) + 6.0 )*(a1*z+4.0*a2)
							- cosh(z) * ( 3.0*a1*z*(z*z+2.0) + a2*(z*z*z*z + 12.0*z*z + 24.0) ) ) / (z*z*z*z*z)
			);
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

//derivatives of components Psi functions
void set_PsiA_first_derivatives_vector(int iu)
{
	double u = u_pts[iu];

	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double prefactor = ds*tauf*tauf*tauf*sh_xi*sh_xi*exp((muf - mT*u)/Tf) / (8.0*M_PI*M_PI*M_PI);
	double lambda = 1.0 + (muf - mT*u)/Tf;
	double s_tilde = sf/Tf;

	dPsiA_dX[iu][0] = -prefactor*s_tilde*((lambda - 1.0)*chi_tilde_mu_mu + chi_tilde_T_mu);
	dPsiA_dX[iu][1] = -prefactor*mT*sh_xi/Tf;
	dPsiA_dX[iu][2] = prefactor*s_tilde*((lambda - 1.0)*chi_tilde_T_mu + chi_tilde_T_T);

	return;
}

void set_PsiA_second_derivatives_array(int iu)
{
	double u = u_pts[iu];

	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;
	double prefactor = ds*tauf*tauf*tauf*sh_xi*sh_xi*exp((muf - mT*u)/Tf) / (8.0*M_PI*M_PI*M_PI);
	double lambda = 1.0 + (muf - mT*u)/Tf;

	dPsiA_dX_dY[iu][0][0] = prefactor*s_tilde*s_tilde*( (lambda*lambda-1.0)*chi_tilde_mu_mu*chi_tilde_mu_mu
													+ 2.0*chi_tilde_mu_mu*chi_tilde_T_mu
													+ chi_tilde_T_mu*chi_tilde_T_mu );
	dPsiA_dX_dY[iu][0][1] = prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu);
	dPsiA_dX_dY[iu][0][2] = -prefactor*s_tilde*s_tilde*(
						chi_tilde_T_mu*((lambda*lambda - 1.0)*chi_tilde_mu_mu + lambda*chi_tilde_T_mu)
						 + chi_tilde_T_T*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu)
					);

	dPsiA_dX_dY[iu][1][0] = prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu);
	dPsiA_dX_dY[iu][1][1] = prefactor*(mT/Tf)*((mT/Tf)*sh_xi*sh_xi-ch_xi);
	dPsiA_dX_dY[iu][1][2] = -prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_T_mu+chi_tilde_T_T);

	dPsiA_dX_dY[iu][2][0] = -prefactor*s_tilde*s_tilde*(
						chi_tilde_T_mu*((lambda*lambda - 1.0)*chi_tilde_mu_mu + lambda*chi_tilde_T_mu)
						 + chi_tilde_T_T*(lambda*chi_tilde_mu_mu+chi_tilde_T_mu)
					);
	dPsiA_dX_dY[iu][2][1] = -prefactor*(mT/Tf)*s_tilde*sh_xi*(lambda*chi_tilde_T_mu+chi_tilde_T_T);
	dPsiA_dX_dY[iu][2][2] = prefactor*s_tilde*s_tilde*( (lambda*lambda-1.0)*chi_tilde_T_mu*chi_tilde_T_mu
													+ 2.0*chi_tilde_T_mu*chi_tilde_T_T
													+ chi_tilde_T_T*chi_tilde_T_T );

	return;
}

void set_PsiBk_first_derivatives_vector(int iu)
{
	double u = u_pts[iu];

	//double sh_xi = sinh(x), ch_xi = cosh(xi);
	double Tf2 = Tf*Tf;
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;

	double z_loc = pT*sh_xi/Tf;
	double local_PsiB0a = Psi_a_B0(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB1a = Psi_a_B1(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB2a = Psi_a_B2(z_loc, mT*ch_xi, pT*sh_xi);
	double kappa_S = (pT/Tf)*sh_xi, kappa_C = (pT/Tf)*ch_xi;

	double dz_dv0 = -s_tilde*kappa_S*chi_tilde_mu_mu;
	double dz_dv1 = kappa_C;
	double dz_dv2 = s_tilde*kappa_S*chi_tilde_T_mu;

	dPsiB0_dX[iu][0] = local_PsiB0a*dz_dv0;
	dPsiB0_dX[iu][1] = local_PsiB0a*dz_dv1;
	dPsiB0_dX[iu][2] = local_PsiB0a*dz_dv2;

	dPsiB1_dX[iu][0] = local_PsiB1a*dz_dv0;
	dPsiB1_dX[iu][1] = local_PsiB1a*dz_dv1;
	dPsiB1_dX[iu][2] = local_PsiB1a*dz_dv2;

	dPsiB2_dX[iu][0] = local_PsiB2a*dz_dv0;
	dPsiB2_dX[iu][1] = local_PsiB2a*dz_dv1;
	dPsiB2_dX[iu][2] = local_PsiB2a*dz_dv2;

	return ;
}

void set_PsiBk_second_derivatives_array(int iu)
{
	double Tf2 = Tf*Tf;
	double Tf3 = Tf2*Tf;
	double Tf4 = Tf3*Tf;
	//double sh_xi = sinh(xi), ch_xi = cosh(xi);
	double u = u_pts[iu];
	double sh_xi = sqrt(u*u-1.0), ch_xi = u;

	double z_loc = pT*sh_xi/Tf;

	double local_PsiB0a = Psi_a_B0(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB0b = Psi_b_B0(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB1a = Psi_a_B1(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB1b = Psi_b_B1(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB2a = Psi_a_B2(z_loc, mT*ch_xi, pT*sh_xi);
	double local_PsiB2b = Psi_b_B2(z_loc, mT*ch_xi, pT*sh_xi);

	double dz_dv0 = -s_tilde*kappa_S*chi_tilde_mu_mu;
	double dz_dv1 = kappa_C;
	double dz_dv2 = s_tilde*kappa_S*chi_tilde_T_mu;

	double dz_dv0_dv0 = 2.0*s_tilde*s_tilde*kappa_S*chi_tilde_mu_mu*chi_tilde_mu_mu;
	double dz_dv0_dv1 = -s_tilde*kappa_C*chi_tilde_mu_mu;
	double dz_dv0_dv2 = -2.0*s_tilde*s_tilde*kappa_S*chi_tilde_T_mu*chi_tilde_mu_mu;
	double dz_dv1_dv0 = dz_dv0_dv1;
	double dz_dv1_dv1 = kappa_S;
	double dz_dv1_dv2 = s_tilde*kappa_C*chi_tilde_T_mu;
	double dz_dv2_dv0 = dz_dv0_dv2;
	double dz_dv2_dv1 = dz_dv1_dv2;
	double dz_dv2_dv2 = 2.0*s_tilde*s_tilde*kappa_S*chi_tilde_T_mu*chi_tilde_T_mu;

	//second derivatives of Psi_B0
	dPsiB0_dX_dY[iu][0][0] = local_PsiB0a*dz_dv0_dv0 + local_PsiB0b*dz_dv0*dz_dv0;
	dPsiB0_dX_dY[iu][0][1] = local_PsiB0a*dz_dv0_dv1 + local_PsiB0b*dz_dv0*dz_dv1;
	dPsiB0_dX_dY[iu][0][2] = local_PsiB0a*dz_dv0_dv2 + local_PsiB0b*dz_dv0*dz_dv2;

	dPsiB0_dX_dY[iu][1][0] = local_PsiB0a*dz_dv1_dv0 + local_PsiB0b*dz_dv1*dz_dv0;
	dPsiB0_dX_dY[iu][1][1] = local_PsiB0a*dz_dv1_dv1 + local_PsiB0b*dz_dv1*dz_dv1;
	dPsiB0_dX_dY[iu][1][2] = local_PsiB0a*dz_dv1_dv2 + local_PsiB0b*dz_dv1*dz_dv2;

	dPsiB0_dX_dY[iu][2][0] = local_PsiB0a*dz_dv2_dv0 + local_PsiB0b*dz_dv2*dz_dv0;
	dPsiB0_dX_dY[iu][2][1] = local_PsiB0a*dz_dv2_dv1 + local_PsiB0b*dz_dv2*dz_dv1;
	dPsiB0_dX_dY[iu][2][2] = local_PsiB0a*dz_dv2_dv2 + local_PsiB0b*dz_dv2*dz_dv2;

	//second derivatives of Psi_B1
	dPsiB1_dX_dY[iu][0][0] = local_PsiB1a*dz_dv0_dv0 + local_PsiB1b*dz_dv0*dz_dv0;
	dPsiB1_dX_dY[iu][0][1] = local_PsiB1a*dz_dv0_dv1 + local_PsiB1b*dz_dv0*dz_dv1;
	dPsiB1_dX_dY[iu][0][2] = local_PsiB1a*dz_dv0_dv2 + local_PsiB1b*dz_dv0*dz_dv2;

	dPsiB1_dX_dY[iu][1][0] = local_PsiB1a*dz_dv1_dv0 + local_PsiB1b*dz_dv1*dz_dv0;
	dPsiB1_dX_dY[iu][1][1] = local_PsiB1a*dz_dv1_dv1 + local_PsiB1b*dz_dv1*dz_dv1;
	dPsiB1_dX_dY[iu][1][2] = local_PsiB1a*dz_dv1_dv2 + local_PsiB1b*dz_dv1*dz_dv2;

	dPsiB1_dX_dY[iu][2][0] = local_PsiB1a*dz_dv2_dv0 + local_PsiB1b*dz_dv2*dz_dv0;
	dPsiB1_dX_dY[iu][2][1] = local_PsiB1a*dz_dv2_dv1 + local_PsiB1b*dz_dv2*dz_dv1;
	dPsiB1_dX_dY[iu][2][2] = local_PsiB1a*dz_dv2_dv2 + local_PsiB1b*dz_dv2*dz_dv2;

	//second derivatives of Psi_B2
	dPsiB2_dX_dY[iu][0][0] = local_PsiB2a*dz_dv0_dv0 + local_PsiB2b*dz_dv0*dz_dv0;
	dPsiB2_dX_dY[iu][0][1] = local_PsiB2a*dz_dv0_dv1 + local_PsiB2b*dz_dv0*dz_dv1;
	dPsiB2_dX_dY[iu][0][2] = local_PsiB2a*dz_dv0_dv2 + local_PsiB2b*dz_dv0*dz_dv2;

	dPsiB2_dX_dY[iu][1][0] = local_PsiB2a*dz_dv1_dv0 + local_PsiB2b*dz_dv1*dz_dv0;
	dPsiB2_dX_dY[iu][1][1] = local_PsiB2a*dz_dv1_dv1 + local_PsiB2b*dz_dv1*dz_dv1;
	dPsiB2_dX_dY[iu][1][2] = local_PsiB2a*dz_dv1_dv2 + local_PsiB2b*dz_dv1*dz_dv2;

	dPsiB2_dX_dY[iu][2][0] = local_PsiB2a*dz_dv2_dv0 + local_PsiB2b*dz_dv2*dz_dv0;
	dPsiB2_dX_dY[iu][2][1] = local_PsiB2a*dz_dv2_dv1 + local_PsiB2b*dz_dv2*dz_dv1;
	dPsiB2_dX_dY[iu][2][2] = local_PsiB2a*dz_dv2_dv2 + local_PsiB2b*dz_dv2*dz_dv2;

	return;
}


//derivatives and values of full Psi functions

inline void set_Psi_k()
{
	for (int iu = 0; iu < n_u_pts; ++iu)
	{
		double u_loc = u_pts[iu];
		Psi0[iu] = Psi_0(u_loc);
		Psi1[iu] = Psi_1(u_loc);
		Psi2[iu] = Psi_2(u_loc);
	}
	return;
}

inline void set_dPsik_dX(int iu, int iX)
{
	dPsi0_dX[iX][iu] = dPsiA_dX[iu][iX]*PsiB0[iu] + dPsiB0_dX[iu][iX]*PsiA[iu];
	dPsi1_dX[iX][iu] = dPsiA_dX[iu][iX]*PsiB1[iu] + dPsiB1_dX[iu][iX]*PsiA[iu];
	dPsi2_dX[iX][iu] = dPsiA_dX[iu][iX]*PsiB2[iu] + dPsiB2_dX[iu][iX]*PsiA[iu];
}

inline void set_dPsik_dX_dY(int iu, int iX, int iY)
{
	dPsi0_dX_dY[iX][iY][iu] = 
			dPsiA_dX_dY[iu][iX][iY] * PsiB0[iu]
			+ dPsiA_dX[iu][iX] * dPsiB0_dX[iu][iY]
			+ dPsiB0_dX[iu][iX] * dPsiA_dX[iu][iY]
			+ dPsiB0_dX_dY[iu][iX][iY] * PsiA[iu]
			;

	dPsi1_dX_dY[iX][iY][iu] = 
			dPsiA_dX_dY[iu][iX][iY] * PsiB1[iu]
			+ dPsiA_dX[iu][iX] * dPsiB1_dX[iu][iY]
			+ dPsiB1_dX[iu][iX] * dPsiA_dX[iu][iY]
			+ dPsiB1_dX_dY[iu][iX][iY] * PsiA[iu]
			;

	dPsi2_dX_dY[iX][iY][iu] = 
			dPsiA_dX_dY[iu][iX][iY] * PsiB2[iu]
			+ dPsiA_dX[iu][iX] * dPsiB2_dX[iu][iY]
			+ dPsiB2_dX[iu][iX] * dPsiA_dX[iu][iY]
			+ dPsiB2_dX_dY[iu][iX][iY] * PsiA[iu]
			;
}

// End of file

#endif
