#ifndef FASTRW_HPP
#define FASTRW_HPP

#include <complex>
#include <cuba.h>
#include <boost/math/constants/constants.hpp>

using namespace std::complex_literals;

namespace FastRW {
namespace cst = boost::math::constants;

/// Enum to choose the quadrature method.
enum IntegrationMethod {CubaCuhre};

int interface_to_cuba_real(const int *ndim, const double x[],
	                       const int *fdim, double fval[], void *fdata);

int interface_to_cuba_imag(const int *ndim, const double x[],
	                       const int *fdim, double fval[], void *fdata);

typedef std::complex<double> (*rw_field)(double, double);
typedef double               (*rw_apodization)(double);

/*!
 *  \class RichardsWolf
 *  \author Joey Dumont			<joey.dumont@gmail.com>
 *  \since 2018-03-08
 *  \brief Implements the Richards-Wolf formalism for a single frequency.
 */
class RichardsWolf
{
public:
	RichardsWolf(rw_field       my_Ex,
		         rw_field       my_Ey,
		         rw_apodization my_q,
		         double         my_alpha_max,
		         double         my_k)
	: Ex_f(my_Ex)
	, Ey_f(my_Ey)
	, q_f(my_q)
	, alpha_max(my_alpha_max)
	, k(my_k)
	{
		xmin.resize(2);
		xmin[0] = 0.0;
		xmin[1] = 0.0;

		xmax.resize(2);
		xmax[0] = alpha_max;
		xmax[1] = 2.0*cst::pi<double>();

		abstol  = 1.0e-4;
		reltol  = 1.0e-4;
		mineval = 1;
		maxeval = 1e6;
	}

	std::vector<std::complex<double>> ComputeFieldComponents(double r_in, double th_in, double z_in)
	{
		std::vector<std::complex<double>> field_values(6);

		r  = r_in;
		th = th_in;
		z  = z_in;

		const int ndim = 2;
		const int fdim = 6;
		int nregions, neval, error_flag;
		double val[fdim], err[fdim], prob[fdim];

		Cuhre(ndim,
			  fdim,
			  interface_to_cuba_real,
			  this,
			  1,
			  reltol,abstol,
			  0,
			  mineval,maxeval,
			  0,
			  NULL,
			  NULL,
			  &nregions,&neval,&error_flag,
			  &val[0],&err[0],&prob[0]);

		for (int i=0; i<fdim;i++)
			field_values[i] = val[i];

		Cuhre(ndim,
			  fdim,
			  interface_to_cuba_imag,
			  this,
			  1,
			  reltol,abstol,
			  0,
			  mineval,maxeval,
			  0,
			  NULL,
			  NULL,
			  &nregions,&neval,&error_flag,
			  &val[0],&err[0],&prob[0]);

		for (int i=0; i<fdim;i++)
			field_values[i] += 1i*val[i];

		return field_values;
	}

protected:
	std::vector<std::complex<double>> RW_Integrand(double alpha, double beta)
	{
		// Evaluate some trig constants.
		double cosAlpha  = std::cos(alpha);
		double sinAlpha  = std::sin(alpha);
		double cosBeta   = std::cos(beta);
		double sinBeta   = std::sin(beta);
		double cosBetaSq = std::pow(cosBeta,2);
		double sinBetaSq = std::pow(sinBeta,2);

		// Compute the incident field.
		std::complex<double> Ex = Ex_f(alpha,beta);
		std::complex<double> Ey = Ey_f(alpha,beta);

		std::vector<std::complex<double>> field_values(6);
		std::complex<double>              phase = std::exp(1i*k*z*cosAlpha)*std::exp(1i*k*r*sinAlpha*std::cos(th-beta));

		field_values[0] = q_f(alpha)*(Ex*(cosAlpha*cosBetaSq+sinBetaSq)+Ey*cosBeta*sinBeta*(cosAlpha-1));
		field_values[1] = q_f(alpha)*(Ex*cosBeta*sinBeta*(cosAlpha-1.0)+Ey*(cosAlpha*sinBetaSq+cosBetaSq));
		field_values[2] = q_f(alpha)*(Ex*sinAlpha*cosBeta+Ey*sinAlpha*sinBeta);
		field_values[3] = q_f(alpha)*(Ex*cosBeta*sinBeta*(cosAlpha-1.0)+Ey*(sinBetaSq-cosAlpha*cosBetaSq));
		field_values[4] = q_f(alpha)*(Ex*(cosAlpha*sinBetaSq+cosBetaSq)+Ey*cosBeta*sinBeta*(cosAlpha+1.0));
		field_values[5] = q_f(alpha)*(Ex*sinAlpha*sinBeta-Ey*sinAlpha*cosBeta);

		for (uint i=0; i<6; i++)
			field_values[i] *= phase*sinAlpha;

		return field_values;
	}

	friend int interface_to_cuba_real(const int*, const double[], const int*, double[], void *);
	friend int interface_to_cuba_imag(const int*, const double[], const int*, double[], void *);

	IntegrationMethod method;
	std::vector<double> xmin;
	std::vector<double> xmax;
	double abstol;
	double reltol;
	int mineval;
	int maxeval;

	rw_field Ex_f;
	rw_field Ey_f;
	rw_apodization q_f;
	double alpha_max;
	double k;

	// Coordinates
	double r,th,z;
};

int interface_to_cuba_real(const int *ndim, const double x[],
	                       const int *fdim,       double fval[], void *fdata)
{
	auto obj = (RichardsWolf* ) fdata;

	// We scale the integrand by applying the transformation
	// 	int_a^b f(x)dx --> int_0^1 f[a+(b-a)*y]*(b-a)dy
	// in each dimension.
	double xscaled[3];
	double jacobian = 1.0;
	for (int i=0; i<(*ndim); i++)
	{
		double range =  obj->xmax[i] - obj->xmin[i];
		jacobian     *= range;
		xscaled[i]   =  obj->xmin[i]+range*x[i];
	}

	// We compute the Richards-Wolf integrands.
	auto field = obj->RW_Integrand(xscaled[0],xscaled[1]);

	// We store the field in the proper array.
	for (int i=0; i<(*fdim); i++)
	{
		fval[i] =  std::real(field[i]);
		fval[i] *= jacobian;
	}

	return 0;
}

int interface_to_cuba_imag(const int *ndim, const double x[],
	                       const int *fdim,       double fval[], void *fdata)
{
	auto obj = (RichardsWolf* ) fdata;

	// We scale the integrand by applying the transformation
	// 	int_a^b f(x)dx --> int_0^1 f[a+(b-a)*y]*(b-a)dy
	// in each dimension.
	double xscaled[3];
	double jacobian = 1.0;
	for (int i=0; i<(*ndim); i++)
	{
		double range =  obj->xmax[i] - obj->xmin[i];
		jacobian     *= range;
		xscaled[i]   =  obj->xmin[i]+range*x[i];
	}

	// We compute the Richards-Wolf integrands.
	auto field = obj->RW_Integrand(xscaled[0],xscaled[1]);

	// We store the field in the proper array.
	for (int i=0; i<(*fdim); i++)
	{
		fval[i] =  std::imag(field[i]);
		fval[i] *= jacobian;
	}

	return 0;
}


} // namespace FastRW

#endif // FASTRW_HPP