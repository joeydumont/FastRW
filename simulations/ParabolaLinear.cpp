#include <fastrw>
#include <cuba.h>
#include <boost/math/constants/constants.hpp>
#include <armadillo>

using namespace FastRW;

const double pi             = boost::math::constants::pi<double>();
const double UNIT_LENGTH    = 3.86159e-13;
const double UNIT_TIME      = 1.2880885e-21;
const double SPEED_OF_LIGHT = 299792458;
const double EPSILON_0      = 8.85418782e-12;
const double MU_0           = 4.0*pi*1.0e-7;
const double UNIT_E_FIELD   = 1.3e18;
const double UNIT_B_FIELD   = UNIT_E_FIELD/SPEED_OF_LIGHT;
const double UNIT_MASS      = 9.109382914e-31;

const double f           = 0.04375/UNIT_LENGTH;
const double w0          = 0.075/UNIT_LENGTH;
const double lamb        = 800.0e-9/UNIT_LENGTH;
const double k           = 2.0*pi/lamb;
const double alpha_max   = pi/2.0;

std::complex<double> Ex(double alpha, double beta)
{
	return std::exp(-std::pow(f*std::sin(alpha)/w0,2));
}

std::complex<double> Ey(double alpha, double beta)
{
	return 0.0;
}

double q(double alpha)
{
	return std::pow(std::cos(alpha/2.0),-2);
}

int main(int argc, char* argv[])
{
	RichardsWolf hna_parabola_lin = RichardsWolf(&Ex,&Ey,&q,alpha_max,k);

	auto r = arma::linspace(0.0,2.5e-6/UNIT_LENGTH, 100);
	auto th= arma::linspace(0.0,2.0*pi,100);
	arma::cx_cube field(100,100,6);
	for (int i=0; i<100; i++)
	{
		for (int j=0; j<100; j++)
		{
			auto field_values = hna_parabola_lin.ComputeFieldComponents(r[i],th[j],0.0);
			for (int k=0; k<6; k++)
				field.tube(i,j)[k] = field_values[k];
		}
	}

	for (int i=0; i<6; i++)
	{
		std::stringstream file_name;
		file_name << "field-component-" << i << ".txt";
		arma::mat field_abs = arma::abs(field.slice(i));
		field_abs.save(file_name.str().c_str(), arma::arma_ascii);
	}

	return 0;
}