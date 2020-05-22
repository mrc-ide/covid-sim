#include "InverseCdf.hpp"
#include <math.h>

void InverseCdf::apply_exponent()
{
	for (int i = 0; i <= CDF_RES; i++)
	{
		cdf_values_[i] = exp(-cdf_values_[i]);
	}
}

void InverseCdf::apply_exponent(double value)
{
	for (int i = 0; i <= CDF_RES; i++)
	{
		cdf_values_[i] = exp(value);
	}
}

// set defaults if not present in parameter file
void InverseCdf::set_defaults(int startValue)
{
	cdf_values_[CDF_RES] = startValue;
	for (int i = 0; i < CDF_RES; i++)
		cdf_values_[i] = -log(1 - ((double)i) / CDF_RES);
}


