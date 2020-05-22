#include <math.h>
#include "InverseCdf.hpp"

void InverseCdf::apply_exponent()
{
	for (int i = 0; i <= CDF_RES; i++)
	{
		cdf_values_[i] = exp(-cdf_values_[i]);
	}
}

// set defaults if not present in parameter file
void InverseCdf::set_defaults(int startValue)
{
	cdf_values_[CDF_RES] = startValue;
	for (int i = 0; i < CDF_RES; i++)
		cdf_values_[i] = -log(1 - ((double)i) / CDF_RES);
}


