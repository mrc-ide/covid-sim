#include "InverseCdf.h"
#include <math.h>
#include "Rand.h"


// set defaults if not present in parameter file
void InverseCdf::set_defaults(double start_value)
{
	cdf_values_[CDF_RES] = start_value;
	for (int i = 0; i < CDF_RES; i++)
		cdf_values_[i] = -log(1 - ((double)i) / CDF_RES);
}

void InverseCdf::apply_exponent()
{
	for (int quantile = 0; quantile <= CDF_RES; quantile++)
	{
		cdf_values_[quantile] = exp(-cdf_values_[quantile]);
	}
}

void InverseCdf::apply_exponent(double value)
{
	for (int i = 0; i <= CDF_RES; i++)
	{
		cdf_values_[i] = exp(value);
	}
}

unsigned short int InverseCdf::choose(double Mean, int tn, double timesteps_per_day)
{
	unsigned short int Value;
	int i;
	double q, ti;

	i = (int)floor(q = ranf_mt(tn) * CDF_RES); //// note q defined here as well as i.
	q -= ((double)i); //// remainder

	//// weighted average (sort of) between quartile values from CDF_RES. logged as it was previously exponentiated in ReadParams. Minus as exp(-cdf) was done in ReadParaams. Sort of
	ti = -Mean * log(q * cdf_values_[i + 1] + (1.0 - q) * cdf_values_[i]); 
	Value = (unsigned short int) floor(0.5 + (ti * timesteps_per_day));

	return Value;
}



