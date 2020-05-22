#pragma once

#include "Constants.h"

class InverseCdf
{
	double cdf_values_[CDF_RES + 1];

public:

	void set_defaults(int startValue);

	void apply_exponent();

	// Getter
	double* get_values()
	{
		return cdf_values_;
	}
	

	
};




