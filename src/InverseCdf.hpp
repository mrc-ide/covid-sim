#pragma once

#include "Constants.h"

class InverseCdf
{
	double cdf_values_[CDF_RES + 1];

public:

	void set_defaults(double startValue);

	void apply_exponent();
	void apply_exponent(double value);

	// Getter
	double* get_values()
	{
		return cdf_values_;
	}

	double get_value(int i)
	{
		return cdf_values_[i];
	}
};




