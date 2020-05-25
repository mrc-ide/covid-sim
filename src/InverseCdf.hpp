#pragma once

#include "Constants.h"

class InverseCdf
{
	double cdf_values_[CDF_RES + 1];

public:

	void set_defaults(double start_value);

	void apply_exponent();
	void apply_exponent(double value);

	// Getter
	double* get_values()
	{
		return cdf_values_;
	}

	// Overloading [] operator to access elements in array style 
	double& InverseCdf::operator[](int i);
};




