#pragma once

#include "Constants.h"

class InverseCdf
{
	double cdf_values_[CDF_RES + 1];

public:
	// default start value for icdf double arrays (was hardcoded as 100)
	static constexpr double start_value_ = 100.0;

	void set_neg_log(double start_value);

	void assign_exponent();

	void assign_exponent(double value);

	unsigned short int choose(double Mean, int tn, double timesteps_per_day);

	// Getter
	double* get_values()
	{
		return cdf_values_;
	}

	// Overloading [] operator to access elements in array style 
	double& operator[](int i)
	{
		return cdf_values_[i];
	}
};




