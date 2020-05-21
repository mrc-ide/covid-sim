#pragma once


#include "Constants.h"

class ICDF
{
	double cdf_values_[CDF_RES + 1];

public:

	void SetDefaults(int startValue);

	void ApplyExponent();

	// Getter
	double* GetCDFValues()
	{
		return cdf_values_;
	}
	

	
};




