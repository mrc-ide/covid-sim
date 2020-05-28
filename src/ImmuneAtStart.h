#pragma once
#include "IInfectionState.h"

class ImmuneAtStart : public IInfectionState
{
public:
	virtual void GetsWorse(int ai, double t, int tn, int run);
	virtual void GetsBetter(int ai, double t, int tn, int run);
};

