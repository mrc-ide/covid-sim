#pragma once
#include "Param.h"
class IInfectionState
{
public:
	virtual void GetsWorse(int ai, double t, int tn, int run) = 0;
	virtual void GetsBetter(int ai, double t, int tn, int run) = 0;

	Param* P;
public:
	IInfectionState(Param* p)
	{
		P = p;
	}
	IInfectionState()
	{
			
	}

};
